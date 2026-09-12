/* Layer 4: the library contract, from a consumer's seat.
 *
 * methscope links libyame.a and calls these directly, so a change in what
 * they RETURN is as breaking as a changed signature and nothing in the CLI
 * would notice. Both shipped breakages this month were of that kind: the
 * bounded read that did not exist, and free_cdata() starting to zero n, which
 * made methscope report 0 CpGs from a record it had already freed.
 *
 * Exits 0 on success; on failure prints which assertion and exits 1.
 */
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include "cfile.h"
#include <unistd.h>
#include <stdlib.h>

static int fails = 0;
#define CHECK(cond, ...) do {                                   \
    if (!(cond)) { fprintf(stderr, "  FAIL %s:%d: ", __func__, __LINE__); \
                   fprintf(stderr, __VA_ARGS__);                \
                   fputc('\n', stderr); fails++; }              \
  } while (0)

/* A record read back is the record that was written. */
static void t_read(const char *path, uint64_t want_n, char want_fmt) {
  cfile_t cf = open_cfile((char *) path);
  cdata_t c = {0};
  char err[CX_ERRBUF];
  cx_read_t st = cx_read_record(&cf, &c, CX_NO_LIMIT, err, sizeof err);
  CHECK(st == CX_READ_OK, "first read of %s returned %d (%s)", path, (int) st, err);
  CHECK(c.fmt == want_fmt, "fmt is '%c', want '%c'", c.fmt, want_fmt);

  /* A record arrives COMPRESSED, and while it is, n counts BYTES, not rows.
   * It becomes a row count only after decompress_in_situ(). Reading n as a
   * row count off a freshly read record is a mistake a consumer makes once. */
  CHECK(c.compressed == 1, "a freshly read record is not marked compressed");
  CHECK(c.n > 0, "a freshly read record has n == 0");
  decompress_in_situ(&c);
  CHECK(c.compressed == 0, "decompress_in_situ left the record marked compressed");
  CHECK(c.n == want_n, "decompressed n is %" PRIu64 ", want %" PRIu64, c.n, want_n);

  /* End of stream is CX_READ_END, not an error and not another record. */
  st = cx_read_record(&cf, &c, CX_NO_LIMIT, err, sizeof err);
  CHECK(st == CX_READ_END, "past the last record returned %d (%s)", (int) st, err);

  /* free_cdata() zeroes n and nulls the pointers, and is idempotent. A
   * consumer that reads a field after freeing gets 0, not a stale value --
   * that is a contract, not an accident, and calling it twice is safe. */
  free_cdata(&c);
  CHECK(c.n == 0, "free_cdata left n at %" PRIu64, c.n);
  CHECK(c.s == NULL, "free_cdata left s non-NULL");
  CHECK(c.aux == NULL, "free_cdata left aux non-NULL");
  free_cdata(&c);               /* must not double free */
  bgzf_close(cf.fh);
}

/* A CX stream that is a PREFIX of a larger file: bounded, it ends cleanly at
 * the limit; unbounded, it reports NOT_CX and does NOT exit the process. */
static void t_bounded(const char *path, int64_t limit) {
  cfile_t cf = open_cfile((char *) path);
  cdata_t c = {0};
  char err[CX_ERRBUF];
  int n = 0;
  for (;;) {
    cx_read_t st = cx_read_record(&cf, &c, limit, err, sizeof err);
    if (st == CX_READ_END) break;
    CHECK(st == CX_READ_OK, "bounded read returned %d (%s)", (int) st, err);
    if (st != CX_READ_OK) break;
    n++;
  }
  CHECK(n == 1, "bounded read saw %d records, want 1", n);
  free_cdata(&c);
  bgzf_close(cf.fh);

  cf = open_cfile((char *) path);
  memset(&c, 0, sizeof c);
  cx_read_t st = cx_read_record(&cf, &c, CX_NO_LIMIT, err, sizeof err);
  CHECK(st == CX_READ_OK, "unbounded first read returned %d", (int) st);
  st = cx_read_record(&cf, &c, CX_NO_LIMIT, err, sizeof err);
  CHECK(st == CX_READ_NOT_CX, "unbounded past the prefix returned %d, want NOT_CX", (int) st);
  CHECK(err[0] != '\0', "a failure status came back with no message");
  free_cdata(&c);
  bgzf_close(cf.fh);
}

/* read_cdata1() is the older entry point methscope still uses: a zeroed
 * record means end of stream. */
static void t_read_cdata1(const char *path, int want_records) {
  cfile_t cf = open_cfile((char *) path);
  int n = 0;
  for (;;) {
    cdata_t c = read_cdata1(&cf);
    if (c.n == 0) { free_cdata(&c); break; }
    n++;
    free_cdata(&c);
  }
  CHECK(n == want_records, "read_cdata1 walked %d records, want %d", n, want_records);
  bgzf_close(cf.fh);
}

/* The accessors a consumer does arithmetic with. */
static void t_accessors(const char *path) {
  cfile_t cf = open_cfile((char *) path);
  cdata_t c = read_cdata1(&cf);
  decompress_in_situ(&c);
  CHECK(c.fmt == '3', "fixture is not format 3");
  uint64_t mu = f3_get_mu(&c, 0);            /* fixture row 0 is 3 1 */
  CHECK((mu >> 32) == 3, "M is %" PRIu64 ", want 3", mu >> 32);
  CHECK((mu & 0xffffffff) == 1, "U is %" PRIu64 ", want 1", mu & 0xffffffff);
  CHECK(MU2cov(mu) == 4, "coverage is %" PRIu64 ", want 4", (uint64_t) MU2cov(mu));
  double b = MU2beta(mu);
  CHECK(b > 0.749 && b < 0.751, "beta is %f, want 0.75", b);
  mu = f3_get_mu(&c, 1);                     /* row 1 is 0 0: uncovered */
  CHECK(MU2cov(mu) == 0, "an uncovered site reports coverage %" PRIu64, (uint64_t) MU2cov(mu));
  free_cdata(&c);
  bgzf_close(cf.fh);
}

/* ---- the store layer: pure functions a fetch depends on ------------------- */
#include "assets.h"
#include <stdlib.h>
#include <sys/stat.h>
#include <stdio.h>

/* SHA-256 against the published vectors: the whole download path trusts
 * this one function, and a wrong digest either rejects every file or accepts
 * any. */
static void t_sha256(void) {
  char out[65];
  yame_assets_sha256_buf("", 0, out);
  CHECK(strcmp(out, "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855") == 0,
        "sha256(\"\") is %s", out);
  yame_assets_sha256_buf("abc", 3, out);
  CHECK(strcmp(out, "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad") == 0,
        "sha256(\"abc\") is %s", out);
  /* a string that crosses one 64-byte block boundary */
  const char *two = "abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq";
  yame_assets_sha256_buf(two, strlen(two), out);
  CHECK(strcmp(out, "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1") == 0,
        "two-block sha256 is %s", out);
  CHECK(yame_assets_digest_equal(out, "248D6A61D20638B8E5C026930C3E6039A33CE45964FF2167F6ECEDD419DB06C1"),
        "digest_equal is case-sensitive");
  CHECK(!yame_assets_digest_equal(out, "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c2"),
        "digest_equal accepted a one-nibble difference");
}

/* The manifest parser: what it keeps, what it drops, and that a path trying
 * to escape the store directory is dropped rather than followed. */
static void t_parse_sums(void) {
  size_t n = 0;
  const char *text =
    "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855  a.cm\n"
    "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad  sub/b.cm\n"
    "\n"
    "# a comment line, or at least a line that is not an entry\n"
    "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1  ../escape.cm\n"
    "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1  /abs/olute.cm\n";
  yame_sums_ent_t *e = yame_assets_parse_sums(text, &n);
  CHECK(e != NULL, "parse_sums returned NULL");
  CHECK(n == 2, "parse_sums kept %zu entries, want 2 (the two safe ones)", n);
  if (e && n >= 2) {
    CHECK(strcmp(e[0].name, "a.cm") == 0, "entry 0 name is %s", e[0].name);
    CHECK(strcmp(e[1].name, "sub/b.cm") == 0, "entry 1 name is %s", e[1].name);
    CHECK(strcmp(e[1].sha, "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad") == 0,
          "entry 1 sha is %s", e[1].sha);
  }
  free(e);
  n = 99;
  e = yame_assets_parse_sums("", &n);
  CHECK(n == 0, "an empty manifest parsed to %zu entries", n);
  free(e);
  CHECK(!yame_assets_safe_relpath("../x"), "safe_relpath accepted ../x");
  CHECK(!yame_assets_safe_relpath("/x"), "safe_relpath accepted /x");
  CHECK(yame_assets_safe_relpath("a/b.cm"), "safe_relpath rejected a/b.cm");
  CHECK(!yame_assets_safe_name("a/b"), "safe_name accepted a slash");
  CHECK(yame_assets_safe_name("EPIC.hg38.mask.cm"), "safe_name rejected a real name");
}

/* The pin: a stored manifest that hashes to the anchor is MATCH, a missing one
 * is ABSENT, a stranger's is CONFLICT, and a known earlier tag's is ANCESTOR.
 * The conflict is what stops two tools re-downloading over each other forever;
 * the ancestor is what lets a build upgrade its own store without -f. */
static void t_pin(const char *dir) {
  char sums[4096];
  snprintf(sums, sizeof sums, "%s/SHA256SUMS", dir);
  unlink(sums);
  CHECK(yame_assets_pin_check(dir, "0000") == YAME_PIN_ABSENT, "no manifest is not ABSENT");

  const char *text = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855  a.cm\n";
  FILE *f = fopen(sums, "w"); fputs(text, f); fclose(f);
  char anchor[65];
  yame_assets_sha256_buf(text, strlen(text), anchor);
  CHECK(yame_assets_pin_check(dir, anchor) == YAME_PIN_MATCH, "own manifest is not MATCH");
  CHECK(yame_assets_pin_check(dir, NULL) == YAME_PIN_UNKNOWN, "no anchor is not UNKNOWN");
  CHECK(yame_assets_pin_check(dir, "ffff") == YAME_PIN_CONFLICT, "a stranger's manifest is not CONFLICT");

  yame_pin_prior_t prior = { "v0.9", anchor };
  CHECK(yame_assets_pin_state(dir, "ffff", &prior, 1) == YAME_PIN_ANCESTOR,
        "a known earlier tag is not ANCESTOR");
  const char *t = yame_assets_pin_prior_tag(dir, &prior, 1);
  CHECK(t && strcmp(t, "v0.9") == 0, "prior tag is %s, want v0.9", t ? t : "(null)");
  CHECK(yame_assets_pin_state(dir, "ffff", NULL, 0) == YAME_PIN_CONFLICT,
        "with no priors, pin_state is not plain pin_check");
  CHECK(yame_assets_pin_prior_tag(dir, NULL, 0) == NULL, "prior tag without priors is not NULL");
  unlink(sums);
}

/* Reference resolution: a row count identifies the row space outright
 * (29,401,795 rows is hg38 and nothing else), and a -R/-m argument is either
 * a path used as given or a name looked up in the store the row space owns.
 * It resolves only; it never downloads. */
static void t_refstore(const char *store) {
  char path[4096]; const char *name = NULL, *fetch = NULL;

  /* hg38 by row count, with nothing fetched: identified, missing, and the
   * message would know which fetch supplies it */
  int st = yame_ref_for_rows(29401795, store, "genome", path, sizeof path, &name, &fetch);
  CHECK(st == YAME_REF_MISSING, "hg38 by rows in an empty store returned %d, want MISSING", st);
  CHECK(name && strcmp(name, "hg38") == 0, "row space name is %s, want hg38", name ? name : "(null)");
  CHECK(fetch && *fetch, "no fetch hint for a missing reference");

  /* a count matching nothing */
  st = yame_ref_for_rows(12345, store, NULL, path, sizeof path, &name, &fetch);
  CHECK(st == YAME_REF_UNKNOWN, "an unknown row count returned %d, want UNKNOWN", st);

  /* the wrong kind: hg38 is a genome, not an array */
  st = yame_ref_for_rows(29401795, store, "array", path, sizeof path, &name, &fetch);
  CHECK(st == YAME_REF_WRONG_KIND || st == YAME_REF_MISSING,
        "hg38 asked for as an array returned %d", st);

  /* an existing path is used as given, whatever the row count says */
  char own[4096]; snprintf(own, sizeof own, "%s/own.cr", store);
  FILE *f = fopen(own, "w"); fputs("x", f); fclose(f);
  st = yame_ref_resolve(own, 29401795, store, NULL, path, sizeof path, &name, &fetch);
  CHECK(st == YAME_REF_OK, "an existing path returned %d, want OK", st);
  CHECK(strcmp(path, own) == 0, "an existing path was rewritten to %s", path);

  /* a name in the row space's own directory resolves once the file exists */
  char dir[4096]; snprintf(dir, sizeof dir, "%s/hg38", store);
  mkdir(dir, 0755);
  char cr[4096]; snprintf(cr, sizeof cr, "%s/cpg_nocontig.cr", dir);
  f = fopen(cr, "w"); fputs("x", f); fclose(f);
  st = yame_ref_for_rows(29401795, store, "genome", path, sizeof path, &name, &fetch);
  CHECK(st == YAME_REF_OK, "hg38 by rows with its .cr present returned %d, want OK", st);
  CHECK(strcmp(path, cr) == 0, "hg38 resolved to %s, want %s", path, cr);
  st = yame_ref_resolve("hg38", 29401795, store, "genome", path, sizeof path, &name, &fetch);
  CHECK(st == YAME_REF_OK, "the name hg38 returned %d, want OK", st);

  /* a name nothing in that row space carries */
  st = yame_ref_resolve("NoSuchMask", 29401795, store, NULL, path, sizeof path, &name, &fetch);
  CHECK(st == YAME_REF_NO_NAME || st == YAME_REF_MISSING,
        "an unknown name returned %d, want NO_NAME", st);

  /* resolve_multi: one name can stand for SEVERAL files -- a knowledgebase
   * directory's worth -- and the caller frees the array. */
  char **paths = NULL; size_t np = 0;
  st = yame_ref_resolve_multi("hg38", 29401795, store, "genome", &paths, &np, &name, &fetch);
  if (st == YAME_REF_OK) {
    CHECK(np >= 1, "resolve_multi returned OK with %zu paths", np);
    CHECK(paths != NULL, "resolve_multi returned OK with a NULL array");
    for (size_t i = 0; i < np; ++i)
      CHECK(paths[i] && paths[i][0] == '/', "path %zu is not absolute: %s",
            i, paths[i] ? paths[i] : "(null)");
    yame_ref_paths_free(paths, np);
  } else {
    CHECK(np == 0, "resolve_multi failed but reported %zu paths", np);
    yame_ref_paths_free(paths, np);      /* must tolerate the empty case */
  }
  paths = NULL; np = 0;
  st = yame_ref_resolve_multi("NoSuchName", 29401795, store, NULL, &paths, &np, &name, &fetch);
  CHECK(st != YAME_REF_OK, "resolve_multi accepted a name nothing carries");
  yame_ref_paths_free(paths, np);

  /* The explain functions write the message a command prints when a reference
   * is missing. They must say something for every status, and name the flag
   * the caller was using. */
  FILE *devnull = fopen("/dev/null", "w");
  if (devnull) {
    int codes[] = { YAME_REF_OK, YAME_REF_MISSING, YAME_REF_WRONG_KIND,
                    YAME_REF_NO_NAME, YAME_REF_UNKNOWN };
    for (size_t i = 0; i < sizeof codes / sizeof *codes; ++i) {
      yame_ref_explain(devnull, 29401795, codes[i], "hg38", "hg38/cpg_nocontig.cr", "-R");
      yame_ref_explain_name(devnull, "ChromHMM", 29401795, codes[i], "hg38",
                            "hg38/KYCG", "-m");
    }
    fclose(devnull);
  }
  /* and a real message goes somewhere a test can read */
  char msg[4096]; FILE *mf = fmemopen(msg, sizeof msg, "w");
  if (mf) {
    yame_ref_explain(mf, 12345, YAME_REF_UNKNOWN, NULL, NULL, "-R");
    fclose(mf);
    CHECK(msg[0] != '\0', "explain printed nothing for an unknown row count");
  }

  /* the row count of a real file */
  unlink(own); unlink(cr);
}

/* The store-state helper, against a registry the PROBE defines -- which is
 * the point: a downstream tool passes its own, and gets advice spelled with
 * its own verb. Every state is reached by writing a manifest whose digest is
 * (or is not) one this registry knows. */
static void t_store_state(const char *store) {
  /* a manifest text, and its digest, which the registry will call current */
  const char *cur_text = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855  a.cm\n";
  char cur_anchor[65]; yame_assets_sha256_buf(cur_text, strlen(cur_text), cur_anchor);
  const char *old_text = "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad  a.cm\n";
  char old_anchor[65]; yame_assets_sha256_buf(old_text, strlen(old_text), old_anchor);

  static yame_asset_file_t files[] = {
    { "a.cm", "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855", 0, NULL },
    { NULL, NULL, 0, NULL } };
  yame_pin_prior_t prior[] = { { "v1", old_anchor } };
  yame_asset_reg_t reg[] = {
    { "probe", "hg38/probe", "http://x", "v2", "", "hg38/probe", cur_anchor,
      files, YAME_NFILES(files), prior, 1 } };
  const size_t n_reg = 1;

  char dir[4096]; snprintf(dir, sizeof dir, "%s/hg38/probe", store);
  char sums[4096]; snprintf(sums, sizeof sums, "%s/SHA256SUMS", dir);
  char afile[4096]; snprintf(afile, sizeof afile, "%s/a.cm", dir);
  mkdir(dir, 0755);
  char adv[1024];
  yame_store_state_t st;

  /* the helper resolves the root from tool_env; point it at the probe's store */
  setenv("PROBE_DATA_HOME", store, 1);

  unlink(sums); unlink(afile);
  st = yame_store_state(reg, n_reg, "probe", "PROBE_DATA_HOME", dir, adv, sizeof adv);
  CHECK(st == YAME_STORE_ABSENT, "empty directory is %d, want ABSENT", (int) st);
  CHECK(strstr(adv, "probe fetch hg38/probe") != NULL, "ABSENT advice is: %s", adv);

  FILE *f = fopen(sums, "w"); fputs(old_text, f); fclose(f);
  st = yame_store_state(reg, n_reg, "probe", "PROBE_DATA_HOME", dir, adv, sizeof adv);
  CHECK(st == YAME_STORE_OLD_TAG, "a prior tag's manifest is %d, want OLD_TAG", (int) st);
  CHECK(strstr(adv, "earlier tag") && strstr(adv, "probe fetch hg38/probe"),
        "OLD_TAG advice is: %s", adv);

  f = fopen(sums, "w"); fputs("ffff  z.cm\n", f); fclose(f);
  st = yame_store_state(reg, n_reg, "probe", "PROBE_DATA_HOME", dir, adv, sizeof adv);
  CHECK(st == YAME_STORE_OTHER_TAG, "an unknown manifest is %d, want OTHER_TAG", (int) st);
  CHECK(strstr(adv, "update probe") != NULL, "OTHER_TAG advice does not say to update the tool: %s", adv);

  f = fopen(sums, "w"); fputs(cur_text, f); fclose(f);
  st = yame_store_state(reg, n_reg, "probe", "PROBE_DATA_HOME", dir, adv, sizeof adv);
  CHECK(st == YAME_STORE_CURRENT, "the current manifest is %d, want CURRENT", (int) st);
  CHECK(adv[0] == '\0', "CURRENT left advice: %s", adv);

  /* now a FILE in a current directory: missing, wrong, right */
  st = yame_store_state(reg, n_reg, "probe", "PROBE_DATA_HOME", afile, adv, sizeof adv);
  CHECK(st == YAME_STORE_MISSING_FILE, "a listed file not on disk is %d, want MISSING_FILE", (int) st);
  f = fopen(afile, "w"); fputs("not empty", f); fclose(f);
  st = yame_store_state(reg, n_reg, "probe", "PROBE_DATA_HOME", afile, adv, sizeof adv);
  CHECK(st == YAME_STORE_STALE_FILE, "a file with the wrong digest is %d, want STALE_FILE", (int) st);
  f = fopen(afile, "w"); fclose(f);          /* the empty file: its sha256 is e3b0...  */
  st = yame_store_state(reg, n_reg, "probe", "PROBE_DATA_HOME", afile, adv, sizeof adv);
  CHECK(st == YAME_STORE_CURRENT, "a file with the right digest is %d, want CURRENT", (int) st);

  /* a path the registry knows nothing about */
  st = yame_store_state(reg, n_reg, "probe", "PROBE_DATA_HOME", "/nowhere/else", adv, sizeof adv);
  CHECK(st == YAME_STORE_NOT_CATALOGUED, "an uncatalogued path is %d", (int) st);

  /* A file in a CURRENT directory that this registry does not list: what a
   * model withdrawn upstream looks like to a store that still holds it. It
   * used to read as CURRENT, so a downstream load-time check called a
   * withdrawn file healthy. The directory itself stays current. */
  char gone[4096]; snprintf(gone, sizeof gone, "%s/withdrawn.clfx", dir);
  f = fopen(gone, "w"); fputs("bytes from an earlier tag", f); fclose(f);
  st = yame_store_state(reg, n_reg, "probe", "PROBE_DATA_HOME", gone, adv, sizeof adv);
  CHECK(st == YAME_STORE_NOT_LISTED, "an unlisted file is %d, want NOT_LISTED", (int) st);
  CHECK(strstr(adv, "withdrawn.clfx") && strstr(adv, "does not list it"),
        "NOT_LISTED advice is: %s", adv);
  CHECK(strstr(adv, "nothing will remove it") != NULL,
        "NOT_LISTED advice does not say the file is left alone: %s", adv);
  /* and it is left alone, by this call and by anything else */
  CHECK(yame_assets_is_file(gone), "the unlisted file was removed");
  st = yame_store_state(reg, n_reg, "probe", "PROBE_DATA_HOME", dir, adv, sizeof adv);
  CHECK(st == YAME_STORE_CURRENT, "the directory holding it is %d, want CURRENT", (int) st);
  unlink(gone);

  /* the report: one line for the old-tag case, none when current */
  f = fopen(sums, "w"); fputs(old_text, f); fclose(f);
  char buf[2048]; FILE *mf = fmemopen(buf, sizeof buf, "w");
  int said = yame_store_report(reg, n_reg, "probe", "PROBE_DATA_HOME", NULL, mf); fclose(mf);
  CHECK(said == 1, "report printed %d lines for one old-tag directory, want 1", said);
  CHECK(strstr(buf, "[probe fetch]") != NULL, "report line is not prefixed with the tool: %s", buf);
  f = fopen(sums, "w"); fputs(cur_text, f); fclose(f);
  mf = fmemopen(buf, sizeof buf, "w");
  said = yame_store_report(reg, n_reg, "probe", "PROBE_DATA_HOME", NULL, mf); fclose(mf);
  CHECK(said == 0, "report printed %d lines for a current store, want 0", said);

  unlink(sums); unlink(afile); rmdir(dir);
}

int main(int argc, char **argv) {
  if (argc < 4) { fprintf(stderr, "usage: probe <one.cg> <three.cg> <bundle> <limit> [dir]\n"); return 2; }
  if (argc > 5) { t_refstore(argv[5]); t_store_state(argv[5]); }
  t_read(argv[1], 4, '3');
  t_read_cdata1(argv[2], 3);
  t_accessors(argv[1]);
  t_sha256();
  {   /* the directory the running binary sits in -- used to find a store
       * beside an installed yame. It must produce an absolute path or say it
       * could not, never a half-filled buffer. */
    char exe[4096]; exe[0] = '\0';
    int ok = yame_assets_exe_dir(exe, sizeof exe);
    CHECK(ok == 0 || exe[0] == '/', "exe_dir returned %d with %s", ok, exe);
    char tiny[2]; tiny[0] = 'x';
    yame_assets_exe_dir(tiny, sizeof tiny);   /* must not overrun */
  }
  t_parse_sums();
  if (argc > 5) t_pin(argv[5]);
  if (argc > 4) t_bounded(argv[3], (int64_t) strtoll(argv[4], NULL, 10));
  if (fails) { fprintf(stderr, "%d library assertion(s) failed\n", fails); return 1; }
  return 0;
}
