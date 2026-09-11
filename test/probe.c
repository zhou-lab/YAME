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

int main(int argc, char **argv) {
  if (argc < 4) { fprintf(stderr, "usage: probe <one.cg> <three.cg> <bundle> <limit>\n"); return 2; }
  t_read(argv[1], 4, '3');
  t_read_cdata1(argv[2], 3);
  t_accessors(argv[1]);
  if (argc > 4) t_bounded(argv[3], (int64_t) strtoll(argv[4], NULL, 10));
  if (fails) { fprintf(stderr, "%d library assertion(s) failed\n", fails); return 1; }
  return 0;
}
