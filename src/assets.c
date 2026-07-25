// SPDX-License-Identifier: AGPL-3.0-or-later
/**
 * This file is part of YAME.
 *
 * Copyright (C) 2021-present Wanding Zhou
 *
 * YAME is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * YAME is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with YAME.  If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * The shared asset store: resolution, verification, and transfer.
 *
 * Ported from kycg's src/store.c and the networking half of its src/fetch.c,
 * generalized so sesame-cli and methscope-cli can drop their own copies. See
 * assets.h for the store layout and the trust chain.
 */

#include "assets.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>
#include <dirent.h>
#include <fcntl.h>
#include <time.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

#ifdef YAME_HAVE_CURL
#include <curl/curl.h>
#include "yame_version.h"
#endif

/* A store path. Deliberately generous: the longest real one is a nested
 * manifest entry under a platform dir under the root. */
#define YAME_PATH_MAX 4096

static void set_err(char **err, const char *fmt, ...) {
  if (!err) return;
  char buf[1024];
  va_list ap;
  va_start(ap, fmt);
  vsnprintf(buf, sizeof(buf), fmt, ap);
  va_end(ap);
  free(*err);
  *err = strdup(buf);
}

/* ------------------------------------------------------------- store paths */

int yame_assets_safe_name(const char *s) {
  if (!s || !*s) return 0;
  if (s[0] == '.') return 0;               /* ".", "..", hidden files */
  for (const char *p = s; *p; ++p) {
    if (*p == '/' || *p == '\\') return 0;
    if ((unsigned char)*p < 0x20) return 0;   /* control chars, incl. newline */
  }
  return 1;
}

int yame_assets_safe_relpath(const char *s) {
  if (!s || !*s) return 0;
  if (*s == '/') return 0;                 /* never absolute */

  char tmp[YAME_PATH_MAX];
  if (snprintf(tmp, sizeof(tmp), "%s", s) >= (int)sizeof(tmp)) return 0;

  /* Every component must stand on its own: one ".." anywhere in the path is
   * enough to escape the store, and a manifest is remote input. */
  char *save = NULL;
  for (char *tok = strtok_r(tmp, "/", &save); tok; tok = strtok_r(NULL, "/", &save))
    if (!yame_assets_safe_name(tok)) return 0;
  return 1;
}

int yame_assets_join(char *out, size_t n, const char *a, const char *b) {
  if (!a || !*a) return snprintf(out, n, "%s", b ? b : "") >= (int)n ? -1 : 0;
  if (!b || !*b) return snprintf(out, n, "%s", a) >= (int)n ? -1 : 0;

  size_t la = strlen(a);
  int wrote = (la && a[la-1] == '/') ? snprintf(out, n, "%s%s", a, b)
                                     : snprintf(out, n, "%s/%s", a, b);
  return wrote >= (int)n ? -1 : 0;
}

int yame_assets_mkdir_p(const char *path) {
  char tmp[YAME_PATH_MAX];
  if (snprintf(tmp, sizeof(tmp), "%s", path) >= (int)sizeof(tmp)) return -1;

  size_t n = strlen(tmp);
  if (n && tmp[n-1] == '/') tmp[n-1] = '\0';

  for (char *p = tmp + 1; *p; ++p) {
    if (*p != '/') continue;
    *p = '\0';
    if (mkdir(tmp, 0755) != 0 && errno != EEXIST) return -1;
    *p = '/';
  }
  if (mkdir(tmp, 0755) != 0 && errno != EEXIST) return -1;
  return 0;
}

int yame_assets_mkdir_parent(const char *path) {
  char tmp[YAME_PATH_MAX];
  if (snprintf(tmp, sizeof(tmp), "%s", path) >= (int)sizeof(tmp)) return -1;

  char *slash = strrchr(tmp, '/');
  if (!slash) return 0;               /* no parent component to create */
  *slash = '\0';
  return *tmp ? yame_assets_mkdir_p(tmp) : 0;
}

int yame_assets_is_file(const char *path) {
  struct stat st;
  return stat(path, &st) == 0 && S_ISREG(st.st_mode);
}

int yame_assets_is_dir(const char *path) {
  struct stat st;
  return stat(path, &st) == 0 && S_ISDIR(st.st_mode);
}

int yame_assets_exe_dir(char *out, size_t n) {
#if defined(__APPLE__)
  uint32_t sz = (uint32_t)n;
  if (_NSGetExecutablePath(out, &sz) != 0) return -1;
#else
  ssize_t got = readlink("/proc/self/exe", out, n - 1);
  if (got <= 0) return -1;
  out[got] = '\0';
#endif
  char *slash = strrchr(out, '/');
  if (!slash) return -1;
  *slash = '\0';
  return 0;
}

/* Legacy per-tool caches, from before the store was shared. Never read from --
 * see yame_assets_legacy_notice(). */
static void legacy_dirs(char v[][YAME_PATH_MAX], size_t cap, size_t *n) {
  *n = 0;
  const char *home = getenv("HOME");
  const char *xdg  = getenv("XDG_CACHE_HOME");

  if (xdg && *xdg && *n < cap)
    snprintf(v[(*n)++], YAME_PATH_MAX, "%s/sesame", xdg);
  if (home && *home) {
    if (*n < cap) snprintf(v[(*n)++], YAME_PATH_MAX, "%s/.cache/kycg", home);
    if (*n < cap) snprintf(v[(*n)++], YAME_PATH_MAX, "%s/.cache/sesame", home);
    /* sesame-cli uses the Apple cache dir on Darwin, so a mac user would
     * otherwise get no notice at all. Checked on every platform: it costs a
     * stat, and a store copied from a mac is a real thing. */
    if (*n < cap) snprintf(v[(*n)++], YAME_PATH_MAX, "%s/Library/Caches/sesame", home);
  }
}

void yame_assets_legacy_notice(const char *root) {
  static int said = 0;
  if (said || !root || yame_assets_is_dir(root)) return;

  char v[4][YAME_PATH_MAX];
  size_t n = 0;
  legacy_dirs(v, 4, &n);

  for (size_t i = 0; i < n; ++i) {
    if (!yame_assets_is_dir(v[i])) continue;
    fprintf(stderr,
            "[yame] note: assets now live in %s; the older per-tool cache at %s\n"
            "       is left untouched and is no longer read. Set YAME_DATA_HOME to\n"
            "       override, or delete the old cache once you are happy.\n",
            root, v[i]);
    said = 1;
    return;
  }
}

const char *yame_assets_root(const char *override, const char *tool_env,
                             char *buf, size_t n) {
  const char *env;

  if (override && *override) {
    snprintf(buf, n, "%s", override);
    return buf;
  }
  if (tool_env && *tool_env && (env = getenv(tool_env)) && *env) {
    snprintf(buf, n, "%s", env);
    return buf;
  }
  if ((env = getenv("YAME_DATA_HOME")) && *env) {
    snprintf(buf, n, "%s", env);
    yame_assets_legacy_notice(buf);
    return buf;
  }
  if ((env = getenv("XDG_DATA_HOME")) && *env) {
    snprintf(buf, n, "%s/yame", env);
    yame_assets_legacy_notice(buf);
    return buf;
  }

  const char *home = getenv("HOME");
  if (!home || !*home) home = ".";
  snprintf(buf, n, "%s/.local/share/yame", home);
  yame_assets_legacy_notice(buf);
  return buf;
}

int yame_assets_root_writable(const char *root) {
  if (!root || !*root) return 0;
  if (yame_assets_is_dir(root)) return access(root, W_OK) == 0;

  /* Not there yet: writable if some existing ancestor is, since that is where
   * mkdir -p would start. */
  char tmp[YAME_PATH_MAX];
  if (snprintf(tmp, sizeof(tmp), "%s", root) >= (int)sizeof(tmp)) return 0;

  for (;;) {
    char *slash = strrchr(tmp, '/');
    if (!slash || slash == tmp) return access("/", W_OK) == 0;
    *slash = '\0';
    if (yame_assets_is_dir(tmp)) return access(tmp, W_OK) == 0;
  }
}

/* --------------------------------------------------------------- manifests */

yame_sums_ent_t *yame_assets_parse_sums(const char *text, size_t *n) {
  size_t cap = 64, cnt = 0;
  yame_sums_ent_t *v = malloc(cap * sizeof(yame_sums_ent_t));
  if (!v) return NULL;

  const char *p = text;
  while (*p) {
    const char *eol = strchr(p, '\n');
    size_t len = eol ? (size_t)(eol - p) : strlen(p);

    /* "<64 hex>  <name>" -- shorter than that cannot carry both. */
    if (len > 66) {
      if (cnt == cap) {
        cap *= 2;
        yame_sums_ent_t *nv = realloc(v, cap * sizeof(yame_sums_ent_t));
        if (!nv) { free(v); return NULL; }
        v = nv;
      }
      memcpy(v[cnt].sha, p, 64);
      v[cnt].sha[64] = '\0';

      const char *q = p + 64;
      while ((size_t)(q - p) < len && (*q == ' ' || *q == '*' || *q == '\t')) ++q;
      size_t nlen = len - (size_t)(q - p);
      while (nlen && (q[nlen-1] == '\r' || q[nlen-1] == ' ')) --nlen;
      if (nlen >= sizeof(v[cnt].name)) nlen = sizeof(v[cnt].name) - 1;
      memcpy(v[cnt].name, q, nlen);
      v[cnt].name[nlen] = '\0';

      /* The name becomes a path under the store, so it is validated here
       * rather than at every join. sesame-cli's fetch_subtree omitted this
       * check; a "../../x" entry in a manifest wrote outside its store. */
      if (yame_assets_safe_relpath(v[cnt].name)) ++cnt;
    }

    if (!eol) break;
    p = eol + 1;
  }
  *n = cnt;
  return v;
}

yame_sums_ent_t *yame_assets_sums_load_file(const char *path, size_t *n) {
  *n = 0;
  FILE *fp = fopen(path, "rb");
  if (!fp) return NULL;

  if (fseek(fp, 0, SEEK_END) != 0) { fclose(fp); return NULL; }
  long sz = ftell(fp);
  if (sz < 0) { fclose(fp); return NULL; }
  rewind(fp);

  char *text = malloc((size_t)sz + 1);
  if (!text) { fclose(fp); return NULL; }
  size_t got = fread(text, 1, (size_t)sz, fp);
  fclose(fp);
  text[got] = '\0';

  yame_sums_ent_t *v = yame_assets_parse_sums(text, n);
  free(text);
  return v;
}

/* ----------------------------------------------------------------- the pin */

int yame_assets_pin_check(const char *dir, const char *anchor_sha) {
  char sums[YAME_PATH_MAX];
  if (yame_assets_join(sums, sizeof(sums), dir, YAME_ASSETS_SUMS_FILE) != 0)
    return YAME_PIN_ABSENT;
  if (!yame_assets_is_file(sums)) return YAME_PIN_ABSENT;
  if (!anchor_sha || !*anchor_sha) return YAME_PIN_UNKNOWN;

  char got[65];
  if (yame_assets_sha256_file(sums, got) != 0) return YAME_PIN_ABSENT;
  return yame_assets_digest_equal(got, anchor_sha) ? YAME_PIN_MATCH
                                                   : YAME_PIN_CONFLICT;
}

/* --------------------------------------------------------------- fetching */

#ifdef YAME_HAVE_CURL

typedef struct { char *s; size_t n, m; } membuf_t;

static size_t mem_write(void *data, size_t sz, size_t nm, void *ud) {
  membuf_t *b = ud;
  size_t add = sz * nm;
  if (b->n + add + 1 > b->m) {
    size_t want = (b->n + add + 1) * 2;
    char *p = realloc(b->s, want);
    if (!p) return 0;
    b->s = p; b->m = want;
  }
  memcpy(b->s + b->n, data, add);
  b->n += add;
  b->s[b->n] = '\0';
  return add;
}

static CURL *new_handle(const char *url) {
  CURL *h = curl_easy_init();
  if (!h) return NULL;
  curl_easy_setopt(h, CURLOPT_URL, url);
  curl_easy_setopt(h, CURLOPT_FAILONERROR, 1L);
  curl_easy_setopt(h, CURLOPT_NOSIGNAL, 1L);
  curl_easy_setopt(h, CURLOPT_CONNECTTIMEOUT, 30L);
  curl_easy_setopt(h, CURLOPT_USERAGENT, "yame/" YAME_VERSION);

  /* Redirects must be followed -- github.com/<repo>/raw/... legitimately lands
   * on raw.githubusercontent.com -- but are confined to https and bounded.
   * Every byte is digest-checked against a compiled-in anchor, so a downgrade
   * could not substitute content; what it could do is move the transfer to
   * cleartext, exposing which references are being fetched. */
  curl_easy_setopt(h, CURLOPT_FOLLOWLOCATION, 1L);
  curl_easy_setopt(h, CURLOPT_MAXREDIRS, 10L);
#if defined(CURLOPT_REDIR_PROTOCOLS_STR) && LIBCURL_VERSION_NUM >= 0x075500
  curl_easy_setopt(h, CURLOPT_REDIR_PROTOCOLS_STR, "https");
#elif defined(CURLOPT_REDIR_PROTOCOLS)
  /* Pre-7.85 spelling; deprecated, but it is what older libcurl understands
   * and bioconda still builds against those. */
  curl_easy_setopt(h, CURLOPT_REDIR_PROTOCOLS, (long)CURLPROTO_HTTPS);
#endif
  return h;
}

int yame_assets_have_curl(void) { return 1; }

char *yame_assets_http_get_mem(const char *url, size_t *len) {
  CURL *h = new_handle(url);
  if (!h) return NULL;

  membuf_t b = {0};
  curl_easy_setopt(h, CURLOPT_WRITEFUNCTION, mem_write);
  curl_easy_setopt(h, CURLOPT_WRITEDATA, &b);

  CURLcode rc = curl_easy_perform(h);
  curl_easy_cleanup(h);

  if (rc != CURLE_OK) { free(b.s); return NULL; }
  if (len) *len = b.n;
  return b.s;
}

typedef struct { const yame_fetch_opt_t *opt; } xfer_ctx_t;

static int on_xfer(void *ud, curl_off_t dltotal, curl_off_t dlnow,
                   curl_off_t ultotal, curl_off_t ulnow) {
  (void)ultotal; (void)ulnow;
  xfer_ctx_t *x = ud;
  if (x && x->opt && x->opt->on_progress)
    x->opt->on_progress(x->opt->ud, (uint64_t)dlnow, (uint64_t)dltotal);
  return 0;
}

/*
 * Download to `path`, a temp name the caller renames on success.
 *
 * O_EXCL, with the pid in the name: the temp used to be a fixed "<file>.part"
 * in both kycg and sesame, so two concurrent fetches into one store wrote the
 * same file. kycg hit this and fixed it; sesame still has it. Consolidating the
 * store makes it worse, because the two racers can now be different tools, so
 * the fixed version is the one that moved here.
 */
static int http_get_file(const char *url, const char *path,
                         const yame_fetch_opt_t *opt) {
  int fd = open(path, O_WRONLY | O_CREAT | O_EXCL, 0644);
  if (fd < 0) return -1;
  FILE *fp = fdopen(fd, "wb");
  if (!fp) { close(fd); unlink(path); return -1; }

  CURL *h = new_handle(url);
  if (!h) { fclose(fp); unlink(path); return -1; }

  xfer_ctx_t x = { opt };
  curl_easy_setopt(h, CURLOPT_WRITEDATA, fp);
  if (opt && opt->on_progress) {
    curl_easy_setopt(h, CURLOPT_XFERINFOFUNCTION, on_xfer);
    curl_easy_setopt(h, CURLOPT_XFERINFODATA, &x);
    curl_easy_setopt(h, CURLOPT_NOPROGRESS, 0L);
  }

  CURLcode rc = curl_easy_perform(h);
  curl_easy_cleanup(h);
  fclose(fp);

  if (rc != CURLE_OK) { unlink(path); return -1; }
  return 0;
}

/*
 * Temp files from a killed process, which unlink their own on any ordinary
 * failure. Only those older than a day go, so a fetch running right now in
 * another process is never touched.
 */
static void sweep_stale_parts(const char *dir) {
  DIR *d = opendir(dir);
  if (!d) return;

  time_t now = time(NULL);
  struct dirent *e;
  while ((e = readdir(d))) {
    size_t l = strlen(e->d_name);
    if (l < 6 || strcmp(e->d_name + l - 5, ".part") != 0) continue;

    char p[YAME_PATH_MAX];
    if (yame_assets_join(p, sizeof(p), dir, e->d_name) != 0) continue;
    struct stat st;
    if (stat(p, &st) == 0 && S_ISREG(st.st_mode) &&
        now - st.st_mtime > 24 * 60 * 60)
      unlink(p);
  }
  closedir(d);
}

int yame_assets_download_verify(const char *url, const char *want_sha,
                                const char *dest, const yame_fetch_opt_t *opt,
                                int *downloaded, char **err) {
  if (downloaded) *downloaded = 0;

  const char *base = strrchr(dest, '/');
  base = base ? base + 1 : dest;

  /* Already here and correct? Then this is a no-op, which is what makes one
   * tool's fetch satisfy the next one's. */
  if (!(opt && opt->force) && want_sha && *want_sha &&
      yame_assets_is_file(dest)) {
    char got[65];
    if (yame_assets_sha256_file(dest, got) == 0 &&
        yame_assets_digest_equal(got, want_sha))
      return 0;
  }

  if (yame_assets_mkdir_parent(dest) != 0) {
    set_err(err, "cannot create the directory for %s", dest);
    return -1;
  }

  char part[YAME_PATH_MAX];
  if (snprintf(part, sizeof(part), "%s.%ld.part", dest, (long)getpid())
      >= (int)sizeof(part)) {
    set_err(err, "path too long: %s", dest);
    return -1;
  }

  if (opt && opt->on_begin) opt->on_begin(opt->ud, base, 0);

  if (http_get_file(url, part, opt) != 0) {
    set_err(err, "download failed: %s", url);
    if (opt && opt->on_done) opt->on_done(opt->ud, base, 0, 0);
    return -1;
  }

  /* No digest, no rename. Every channel publishes one, so a missing digest is
   * a bug in the caller rather than a permissive case. */
  int ok = 0;
  char got[65];
  if (want_sha && *want_sha) {
    if (yame_assets_sha256_file(part, got) == 0)
      ok = yame_assets_digest_equal(got, want_sha);
  }
  if (!ok) {
    unlink(part);
    set_err(err, "digest mismatch, discarded: %s", base);
    if (opt && opt->on_done) opt->on_done(opt->ud, base, 0, 0);
    return -1;
  }

  struct stat st;
  uint64_t sz = (stat(part, &st) == 0) ? (uint64_t)st.st_size : 0;

  if (rename(part, dest) != 0) {
    unlink(part);
    set_err(err, "cannot move into the store: %s", dest);
    if (opt && opt->on_done) opt->on_done(opt->ud, base, 0, 0);
    return -1;
  }

  if (downloaded) *downloaded = 1;
  if (opt && opt->on_done) opt->on_done(opt->ud, base, sz, 1);
  return 0;
}

int yame_assets_fetch_subtree(const char *base, const char *tag,
                              const char *remote_sub, const char *store_sub,
                              const char *anchor_sha,
                              const yame_fetch_opt_t *opt, char **err) {
  /* Whose tag filled this directory? A conflict means another tool, or an
   * older build of this one, populated it from a tag this build does not pin.
   * Overwriting would start a re-download war between the two binaries, so
   * stop and let the human decide. */
  int pin = yame_assets_pin_check(store_sub, anchor_sha);
  if (pin == YAME_PIN_CONFLICT && !(opt && opt->force)) {
    set_err(err,
            "%s holds a %s this build does not pin: it was populated from a "
            "different upstream tag (this build pins %s). Rebuild against the "
            "same tag, point the tool's own data-dir variable at a private "
            "root, or re-fetch with -f to overwrite.",
            store_sub, YAME_ASSETS_SUMS_FILE, tag ? tag : "(none)");
    return -1;
  }

  if (yame_assets_mkdir_p(store_sub) != 0) {
    set_err(err, "cannot create %s", store_sub);
    return -1;
  }
  sweep_stale_parts(store_sub);

  /* The manifest first: everything else is verified against it, and it is
   * verified against the anchor compiled into the caller. */
  char url[YAME_PATH_MAX];
  int wrote = remote_sub && *remote_sub
    ? snprintf(url, sizeof(url), "%s/%s/%s/%s", base, tag, remote_sub,
               YAME_ASSETS_SUMS_FILE)
    : snprintf(url, sizeof(url), "%s/%s/%s", base, tag, YAME_ASSETS_SUMS_FILE);
  if (wrote >= (int)sizeof(url)) { set_err(err, "URL too long"); return -1; }

  size_t sums_len = 0;
  char *sums_text = yame_assets_http_get_mem(url, &sums_len);
  if (!sums_text) {
    set_err(err, "cannot fetch %s", url);
    return -1;
  }

  if (anchor_sha && *anchor_sha) {
    char got[65];
    yame_assets_sha256_buf(sums_text, sums_len, got);
    if (!yame_assets_digest_equal(got, anchor_sha)) {
      free(sums_text);
      set_err(err,
              "%s at %s does not match the digest this build pins. Either the "
              "upstream tag moved, or this is not the content this build "
              "expects; nothing was written.",
              YAME_ASSETS_SUMS_FILE, url);
      return -1;
    }
  }

  size_t n = 0;
  yame_sums_ent_t *ents = yame_assets_parse_sums(sums_text, &n);
  if (!ents || n == 0) {
    free(sums_text); free(ents);
    set_err(err, "empty or unreadable %s at %s", YAME_ASSETS_SUMS_FILE, url);
    return -1;
  }

  int failed = 0;
  for (size_t i = 0; i < n; ++i) {
    char dest[YAME_PATH_MAX], furl[YAME_PATH_MAX];
    if (yame_assets_join(dest, sizeof(dest), store_sub, ents[i].name) != 0)
      { ++failed; continue; }
    wrote = remote_sub && *remote_sub
      ? snprintf(furl, sizeof(furl), "%s/%s/%s/%s", base, tag, remote_sub,
                 ents[i].name)
      : snprintf(furl, sizeof(furl), "%s/%s/%s", base, tag, ents[i].name);
    if (wrote >= (int)sizeof(furl)) { ++failed; continue; }

    char *ferr = NULL;
    if (yame_assets_download_verify(furl, ents[i].sha, dest, opt, NULL, &ferr)
        != 0) {
      if (!(opt && opt->quiet) && ferr) fprintf(stderr, "[yame] %s\n", ferr);
      ++failed;
    }
    free(ferr);
  }

  /* Keep the manifest verbatim: it is what lets `shasum -a 256 -c SHA256SUMS`
   * re-verify the store with none of this code, and it is what the next tool
   * reads to learn which tag filled this directory. */
  if (!failed) {
    char sp[YAME_PATH_MAX];
    if (yame_assets_join(sp, sizeof(sp), store_sub, YAME_ASSETS_SUMS_FILE) == 0) {
      FILE *fp = fopen(sp, "wb");
      if (fp) { fwrite(sums_text, 1, sums_len, fp); fclose(fp); }
    }
  }

  free(sums_text);
  free(ents);

  if (failed) {
    set_err(err, "%d of %zu files failed; the manifest was not written",
            failed, n);
    return -1;
  }
  return 0;
}

#else  /* !YAME_HAVE_CURL */

int yame_assets_have_curl(void) { return 0; }

char *yame_assets_http_get_mem(const char *url, size_t *len) {
  (void)url; (void)len;
  return NULL;
}

static const char NO_CURL[] =
  "this build has no libcurl, so it cannot download. Rebuild with libcurl "
  "available (see the Makefile's CURL block), or install the files by hand";

int yame_assets_download_verify(const char *url, const char *want_sha,
                                const char *dest, const yame_fetch_opt_t *opt,
                                int *downloaded, char **err) {
  (void)want_sha; (void)opt;
  if (downloaded) *downloaded = 0;
  set_err(err, "%s: %s -> %s", NO_CURL, url, dest);
  return -1;
}

int yame_assets_fetch_subtree(const char *base, const char *tag,
                              const char *remote_sub, const char *store_sub,
                              const char *anchor_sha,
                              const yame_fetch_opt_t *opt, char **err) {
  (void)base; (void)tag; (void)remote_sub; (void)anchor_sha; (void)opt;
  set_err(err, "%s: %s", NO_CURL, store_sub);
  return -1;
}

#endif /* YAME_HAVE_CURL */
