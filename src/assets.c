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

const char *yame_assets_index_suffix(const char *name) {
  static const char *const sfx[] = { ".idx", ".tbi", NULL };
  if (!name) return NULL;
  size_t l = strlen(name);
  for (size_t k = 0; sfx[k]; ++k) {
    size_t sl = strlen(sfx[k]);
    if (l > sl && strcmp(name + l - sl, sfx[k]) == 0) return sfx[k];
  }
  return NULL;
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

/* The store variables the suite's other tools export. The store is one shared
 * tree, so the variable that moves it for one tool moves it for all of them:
 * a reader who exported METHSCOPE_DATA_HOME and then runs `yame hprint` means
 * that store, and finding nothing there is a bug, not a policy. Read after
 * $YAME_DATA_HOME, so yame's own variable still wins for yame. */
const char *const yame_assets_suite_env[] = {
  "METHSCOPE_DATA_HOME",
  NULL
};

const char *yame_assets_root_env(const char *tool_env) {
  const char *env;
  size_t i;

  if (tool_env && *tool_env && (env = getenv(tool_env)) && *env) return tool_env;
  if ((env = getenv("YAME_DATA_HOME")) && *env) return "YAME_DATA_HOME";
  for (i = 0; yame_assets_suite_env[i]; ++i) {
    if (tool_env && strcmp(tool_env, yame_assets_suite_env[i]) == 0) continue;
    if ((env = getenv(yame_assets_suite_env[i])) && *env)
      return yame_assets_suite_env[i];
  }
  return NULL;
}

const char *yame_assets_root(const char *override, const char *tool_env,
                             char *buf, size_t n) {
  const char *env, *var;

  if (override && *override) {
    snprintf(buf, n, "%s", override);
    return buf;
  }
  if ((var = yame_assets_root_env(tool_env)) != NULL) {
    snprintf(buf, n, "%s", getenv(var));
    return buf;
  }
  if ((env = getenv("XDG_DATA_HOME")) && *env) {
    snprintf(buf, n, "%s/yame", env);
    return buf;
  }

  const char *home = getenv("HOME");
  if (!home || !*home) home = ".";
  snprintf(buf, n, "%s/.local/share/yame", home);
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

/* ------------------------------------------------------------------ keys */

/* source@tag:remote_path, split on the LAST @ and the LAST colon. */
static size_t put(char *out, size_t n, const char *p, size_t l) {
  if (n) { size_t c = l < n - 1 ? l : n - 1; memcpy(out, p, c); out[c] = '\0'; }
  return l;
}
size_t yame_key_path(const char *key, char *out, size_t n) {
  const char *c = strrchr(key, ':');
  const char *p = c ? c + 1 : key;
  return put(out, n, p, strlen(p));
}
size_t yame_key_tag(const char *key, char *out, size_t n) {
  const char *c = strrchr(key, ':');
  size_t end = c ? (size_t)(c - key) : strlen(key);
  size_t at = end;
  while (at > 0 && key[at - 1] != '@') --at;
  if (at == 0) return put(out, n, "", 0);
  return put(out, n, key + at, end - at);
}
size_t yame_key_source(const char *key, char *out, size_t n) {
  const char *c = strrchr(key, ':');
  size_t end = c ? (size_t)(c - key) : strlen(key);
  size_t at = end;
  while (at > 0 && key[at - 1] != '@') --at;
  if (at == 0) return put(out, n, key, end);
  return put(out, n, key, at - 1);
}

const char *yame_file_name(const yame_asset_file_t *f) {
  const char *s = strrchr(f->store_path, '/');
  return s ? s + 1 : f->store_path;
}
size_t yame_file_dirlen(const yame_asset_file_t *f) {
  const char *s = strrchr(f->store_path, '/');
  return s ? (size_t)(s - f->store_path) : 0;
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

/* YAME_ASSETS_MIRROR=<scheme://host[:port]> replaces the scheme and host of
 * every download URL and keeps the path, so
 *
 *   https://raw.githubusercontent.com/zhou-lab/InfiniumAnnotation/v8.1/EPIC/SHA256SUMS
 *
 * becomes <mirror>/zhou-lab/InfiniumAnnotation/v8.1/EPIC/SHA256SUMS. For a
 * site that mirrors the public repositories -- compute nodes that cannot
 * reach GitHub are the usual case -- and for the test suite, which serves a
 * loopback mirror. Verification is untouched: every byte is still checked
 * against the compiled-in digests, so a mirror can fail to serve the right
 * bytes but cannot pass wrong ones. */
static const char *mirror_url(const char *url, char *buf, size_t n) {
  const char *m = getenv("YAME_ASSETS_MIRROR");
  if (!m || !*m) return url;
  const char *p = strstr(url, "://");
  if (!p) return url;
  const char *path = strchr(p + 3, '/');
  if (!path) return url;
  size_t ml = strlen(m);
  while (ml && m[ml - 1] == '/') ml--;          /* a trailing slash is fine */
  int k = snprintf(buf, n, "%.*s%s", (int) ml, m, path);
  return (k < 0 || (size_t) k >= n) ? url : buf;
}

static CURL *new_handle(const char *url) {
  CURL *h = curl_easy_init();
  if (!h) return NULL;
  char mirrored[4096];                          /* libcurl copies the string */
  url = mirror_url(url, mirrored, sizeof mirrored);
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

/*
 * What a transfer came back with, kept so a failure can say WHY. A bare
 * "download failed: <url>" once cost a user a `curl -I` by hand to learn that
 * HuggingFace was rate-limiting them (429) -- "wait and retry" -- and not that
 * the registry was broken; nothing in the message told a 429 from a 404, a
 * moved tag, a DNS miss or a TLS problem.
 */
typedef struct {
  CURLcode rc;               /* CURLE_HTTP_RETURNED_ERROR when `http` speaks */
  long     http;             /* response code, 0 when the request never got one */
  long     retry_after;      /* seconds the server asked for, 0 when it did not */
  int      tries;            /* performs made, 1 when nothing was retried */
} http_res_t;

static const char *http_phrase(long code) {
  switch (code) {
  case 400: return "Bad Request";       case 401: return "Unauthorized";
  case 403: return "Forbidden";         case 404: return "Not Found";
  case 408: return "Request Timeout";   case 410: return "Gone";
  case 429: return "Too Many Requests"; case 500: return "Internal Server Error";
  case 502: return "Bad Gateway";       case 503: return "Service Unavailable";
  case 504: return "Gateway Timeout";   default:  return NULL;
  }
}

/* "HTTP 429 Too Many Requests (Retry-After: 30 s)", or the libcurl text for
 * a failure below HTTP: "curl: Couldn't resolve host name". */
static void http_why(const http_res_t *r, char *buf, size_t n) {
  if (r->rc == CURLE_HTTP_RETURNED_ERROR || (r->rc == CURLE_OK && r->http >= 400)) {
    const char *ph = http_phrase(r->http);
    int k = snprintf(buf, n, "HTTP %ld%s%s", r->http, ph ? " " : "", ph ? ph : "");
    if (r->retry_after > 0 && k > 0 && (size_t) k < n)
      snprintf(buf + k, n - (size_t) k, " (Retry-After: %ld s)", r->retry_after);
  } else {
    snprintf(buf, n, "curl: %s", curl_easy_strerror(r->rc));
  }
}

/* A transient status, one the same request will likely clear on its own. */
static int http_transient(long code) { return code == 429 || code == 503; }

/*
 * Bounded automatic retry, for 429 and 503 only. Three tries, waiting what
 * the server asked for in Retry-After when it said, else 30 s then 60 s. A
 * Retry-After beyond two minutes is reported, not slept: a user can decide
 * to wait ten minutes, a fetch should not decide it for them. Every other
 * status fails at once -- a 404 does not get better by asking again.
 *
 * `reset` puts the sink back to empty before the next try. libcurl under
 * FAILONERROR stops at the status line and delivers no error body, so
 * nothing should have landed, but the sink is truncated regardless.
 */
static CURLcode http_perform(CURL *h, http_res_t *r,
                             void (*reset)(void *), void *sink) {
  const int TRIES = 3;
  const long backoff[] = { 30, 60 };
  memset(r, 0, sizeof *r);
  for (;;) {
    r->tries++;
    r->rc = curl_easy_perform(h);
    r->http = 0; r->retry_after = 0;
    curl_easy_getinfo(h, CURLINFO_RESPONSE_CODE, &r->http);
#if LIBCURL_VERSION_NUM >= 0x074200                    /* 7.66: CURLINFO_RETRY_AFTER */
    { curl_off_t ra = 0;
      if (curl_easy_getinfo(h, CURLINFO_RETRY_AFTER, &ra) == CURLE_OK && ra > 0)
        r->retry_after = (long) (ra > 86400 ? 86400 : ra); }
#endif
    if (r->rc == CURLE_OK) return r->rc;
    if (r->rc != CURLE_HTTP_RETURNED_ERROR || !http_transient(r->http)) return r->rc;
    if (r->tries >= TRIES) return r->rc;

    long wait = r->retry_after > 0 ? r->retry_after : backoff[r->tries - 1];
    if (wait > 120) return r->rc;                      /* said, not slept */

    char why[128];
    http_why(r, why, sizeof why);
    /* On a terminal the progress bar owns the line; take it over cleanly. */
    fprintf(stderr, "%s[yame] %s; retrying in %ld s (try %d of %d)\n",
            isatty(STDERR_FILENO) ? "\r\033[K" : "", why, wait, r->tries + 1, TRIES);
    fflush(stderr);
    sleep((unsigned) wait);
    if (reset) reset(sink);
  }
}

static void mem_reset(void *ud) { ((membuf_t *) ud)->n = 0; }

static char *http_get_mem_why(const char *url, size_t *len, http_res_t *res) {
  CURL *h = new_handle(url);
  if (!h) return NULL;

  membuf_t b = {0};
  curl_easy_setopt(h, CURLOPT_WRITEFUNCTION, mem_write);
  curl_easy_setopt(h, CURLOPT_WRITEDATA, &b);

  http_res_t r;
  CURLcode rc = http_perform(h, &r, mem_reset, &b);
  curl_easy_cleanup(h);
  if (res) *res = r;

  if (rc != CURLE_OK) { free(b.s); return NULL; }
  if (len) *len = b.n;
  return b.s;
}

char *yame_assets_http_get_mem(const char *url, size_t *len) {
  return http_get_mem_why(url, len, NULL);
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
static void file_reset(void *ud) {
  FILE *fp = ud;
  fflush(fp);
  if (ftruncate(fileno(fp), 0) == 0) rewind(fp);
}

/* Returns 0, or -1 with `res` saying why when the transfer itself failed
 * (res->tries == 0 means it never started: a local file problem). */
static int http_get_file(const char *url, const char *path,
                         const yame_fetch_opt_t *opt, http_res_t *res) {
  memset(res, 0, sizeof *res);
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

  CURLcode rc = http_perform(h, res, file_reset, fp);
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

  http_res_t res;
  if (http_get_file(url, part, opt, &res) != 0) {
    /* Name the URL actually asked, so a mirror that is down is blamed and
     * not the upstream it stands in for. */
    char mbuf[4096];
    url = mirror_url(url, mbuf, sizeof mbuf);
    if (res.tries == 0) {
      set_err(err, "cannot create %s: %s", part, strerror(errno));
    } else {
      char why[128];
      http_why(&res, why, sizeof why);
      if (res.rc == CURLE_HTTP_RETURNED_ERROR && http_transient(res.http)) {
        /* The one case where "try again later" is the whole answer. */
        const char *ph = http_phrase(res.http);
        if (res.tries > 1)
          set_err(err, "download failed after %d tries: HTTP %ld %s, the "
                  "server is rate-limiting this client; wait a few minutes "
                  "and re-run: %s", res.tries, res.http, ph, url);
        else
          set_err(err, "download failed: HTTP %ld %s, the server asks for a "
                  "%ld s wait, longer than a fetch waits on its own (120 s); "
                  "re-run after that: %s", res.http, ph, res.retry_after, url);
      } else {
        set_err(err, "download failed: %s: %s", why, url);
      }
    }
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

/* Write a manifest without ever exposing a truncated one. The manifest is the
 * directory's tag identity, so preserve the old bytes until its replacement
 * is complete. */
static int write_manifest(const char *path, const char *text, size_t len) {
  char part[YAME_PATH_MAX];
  if (snprintf(part, sizeof(part), "%s.%ld.part", path, (long)getpid())
      >= (int)sizeof(part)) return -1;

  int fd = open(part, O_WRONLY | O_CREAT | O_EXCL, 0644);
  if (fd < 0) return -1;
  FILE *fp = fdopen(fd, "wb");
  if (!fp) { close(fd); unlink(part); return -1; }

  int ok = fwrite(text, 1, len, fp) == len;
  if (fclose(fp) != 0) ok = 0;
  if (!ok || rename(part, path) != 0) {
    unlink(part);
    return -1;
  }
  return 0;
}

/* One line of a directory's manifest, by file name: the digest recorded for
 * it, or NULL when the manifest has no line for it (or does not exist). */
static const char *manifest_digest(const yame_sums_ent_t *ents, size_t n,
                                   const char *name) {
  for (size_t i = 0; i < n; ++i)
    if (strcmp(ents[i].name, name) == 0) return ents[i].sha;
  return NULL;
}

int yame_assets_fetch_files(const yame_fetch_cfg_t *cfg, const char *root,
                            const yame_asset_file_t *const *want, size_t n_want,
                            const yame_fetch_opt_t *opt, char **err) {
  int failed = 0;
  /* verified[i] for want[i]: 1 once its bytes are known to carry the pinned
   * digest, whether they just arrived or were already there. */
  unsigned char *verified = calloc(n_want ? n_want : 1, 1);
  if (!verified) { set_err(err, "out of memory"); return -1; }

  /* The manifest of the directory last looked at, so a directory's files are
   * classified from one parse. */
  char cur_dir[YAME_PATH_MAX] = "";
  yame_sums_ent_t *ents = NULL; size_t n_ents = 0;

  for (size_t i = 0; i < n_want; ++i) {
    const yame_asset_file_t *f = want[i];
    char dest[YAME_PATH_MAX], dir[YAME_PATH_MAX];
    if (yame_assets_join(dest, sizeof dest, root, f->store_path) != 0) { ++failed; continue; }
    size_t dl = yame_file_dirlen(f);
    if (dl >= sizeof dir) { ++failed; continue; }
    if (yame_assets_join(dir, sizeof dir, root, "") != 0) { ++failed; continue; }
    { char sub[YAME_PATH_MAX]; memcpy(sub, f->store_path, dl); sub[dl] = '\0';
      if (yame_assets_join(dir, sizeof dir, root, sub) != 0) { ++failed; continue; } }

    if (strcmp(dir, cur_dir) != 0) {
      free(ents); ents = NULL; n_ents = 0;
      char mp[YAME_PATH_MAX];
      if (yame_assets_join(mp, sizeof mp, dir, YAME_ASSETS_SUMS_FILE) == 0)
        ents = yame_assets_sums_load_file(mp, &n_ents);
      snprintf(cur_dir, sizeof cur_dir, "%s", dir);
      if (yame_assets_mkdir_p(dir) != 0) { set_err(err, "cannot create %s", dir); ++failed; continue; }
      sweep_stale_parts(dir);
    }

    /* Recorded as current, and on disk: nothing to move, nothing to hash. */
    const char *rec = manifest_digest(ents, n_ents, yame_file_name(f));
    if (!(opt && opt->force) && rec && yame_assets_digest_equal(rec, f->sha256)
        && yame_assets_is_file(dest)) {
      verified[i] = 1;
      continue;
    }
    /* Recorded at ANOTHER digest, and on disk: stale -- fetched by a build
     * that pinned something else, or by another tool. Replacing it is what
     * -f is for; without it, say so and leave the file alone. */
    if (!(opt && opt->force) && rec && !yame_assets_digest_equal(rec, f->sha256)
        && yame_assets_is_file(dest)) {
      set_err(err, "stale: %s is recorded at a different digest than this build "
              "pins; re-run with -f to replace it", f->store_path);
      ++failed;
      continue;
    }
    /* Otherwise download_verify decides: it hashes a present file and skips
     * a match, and replaces anything else only once the new bytes verify. */
    char *ferr = NULL;
    if (yame_assets_download_verify(f->url, f->sha256, dest, opt, NULL, &ferr) != 0) {
      /* One file asked for, one reason to give: it is the whole answer, so
       * it goes back to the caller rather than to stderr with a count after
       * it. Several files each report on their own line and the count sums
       * them up. */
      if (n_want == 1 && ferr && err && !*err) { *err = ferr; ferr = NULL; }
      else if (!(opt && opt->quiet) && ferr) fprintf(stderr, "[yame] %s\n", ferr);
      ++failed;
    } else {
      verified[i] = 1;
    }
    free(ferr);
  }
  free(ents); ents = NULL; n_ents = 0;

  /* The manifests, one per directory touched: every registry file in the
   * directory, in registry order; the pinned digest where this call
   * verified the file, the old line where it did not and one existed, the
   * pinned digest where the file is not on disk at all. */
  for (size_t i = 0; i < n_want; ++i) {
    size_t dl = yame_file_dirlen(want[i]);
    int done = 0;
    for (size_t k = 0; k < i && !done; ++k)
      if (yame_file_dirlen(want[k]) == dl &&
          strncmp(want[k]->store_path, want[i]->store_path, dl) == 0) done = 1;
    if (done) continue;

    char sub[YAME_PATH_MAX], dir[YAME_PATH_MAX], mp[YAME_PATH_MAX];
    memcpy(sub, want[i]->store_path, dl); sub[dl] = '\0';
    if (yame_assets_join(dir, sizeof dir, root, sub) != 0 ||
        yame_assets_join(mp, sizeof mp, dir, YAME_ASSETS_SUMS_FILE) != 0) { ++failed; continue; }
    size_t n_old = 0;
    yame_sums_ent_t *old = yame_assets_sums_load_file(mp, &n_old);

    size_t cap = 4096, len = 0;
    char *text = malloc(cap);
    if (!text) { free(old); ++failed; continue; }
    for (size_t j = 0; j < cfg->n_files; ++j) {
      const yame_asset_file_t *g = &cfg->files[j];
      if (yame_file_dirlen(g) != dl || strncmp(g->store_path, sub, dl) != 0) continue;
      const char *name = yame_file_name(g);
      const char *sha = g->sha256;
      int v = 0;
      for (size_t k = 0; k < n_want; ++k) if (want[k] == g) { v = verified[k]; break; }
      if (!v) {
        const char *rec = manifest_digest(old, n_old, name);
        char full[YAME_PATH_MAX];
        if (rec && yame_assets_join(full, sizeof full, dir, name) == 0 &&
            yame_assets_is_file(full)) sha = rec;
      }
      size_t need = len + 64 + 2 + strlen(name) + 2;
      if (need > cap) { cap = need * 2; char *t = realloc(text, cap); if (!t) { free(text); text = NULL; break; } text = t; }
      len += (size_t)snprintf(text + len, cap - len, "%s  %s\n", sha, name);
    }
    free(old);
    if (!text) { ++failed; continue; }
    if (write_manifest(mp, text, len) != 0) {
      set_err(err, "files arrived, but cannot write %s", mp);
      ++failed;
    }
    free(text);
  }

  free(verified);
  if (failed) {
    if (!err || !*err) set_err(err, "%d of %zu files failed", failed, n_want);
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

int yame_assets_fetch_files(const yame_fetch_cfg_t *cfg, const char *root,
                            const yame_asset_file_t *const *want, size_t n_want,
                            const yame_fetch_opt_t *opt, char **err) {
  (void)cfg; (void)want; (void)n_want; (void)opt;
  set_err(err, "%s: %s", NO_CURL, root);
  return -1;
}

#endif /* YAME_HAVE_CURL */

/* ------------------------------------------ is the store what this build pins? */

yame_store_state_t yame_file_state(const char *root, const yame_asset_file_t *f) {
  char full[YAME_PATH_MAX], dir[YAME_PATH_MAX], mp[YAME_PATH_MAX], sub[YAME_PATH_MAX];
  if (yame_assets_join(full, sizeof full, root, f->store_path) != 0) return YAME_STORE_ABSENT;
  if (!yame_assets_is_file(full)) return YAME_STORE_ABSENT;
  size_t dl = yame_file_dirlen(f);
  memcpy(sub, f->store_path, dl); sub[dl] = '\0';
  if (yame_assets_join(dir, sizeof dir, root, sub) != 0 ||
      yame_assets_join(mp, sizeof mp, dir, YAME_ASSETS_SUMS_FILE) != 0) return YAME_STORE_CURRENT;
  size_t n = 0;
  yame_sums_ent_t *ents = yame_assets_sums_load_file(mp, &n);
  yame_store_state_t st = YAME_STORE_CURRENT;
  const char *name = yame_file_name(f);
  for (size_t i = 0; i < n; ++i)
    if (strcmp(ents[i].name, name) == 0) {
      if (!yame_assets_digest_equal(ents[i].sha, f->sha256)) st = YAME_STORE_STALE;
      break;
    }
  free(ents);
  return st;
}

/* The state of one store directory: how many of its registry files are on
 * disk, and how many of those are stale. */
static void dir_states(const yame_fetch_cfg_t *cfg, const char *root,
                       const char *sub, size_t *n, size_t *present, size_t *stale) {
  size_t dl = strlen(sub);
  *n = *present = *stale = 0;
  for (size_t j = 0; j < cfg->n_files; ++j) {
    const yame_asset_file_t *g = &cfg->files[j];
    if (yame_file_dirlen(g) != dl || strncmp(g->store_path, sub, dl) != 0) continue;
    ++*n;
    switch (yame_file_state(root, g)) {
    case YAME_STORE_CURRENT: ++*present; break;
    case YAME_STORE_STALE:   ++*present; ++*stale; break;
    default: break;
    }
  }
}

/* The upstream a directory's files come from, for prose: the source out of
 * the key, without the hf: host prefix. */
static void dir_upstream(const yame_fetch_cfg_t *cfg, const char *sub,
                         char *out, size_t n) {
  size_t dl = strlen(sub);
  out[0] = '\0';
  for (size_t j = 0; j < cfg->n_files; ++j) {
    const yame_asset_file_t *g = &cfg->files[j];
    if (yame_file_dirlen(g) != dl || strncmp(g->store_path, sub, dl) != 0) continue;
    yame_key_source(g->key, out, n);
    if (strncmp(out, "hf:", 3) == 0) memmove(out, out + 3, strlen(out + 3) + 1);
    return;
  }
}

/* What to say about a directory holding stale files: which files, from
 * what, and the one command that replaces them. "Differ from this build"
 * was the old text, and nobody could tell from it what had happened. */
static void stale_dir_advice(const yame_fetch_cfg_t *cfg, const char *root,
                             const char *sub, size_t stale, size_t cnt,
                             char *out, size_t n) {
  const char *tool = cfg->tool ? cfg->tool : "yame";
  char up[160], names[512] = "";
  dir_upstream(cfg, sub, up, sizeof up);
  size_t dl = strlen(sub), listed = 0;
  for (size_t j = 0; j < cfg->n_files; ++j) {
    const yame_asset_file_t *g = &cfg->files[j];
    if (yame_file_dirlen(g) != dl || strncmp(g->store_path, sub, dl) != 0) continue;
    if (yame_file_state(root, g) != YAME_STORE_STALE) continue;
    if (listed == 3) { strncat(names, ", ...", sizeof names - strlen(names) - 1); break; }
    if (listed) strncat(names, ", ", sizeof names - strlen(names) - 1);
    strncat(names, yame_file_name(g), sizeof names - strlen(names) - 1);
    ++listed;
  }
  /* "this build", not "this <tool>": a downstream tool pins the digests
   * but yame fetches, so naming one program in both halves was wrong for
   * every tool but yame itself. `tool` names the command only. */
  snprintf(out, n, "%s: %zu of %zu files come from an earlier release of %s "
           "than this build pins (%s); replace them with: %s fetch -y -f %s",
           sub, stale, cnt, up, names, tool, sub);
}

yame_store_state_t yame_store_state(const yame_fetch_cfg_t *cfg, const char *path,
                                    char *advice, size_t n) {
  if (advice && n) advice[0] = '\0';
  const char *tool = cfg->tool ? cfg->tool : "yame";
  char root[YAME_PATH_MAX];
  yame_assets_root(NULL, cfg->tool_env, root, sizeof root);

  const char *rel = path;
  size_t rl = strlen(root);
  if (strncmp(path, root, rl) == 0 && (path[rl] == '/' || path[rl] == '\0'))
    rel = path + rl + (path[rl] == '/');

  /* A file? */
  for (size_t j = 0; j < cfg->n_files; ++j) {
    const yame_asset_file_t *f = &cfg->files[j];
    if (strcmp(f->store_path, rel) != 0) continue;
    yame_store_state_t st = yame_file_state(root, f);
    if (advice) {
      if (st == YAME_STORE_ABSENT)
        snprintf(advice, n, "%s is not in the store; run: %s fetch -y %s", rel, tool, rel);
      else if (st == YAME_STORE_STALE) {
        char sub[YAME_PATH_MAX], up[160];
        size_t dl = yame_file_dirlen(f);
        memcpy(sub, f->store_path, dl); sub[dl] = '\0';
        dir_upstream(cfg, sub, up, sizeof up);
        snprintf(advice, n, "%s comes from an earlier release of %s than this "
                 "build pins; replace it with: %s fetch -y -f %s", rel, up, tool, rel);
      }
    }
    return st;
  }
  /* A directory? */
  size_t cnt, present, stale;
  dir_states(cfg, root, rel, &cnt, &present, &stale);
  if (!cnt) return YAME_STORE_NOT_CATALOGUED;
  if (stale) {
    if (advice) stale_dir_advice(cfg, root, rel, stale, cnt, advice, n);
    return YAME_STORE_STALE;
  }
  if (!present) {
    if (advice) snprintf(advice, n, "%s has not been fetched; run: %s fetch -y %s", rel, tool, rel);
    return YAME_STORE_ABSENT;
  }
  return YAME_STORE_CURRENT;
}

int yame_store_report(const yame_fetch_cfg_t *cfg, const char *root_override, FILE *out) {
  const char *tool = cfg->tool ? cfg->tool : "yame";
  char root[YAME_PATH_MAX];
  yame_assets_root(root_override, cfg->tool_env, root, sizeof root);
  int said = 0;
  for (size_t j = 0; j < cfg->n_files; ++j) {
    const yame_asset_file_t *f = &cfg->files[j];
    size_t dl = yame_file_dirlen(f);
    int seen = 0;
    for (size_t k = 0; k < j && !seen; ++k)
      if (yame_file_dirlen(&cfg->files[k]) == dl &&
          strncmp(cfg->files[k].store_path, f->store_path, dl) == 0) seen = 1;
    if (seen) continue;
    char sub[YAME_PATH_MAX];
    memcpy(sub, f->store_path, dl); sub[dl] = '\0';
    size_t cnt, present, stale;
    dir_states(cfg, root, sub, &cnt, &present, &stale);
    if (!stale) continue;
    char adv[1024];
    stale_dir_advice(cfg, root, sub, stale, cnt, adv, sizeof adv);
    fprintf(out, "[%s fetch] %s\n", tool, adv);
    ++said;
  }
  return said;
}

/* ------------------------------------------------ the tool's own registry */
static const yame_fetch_cfg_t *default_fetch_cfg;
void yame_set_default_fetch_cfg(const yame_fetch_cfg_t *cfg) { default_fetch_cfg = cfg; }
const yame_fetch_cfg_t *yame_default_fetch_cfg(void) { return default_fetch_cfg; }
