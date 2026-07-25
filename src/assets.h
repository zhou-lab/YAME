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
 * Shared asset store: where downloaded reference data lives, and how it gets
 * there.
 *
 * The CLIs that link this library each grew their own copy of the same
 * downloader, and each kept its assets in a private cache, so two tools
 * consuming the identical upstream file downloaded it twice. This is the one
 * implementation they all link, and the one store they all read.
 *
 * THE STORE
 *   $YAME_DATA_HOME/<source>/<platform-or-build>/...
 *
 *   The key convention is the reuse mechanism: a file fetched by one tool is
 *   found by the next because the path is derived from the upstream repo and
 *   the platform, not from who asked for it.
 *
 * THE TRUST CHAIN
 *   Each directory holds a SHA256SUMS that is a byte-identical copy of the
 *   remote one, and the caller supplies an anchor -- sha256(SHA256SUMS) --
 *   compiled into its binary. Fetch verifies the manifest against the anchor,
 *   then every file against a digest from that manifest. Nothing trusts a
 *   digest fetched at run time, and `shasum -a 256 -c SHA256SUMS` re-verifies
 *   a store by hand with none of this code.
 *
 * THE PIN
 *   That same stored manifest is how a directory says which upstream tag filled
 *   it: its hash IS the anchor. Two tools pinned to different tags therefore
 *   cannot silently overwrite each other -- see yame_assets_pin_check().
 */

#ifndef _YAME_ASSETS_H
#define _YAME_ASSETS_H

#include <stddef.h>
#include <stdint.h>

/* The manifest filename, in every directory of every source. */
#define YAME_ASSETS_SUMS_FILE "SHA256SUMS"

/* ------------------------------------------------------------------ digest */
/* SHA-256 (FIPS 180-4), self-contained; see digest.c. */

void yame_assets_sha256_buf(const void *data, size_t len, char out[65]);
int  yame_assets_sha256_file(const char *path, char out[65]);

/* Case-insensitive compare of two hex digests. Returns 1 when equal. */
int  yame_assets_digest_equal(const char *a, const char *b);

/* ------------------------------------------------------------- store paths */

/**
 * Resolve the store root, in this order:
 *
 *   override            (a --store/-d flag; NULL or "" to skip)
 *   $<tool_env>         (KYCG_DATA_DIR, SESAME_INDEX_DIR, ...; NULL to skip)
 *   $YAME_DATA_HOME
 *   ${XDG_DATA_HOME:-$HOME/.local/share}/yame
 *
 * The default is the data tier rather than a cache tier on purpose: these are
 * multi-GB reference sets that should survive somebody clearing ~/.cache.
 *
 * Writes into buf and returns it (never NULL). Asking where the store is says
 * nothing to the user: a caller that is about to USE it calls
 * yame_assets_legacy_notice() itself, so `--help` stays quiet.
 */
const char *yame_assets_root(const char *override, const char *tool_env,
                             char *buf, size_t n);

/**
 * If the resolved root does not exist yet but a pre-consolidation per-tool
 * cache does (the ~/.cache/<tool> directories used before this store existed),
 * print one line to stderr saying so. That cache is never read from and never
 * moved: this exists only so a user whose downloads appear to have vanished is
 * told where they went.
 *
 * Prints at most once per process. Safe to call directly; yame_assets_root()
 * already does.
 */
void yame_assets_legacy_notice(const char *root);

/* Is the root writable (or creatable)? A read-only shared store is legitimate,
 * so this is not an error by itself -- it is what lets fetch say "this store is
 * read-only" up front rather than failing partway through a transfer. */
int yame_assets_root_writable(const char *root);

int yame_assets_join(char *out, size_t n, const char *a, const char *b);
int yame_assets_mkdir_p(const char *path);
int yame_assets_mkdir_parent(const char *path);   /* mkdir -p $(dirname path) */
int yame_assets_is_file(const char *path);
int yame_assets_is_dir(const char *path);

/**
 * Reject a name that must not become a path component: empty, leading '.'
 * (so "." and ".." and hidden files), any '/' or '\\', any control character.
 * Every name read out of a SHA256SUMS goes through this before it is joined to
 * a directory -- a manifest entry is remote input, and "../../x" would
 * otherwise write outside the store.
 *
 * Nested names ARE legitimate in a manifest (e.g. "KYCG/foo.cm"), so callers
 * that support them validate each component; see yame_assets_safe_relpath().
 */
int yame_assets_safe_name(const char *s);

/* Like safe_name but for a relative path: every '/'-separated component must
 * itself be a safe name. Returns 1 when the whole path is safe. */
int yame_assets_safe_relpath(const char *s);

/* Directory holding the running executable, for a tool that ships data beside
 * its binary, e.g. <exe>/data. Returns 0 on success. */
int yame_assets_exe_dir(char *out, size_t n);

/* --------------------------------------------------------------- manifests */

typedef struct {
  char sha[65];
  char name[512];        /* may contain '/' -- validated as a relative path */
} yame_sums_ent_t;

/* Parse SHA256SUMS text. Entries whose name is unsafe are dropped. Caller
 * frees the returned array. */
yame_sums_ent_t *yame_assets_parse_sums(const char *text, size_t *n);
yame_sums_ent_t *yame_assets_sums_load_file(const char *path, size_t *n);

/* ----------------------------------------------------------------- the pin */

enum {
  YAME_PIN_MATCH    =  0,   /* stored manifest hashes to the caller's anchor */
  YAME_PIN_ABSENT   =  1,   /* no manifest stored yet -- nothing to conflict */
  YAME_PIN_UNKNOWN  =  2,   /* manifest present but caller passed no anchor */
  YAME_PIN_CONFLICT = -1    /* populated at a tag this build does not pin */
};

/**
 * Compare <dir>/SHA256SUMS against the caller's compiled anchor.
 *
 * A YAME_PIN_CONFLICT means another tool (or an older build of this one)
 * populated the directory from a different upstream tag. Fetching would
 * overwrite it, and then that other tool would overwrite it back on its next
 * run, forever. So a conflict stops a write; a reader may proceed with a
 * warning, since reading never re-verifies digests anyway.
 */
int yame_assets_pin_check(const char *dir, const char *anchor_sha);

/* --------------------------------------------------------------- fetching */

typedef struct {
  int force;                 /* re-download what is present; overrule a pin */
  int quiet;                 /* no progress reporting */

  /* Optional progress hooks. Left NULL, a fetch is silent apart from errors,
   * which is what a library caller with its own UI wants: a caller passes its
   * own renderer here rather than inheriting one. */
  void (*on_begin)(void *ud, const char *name, uint64_t total);
  void (*on_progress)(void *ud, uint64_t now, uint64_t total);
  void (*on_done)(void *ud, const char *name, uint64_t bytes, int ok);
  void *ud;
} yame_fetch_opt_t;

/* Was this build compiled against libcurl? Everything below fails cleanly with
 * a "built without libcurl" error when it was not. */
int yame_assets_have_curl(void);

/* GET into memory, NUL-terminated. Caller frees. NULL on failure. */
char *yame_assets_http_get_mem(const char *url, size_t *len);

/**
 * Download one URL to `dest`, verifying it against `want_sha` before it is
 * allowed to exist under that name: the transfer lands in a per-process
 * "<dest>.<pid>.part" opened O_EXCL, is hashed there, and is renamed into
 * place only if the digest matches. A mismatch discards it.
 *
 * With *downloaded set, reports whether bytes actually moved (0 means the file
 * was already present and verified). On failure returns non-zero and, if err
 * is non-NULL, stores a malloc'd message there.
 */
int yame_assets_download_verify(const char *url, const char *want_sha,
                                const char *dest, const yame_fetch_opt_t *opt,
                                int *downloaded, char **err);

/**
 * Fetch a whole SHA256SUMS-anchored directory:
 *
 *   <base>/<tag>/<remote_sub>/{SHA256SUMS, files...}  ->  <store_sub>/
 *
 * The manifest is pulled first and checked against `anchor_sha` (pass NULL to
 * accept whatever the remote publishes -- only sensible for a tag this build
 * does not pin). Then every file it lists is fetched and verified, and the
 * manifest is written into the store verbatim so the result re-verifies with
 * shasum alone.
 *
 * Refuses to touch a directory whose stored manifest conflicts with
 * `anchor_sha` unless opt->force; see yame_assets_pin_check().
 */
int yame_assets_fetch_subtree(const char *base, const char *tag,
                              const char *remote_sub, const char *store_sub,
                              const char *anchor_sha,
                              const yame_fetch_opt_t *opt, char **err);

/**
 * As above, but only the entries named in `only` (n_only of them). Passing
 * only == NULL fetches everything, which is exactly what fetch_subtree does.
 *
 * A partial directory is a normal state, not a degraded one: a knowledgebase
 * directory can hold dozens of sets and most callers want a few. The manifest
 * is still written verbatim, because it describes the TAG rather than what was
 * taken from it -- that is what keeps the pin check meaningful, and
 * `shasum -a 256 -c SHA256SUMS` then reports the ones not taken as missing,
 * which is the truth.
 */
int yame_assets_fetch_subset(const char *base, const char *tag,
                             const char *remote_sub, const char *store_sub,
                             const char *anchor_sha,
                             const char *const *only, size_t n_only,
                             const yame_fetch_opt_t *opt, char **err);

#endif /* _YAME_ASSETS_H */
