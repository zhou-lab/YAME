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
 *   $YAME_DATA_HOME/<platform-or-build>/...
 *
 *   The path mirrors the catalogue browser: hg38/data in the tree lands at
 *   <root>/hg38/data, while the upstream repository remains registry metadata.
 *   Every tool uses these same paths, so one tool's fetch satisfies the next.
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
#include <stdio.h>
#include <inttypes.h>

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
 * Prints at most once per process. Callers invoke it when they are about to
 * use the store; yame_assets_root() deliberately only resolves a path.
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

/**
 * Does this filename index another file in the same directory?
 *
 * Two spellings, because two ecosystems: YAME writes <x>.idx, tabix writes
 * <x>.tbi. One definition of it, because everything that treats a pair as one
 * thing has to agree -- the browser folds an index into its data file's row,
 * and a name lookup must not return one when asked for the other. They
 * disagreed once: "genes" resolved to genes.bed.gz.tbi, because .tbi sorts
 * after .bed.gz and only .idx was being skipped.
 *
 * Returns the suffix that matched (a static string), or NULL.
 */
const char *yame_assets_index_suffix(const char *name);

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
  YAME_PIN_ANCESTOR =  3,   /* filled at an earlier tag this build supersedes */
  YAME_PIN_CONFLICT = -1    /* populated at a tag this build does not pin */
};

/* A tag this build knows it supersedes: `anchor` identifies a store filled at
 * that tag, `tag` is what to call it when reporting the upgrade. The registry
 * emits one of these per past tag it has cached a manifest for. */
typedef struct {
  const char *tag;
  const char *anchor;
} yame_pin_prior_t;

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

/**
 * As pin_check, but able to tell an UPGRADE from a conflict.
 *
 * Plain pin_check is a scalar equality, so a store this build simply moved
 * past looks exactly like one some stranger filled. Given the anchors of the
 * tags this build knows it supersedes, a stored manifest matching one of them
 * is YAME_PIN_ANCESTOR: safe to overwrite without -f, because the two tags are
 * the same lineage and this build holds the later one.
 *
 * The re-download war pin_check exists to prevent is still prevented. An older
 * binary meeting a newer store sees a hash it has never heard of -- a newer tag
 * is not among its ancestors -- and still stops. Only a strictly known-earlier
 * tag is adopted, and only in the direction of the newer build.
 *
 * `prior` may be NULL (n_prior 0), which makes this identical to pin_check.
 */
int yame_assets_pin_state(const char *dir, const char *anchor_sha,
                          const yame_pin_prior_t *prior, size_t n_prior);

/* Which known earlier tag filled this directory, for the message that reports
 * the upgrade. NULL when no prior matches -- i.e. whenever pin_state did not
 * say YAME_PIN_ANCESTOR. */
const char *yame_assets_pin_prior_tag(const char *dir,
                                      const yame_pin_prior_t *prior,
                                      size_t n_prior);

/* ---------------------------------------------------- the registry's rows */
/* One file a directory publishes, with the digest it must have.
 *
 * store_sub is normally the unit's, and NULL says so. It is set only where a
 * file belongs somewhere else in the store than the directory that publishes
 * it: a genome's cpg_nocontig.cr comes from the KYCGKB repo but is the
 * genome's index, so it lands at <genome>/ rather than <genome>/KYCG/. The
 * browser renders a file wherever this puts it, so the tree and the store
 * cannot drift apart. */
typedef struct {
    const char *name;
    const char *sha256;
    uint64_t    size;        /* 0 when upstream does not publish one */
    const char *store_sub;   /* NULL: the unit's own store_sub */
} yame_asset_file_t;

#define YAME_NFILES(t) (sizeof(t)/sizeof((t)[0]) - 1)
#define YAME_NPRIOR(t) (sizeof(t)/sizeof((t)[0]))

/* One fetchable directory: where it comes from, where it lands, and the
 * digest its manifest must have. A tool's registry.h is an array of these;
 * the types live HERE, not in that generated file, so library code can take
 * any tool's registry as an argument rather than being bound to yame's at
 * yame's build time. */
typedef struct {
    const char *source;      /* upstream repo family: InfiniumAnnotation, ... */
    const char *target;      /* platform, genome, or "<platform>/KYCG" */
    const char *base_url;    /* <base>/<tag>/<remote_sub>/SHA256SUMS */
    const char *tag;
    const char *remote_sub;  /* "" when the manifest is at the repo root */
    const char *store_sub;   /* path under the store root */
    const char *anchor;      /* sha256 of that directory's SHA256SUMS */
    const yame_asset_file_t *files;  /* what the directory holds */
    size_t      n_files;
    const yame_pin_prior_t *prior;   /* earlier tags this build supersedes */
    size_t      n_prior;
} yame_asset_reg_t;

/* ------------------------------------------ is the store what this build pins? */
/**
 * What a tool should say about a store directory, or a file in one, before it
 * trusts it. This is the one call a downstream tool makes at load time -- and
 * what `fetch` itself prints when a store is behind -- so that every tool
 * describes the same situation in the same words.
 *
 * The registry is an ARGUMENT, not yame's compiled-in one: methscope passes
 * its own, so the verb in the advice is `methscope fetch`, and the tag it is
 * judged against is the one methscope pins.
 */
typedef enum {
  YAME_STORE_CURRENT = 0,   /* at this registry's tag; a named file's digest matches  */
  YAME_STORE_ABSENT,        /* nothing fetched for this directory yet                  */
  YAME_STORE_OLD_TAG,       /* filled at an earlier tag this registry knows -- fetch it */
  YAME_STORE_OTHER_TAG,     /* filled at a tag this registry does not know -- the TOOL
                             * is behind the store, or another tool filled it          */
  YAME_STORE_UNPINNED,      /* a manifest is there but this registry has no anchor     */
  YAME_STORE_STALE_FILE,    /* directory is current but this file's digest differs     */
  YAME_STORE_MISSING_FILE,  /* directory is current but this file is not on disk       */
  YAME_STORE_NOT_CATALOGUED /* the path is under no directory this registry lists      */
} yame_store_state_t;

/**
 * Classify `path` -- a store directory, or a file inside one -- against
 * `reg`. `path` may be absolute or relative to the store root that
 * yame_assets_root() resolves from `tool_env` and $YAME_DATA_HOME.
 *
 * `advice`, if given, receives one sentence naming the situation and the
 * exact command that fixes it, spelled with `tool` ("yame", "methscope"):
 *
 *   OLD_TAG    "hg38/models was fetched at an earlier tag than this yame;
 *               run: yame fetch hg38/models"
 *   OTHER_TAG  "hg38/models is at a tag this yame does not know -- update
 *               yame, then run: yame fetch hg38/models"
 *
 * CURRENT leaves `advice` empty. Nothing here downloads or deletes.
 */
yame_store_state_t yame_store_state(const yame_asset_reg_t *reg, size_t n_reg,
                                    const char *tool, const char *tool_env,
                                    const char *path, char *advice, size_t n);

/**
 * Walk every directory in `reg` and print one advice line to `out` for each
 * that is OLD_TAG or OTHER_TAG. Returns how many lines it printed, so a
 * caller can decide whether to say anything more. This is what a bare
 * `<tool> fetch` and `<tool> fetch -l` print on stderr.
 */
int yame_store_report(const yame_asset_reg_t *reg, size_t n_reg,
                      const char *tool, const char *tool_env,
                      const char *root_override, FILE *out);

/* --------------------------------------------------------------- fetching */

typedef struct {
  int force;                 /* re-download what is present; overrule a pin */
  int quiet;                 /* no progress reporting */

  /* The tags this build supersedes, for the directory being fetched. Left
   * NULL, a fetch keeps the old all-or-nothing behaviour: any manifest that
   * is not the pinned one is a conflict needing -f. Set, an upgrade from one
   * of these proceeds on its own. It rides here rather than in the argument
   * list because it is per-directory registry data, like the anchor it
   * qualifies, and every caller already carries an options struct. */
  const yame_pin_prior_t *prior;
  size_t n_prior;

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

/* ------------------------------------------- inferring a row-space reference */

/**
 * Which reference a file of `rows` rows is written against.
 *
 * The row count is the only handle a CX record gives on what its rows ARE,
 * and it is enough: no two row spaces in the catalogue share a count. This is
 * what lets a command work out -R rather than make the caller repeat what the
 * file already implies.
 *
 * `want_kind` is "genome" for a coordinate stream (.cr) or "array" for a probe
 * ordering; NULL accepts either. On YAME_REF_OK, `path` holds the resolved
 * file. `hit` (optional) receives the matched row, which carries the name and
 * the `yame fetch` argument for an error message.
 *
 * Resolves only -- it never downloads. A command that finds the reference
 * missing should say which fetch would supply it and stop, not start a
 * multi-megabyte transfer in the middle of a print.
 */
enum {
  YAME_REF_OK         =  0,   /* identified and in the store */
  YAME_REF_MISSING    =  1,   /* identified, not downloaded */
  YAME_REF_WRONG_KIND =  2,   /* identified, but not the kind asked for */
  YAME_REF_NO_NAME    =  3,   /* row space known, nothing in it by that name */
  YAME_REF_UNKNOWN    = -1    /* the count matches nothing known */
};

int yame_ref_for_rows(uint64_t rows, const char *store_override,
                      const char *want_kind, char *path, size_t n,
                      const char **name, const char **fetch);

/**
 * Turn a -R / -m argument into a path.
 *
 * An existing path is used as given -- the ordinary spelling costs nothing
 * and cannot change meaning. Anything else is a NAME, resolved against the
 * store for the row space `rows` identifies: the row space's own name
 * ("hg38") gives its reference, and any other name gives the newest file of
 * that name in any directory the row space owns -- its own first, then its
 * knowledgebase -- so `-m ChromHMM` finds ChromHMM.20220303.cm without anyone
 * having to know the date or the directory.
 */
int yame_ref_resolve(const char *spec, uint64_t rows, const char *store_override,
                     const char *want_kind, char *path, size_t n,
                     const char **name, const char **fetch);

/**
 * Resolve a spec that may name more than one file.
 *
 * For a repeatable `-m`, where one argument standing for several sets saves
 * the caller from spelling out a directory listing. Everything
 * yame_ref_resolve() accepts still means the same thing here and yields one
 * path; the additions are:
 *
 *   CGI,ChromHMM       a comma list -- each element resolved on its own, in
 *                      the order written, duplicates dropped
 *   mm10:ChromHMM      resolved in the row space NAMED here rather than the
 *                      one `rows` implies. This is how a caller reaches a set
 *                      for a row space it does not have a file from, and the
 *                      only form that works when `rows` is 0
 *   mm10:  /  mm10:*   every set in that row space
 *   *                  every set in the row space `rows` implies
 *   mm10               likewise every set -- a bare row space name means its
 *                      sets HERE, not its reference. Resolving it to a .cr
 *                      would hand back a coordinate stream where a mask was
 *                      asked for, which is never what the caller wanted.
 *                      yame_ref_resolve() keeps the reference meaning, so
 *                      `-R mm10` is unaffected.
 *
 * "Every set" is the .cm files a row space owns, its own directory before its
 * knowledgebase, newest per set name, first directory winning -- the same
 * precedence a single name follows.
 *
 * On YAME_REF_OK, `*paths` is a malloc'd array of `*n_paths` malloc'd strings;
 * release both with yame_ref_paths_free(). Any element failing to resolve
 * fails the call, so a caller never silently gets a subset of what it asked
 * for; `name`/`fetch` then describe the element that failed.
 */
int yame_ref_resolve_multi(const char *spec, uint64_t rows,
                           const char *store_override, const char *want_kind,
                           char ***paths, size_t *n_paths,
                           const char **name, const char **fetch);

/** Release what yame_ref_resolve_multi() returned. */
void yame_ref_paths_free(char **paths, size_t n);

/** Rows in a CX file's first record, or 0 if it cannot be read. */
uint64_t yame_ref_file_rows(const char *path);

/** Print why a name did not resolve. */
void yame_ref_explain_name(FILE *out, const char *spec, uint64_t rows,
                           int status, const char *name, const char *fetch,
                           const char *flag);

/**
 * Print why an inference did not produce a usable reference: which reference
 * the row count named, and the `yame fetch` that would supply it. `name` and
 * `fetch` are what yame_ref_for_rows returned; either may be NULL.
 */
void yame_ref_explain(FILE *out, uint64_t rows, int status, const char *name,
                      const char *fetch, const char *flag);

/* ------------------------------------------------------ choosing in a tree */

/**
 * Browse the catalogue and hand back what was chosen, fetching it if needed.
 *
 * The `yame fetch` screen, opened at `open_unit` (a row space name such as
 * "hg38", or NULL), with one extra verb: `u` ends the session and returns the
 * selection. It lets a command offer "show me what there is" instead of
 * requiring a path to a file the user has not downloaded yet.
 *
 * Returns the number of paths, malloc'd into *paths (caller frees the strings
 * and the array), or 0 if nothing was chosen or the terminal cannot host a
 * tree -- in which case the caller should say what it needs and stop.
 */
/* ------------------------------------- fetch, over the caller's registry */
/**
 * What `fetch` needs to know about the tool it is running inside. yame's own
 * entry point builds one of these from its compiled-in registry; a downstream
 * tool that ships fetch over the YAME code it bundles builds one from ITS
 * registry (make_registry.sh --tool=<name>), so the models its docs name are,
 * by construction, the ones its own binary pins. Store conflicts between
 * tools stay handled by yame_assets_pin_state.
 */
typedef struct {
  const yame_asset_reg_t *reg;  /* the catalogue: one row per directory       */
  size_t n_reg;
  const char *tool;             /* "yame", "methscope": usage text, messages   */
  const char *tool_env;         /* NULL, or e.g. "METHSCOPE_DATA_HOME" -- read
                                 * ahead of $YAME_DATA_HOME for the store root */
} yame_fetch_cfg_t;

/* `<tool> fetch ...` -- the whole subcommand: browser, -l, named targets,
 * the single-file form, mirrors, -t/-k. argv[0] is ignored. */
int yame_fetch_main(const yame_fetch_cfg_t *cfg, int argc, char *argv[]);

/* The mask picker `summary -b` uses, over the same registry. */
size_t yame_browse_pick(const yame_fetch_cfg_t *cfg, const char *open_unit,
                        char ***paths);

#endif /* _YAME_ASSETS_H */
