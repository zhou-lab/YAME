// SPDX-License-Identifier: LicenseRef-CHOP-Academic-BSD-2-Clause
// 2-Clause BSD for academic and non-profit research use; commercial use by
// inquiry to zhouw3@chop.edu (see LICENSE).
/**
 * This file is part of YAME.
 *
 * Copyright (C) 2021-present The Children's Hospital of Philadelphia
 *
 * Use of this software is available to academic and non-profit institutions
 * for research purposes under the 2-Clause BSD License; for use or transfers
 * to commercial entities, inquire with Dr. Wanding Zhou at zhouw3@chop.edu.
 * See the LICENSE file at the root of the repository for the full terms.
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
 *   $<each yame_assets_suite_env>   (METHSCOPE_DATA_HOME, ...)
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
 * The store variables the suite's other tools export, NULL-terminated. The
 * store is one shared tree, so whichever tool a reader used to move it moved
 * it for all of them -- `yame hprint` under a store set with
 * $METHSCOPE_DATA_HOME must find the coordinate track that is sitting in it.
 * Read after $YAME_DATA_HOME, so yame's own variable still wins for yame.
 */
extern const char *const yame_assets_suite_env[];

/**
 * The name of the environment variable that resolved the store, for printing
 * beside the path, or NULL when nothing in the environment did (the XDG or
 * ~/.local/share default). Same order as yame_assets_root(), minus `override`:
 * a caller that passed -d already knows it won.
 */
const char *yame_assets_root_env(const char *tool_env);

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

/* ---------------------------------------------------- the registry's rows */
/**
 * One file the suite can fetch. A tool's registry.h is an array of these and
 * nothing else -- no directory rows, no anchors, no ancestry. The type lives
 * HERE, not in the generated file, so library code can take any tool's
 * registry as an argument rather than being bound to yame's at build time.
 *
 * `key` is the file's identity, source@tag:remote_path: which repo, which
 * immutable tag, which path within it. The repo and the tag are read out of
 * it with yame_key_*() below. `store_path` is where the file lands under the
 * store root and its place in the browser tree; it is decoupled from the
 * remote path. `url` is compiled rather than derived so the host rule exists
 * in one place (tools/registry/lib.sh). The four prose fields are what the
 * browser shows; `source` here is the upstream database, not the repo.
 */
typedef struct {
    const char *key;
    const char *store_path;
    const char *url;
    const char *sha256;
    uint64_t    size;        /* 0 when upstream does not publish one */
    const char *title;
    const char *description;
    const char *source;
    const char *citation;
} yame_asset_file_t;

/* The array MUST end with a zeroed row; YAME_NFILES subtracts it. An emitter
 * that omits the terminator does not fail to compile or to run -- it silently
 * reports one file too few, and the file that vanishes is whichever the
 * table's order put last. */
#define YAME_NFILES(t) (sizeof(t)/sizeof((t)[0]) - 1)

/* The parts of a key, split on the LAST @ and the LAST colon: a source may
 * carry a colon (hf:org/repo), a path never does. Each writes at most n-1
 * bytes plus NUL and returns the length it wanted. */
size_t yame_key_source(const char *key, char *out, size_t n);
size_t yame_key_tag(const char *key, char *out, size_t n);
size_t yame_key_path(const char *key, char *out, size_t n);

/* The file name within its store directory, and that directory (as a
 * length: everything before the last slash, 0 for the root). */
const char *yame_file_name(const yame_asset_file_t *f);
size_t yame_file_dirlen(const yame_asset_file_t *f);

/* ------------------------------------- fetch, over the caller's registry */
/**
 * What `fetch` needs to know about the tool it is running inside. yame's own
 * entry point builds one of these from its compiled-in registry; a downstream
 * tool that ships fetch over the YAME code it bundles builds one from ITS
 * registry, which it generates itself from the shared table, so the models
 * its docs name are, by construction, the ones its own binary pins.
 */
typedef struct {
  const yame_asset_file_t *files; /* the registry: one row per file          */
  size_t n_files;
  const char *tool;             /* "yame", "methscope": usage text, messages   */
  const char *tool_env;         /* NULL, or e.g. "METHSCOPE_DATA_HOME" -- read
                                 * ahead of $YAME_DATA_HOME for the store root */
  int no_prompt;                /* never ask on a terminal: yame_store_resolve
                                 * refuses an ambiguous name instead of showing
                                 * a picker. Affects yame_store_resolve ONLY:
                                 * the fetch browser and its confirmations are
                                 * unaffected, so one cfg serves both. A tool whose commands must not turn
                                 * interactive (kycg annotate) sets it once.
                                 * `fetch` itself never prompts over a name:
                                 * one claimed by several directories is an
                                 * error listing them, -y or not. Zero in a
                                 * four-field initializer, so existing callers
                                 * are unchanged. */
} yame_fetch_cfg_t;

/* ------------------------------------------ is the store what this build pins? */
/**
 * A file's state in the store, judged by the line that names it in the
 * directory's SHA256SUMS -- the record of what was downloaded -- against the
 * digest this registry compiles in. No file is hashed here: a 3 GB model is
 * classified by one line of text.
 *
 *   CURRENT   on disk, and the manifest records the digest this build pins
 *             (or records nothing about it -- a file that was put there by
 *             hand is taken at its word until a fetch verifies it)
 *   ABSENT    not on disk
 *   STALE     on disk, but the manifest records a different digest: fetched
 *             by an earlier or later build, or by another tool at another
 *             tag. A fetch replaces it; -f is asked for first.
 *
 * There is no directory-level state: a directory is the sum of its files.
 */
typedef enum {
  YAME_STORE_CURRENT = 0,
  YAME_STORE_ABSENT,
  YAME_STORE_STALE,
  YAME_STORE_NOT_CATALOGUED   /* the path names nothing this registry lists */
} yame_store_state_t;

yame_store_state_t yame_file_state(const char *root, const yame_asset_file_t *f);

/**
 * Classify `path` -- a store path, absolute or relative to the store root
 * resolved from `tool_env` and $YAME_DATA_HOME -- against `cfg`. A directory
 * is STALE if any file in it is, ABSENT if none of them is on disk, else
 * CURRENT. `advice`, if given, receives one sentence naming the situation and
 * the exact command that fixes it, spelled with cfg->tool and carrying -y so
 * it also works where nobody can answer a prompt.
 */
yame_store_state_t yame_store_state(const yame_fetch_cfg_t *cfg, const char *path,
                                    char *advice, size_t n);

/**
 * Walk every directory the registry knows and print one line to `out` for
 * each that holds stale files:
 *
 *   [methscope fetch] hg38/models: 3 of 7 files differ from this build;
 *                     run: methscope fetch -y -f hg38/models
 *
 * Returns how many lines it printed. This is what a bare `<tool> fetch` and
 * `<tool> fetch -l` print on stderr.
 */
int yame_store_report(const yame_fetch_cfg_t *cfg, const char *root_override,
                      FILE *out);

/**
 * Turn a command-line argument into a store file, the one way for every tool
 * in the suite -- a model for methscope, a platform's ordering for sesame, a
 * knowledgebase set for kycg.
 *
 * An existing path is used as given: `path` is that spelling, `*rec` is NULL,
 * the result is CURRENT -- the ordinary case, and it costs nothing. Anything
 * else is a NAME looked up in cfg->files: a store path ("hg38/models/x.clfx")
 * exactly; or a file name ("x.clfx") or SET name -- the part before the first
 * dot, case-insensitive ("CGI" for CGI.20220904.cm, "hg38_celltype_lite" for
 * the .clfx), the shorthand `-m` takes -- with an index never answering for
 * its data file; or a glob over store paths and file names (`hg38/models/` followed by `*`).
 * Reaching ONE file, `path` is <store root>/<store_path>, `*rec` the record,
 * and the result is that file's state: CURRENT (use it), ABSENT (not
 * fetched) or STALE (on disk at a digest this build does not pin); for those
 * two, `advice` carries the sentence to print, with the `<tool> fetch`
 * command that repairs it. Reaching SEVERAL, nothing is guessed -- not the
 * newest, not the first: on a terminal the person is shown the candidates
 * and picks one; off a terminal, or with cfg->no_prompt set, the call
 * refuses, `path` empty and `advice` naming every candidate, so a script
 * must spell the file. The ABSENT advice for one file names it without -y
 * (a single file is its own confirmation); a directory's carries -y. NOT_CATALOGUED
 * also for no such name, with `path` the spec unchanged so the caller's
 * opener can report it in its own words. The state is returned, not
 * enforced: whether to stop on STALE or warn and go on is the caller's
 * policy (the suite's: warn and go on; stop on ABSENT). Resolves only;
 * never downloads.
 */
yame_store_state_t yame_store_resolve(const yame_fetch_cfg_t *cfg, const char *spec,
                                      const char *root_override,
                                      char *path, size_t n,
                                      const yame_asset_file_t **rec,
                                      char *advice, size_t adv_n);

/**
 * Where several files are allowed: everything a spec names, expanded. An
 * existing path gives itself (one entry, its record NULL); a name or glob
 * gives every file it reaches, in registry order, each with its store path
 * and record. Returns 0 with `*paths`/`*recs` malloc'd arrays of `*n`
 * entries -- release the paths with yame_ref_paths_free() and free(*recs)
 * -- or -1 when nothing is called that, `advice` saying so. State per file
 * is yame_file_state(); nothing is asked, nothing downloaded.
 */
int yame_store_resolve_multi(const yame_fetch_cfg_t *cfg, const char *spec,
                             const char *root_override,
                             char ***paths, const yame_asset_file_t ***recs,
                             size_t *n, char *advice, size_t adv_n);

/* The stale files themselves, in registry order, so a caller can offer to
 * replace exactly those: `yame fetch` asks before its browser opens. Returns
 * how many there are, filling up to `cap`. */
size_t yame_store_stale(const yame_fetch_cfg_t *cfg, const char *root_override,
                        const yame_asset_file_t **out, size_t cap);

/* --------------------------------------------------------------- fetching */

typedef struct {
  int force;                 /* re-download what is present; replace what is stale */
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
 * Fetch the files in `want` (n_want of them, rows of cfg->files) into the
 * store under `root`, each verified against the digest the registry compiles
 * in, and then rewrite the SHA256SUMS of every directory touched.
 *
 * A file whose manifest line already records the pinned digest is skipped
 * without being hashed; one with no line is hashed and skipped if it
 * matches; anything else is downloaded to a .part and renamed into place
 * only when its digest matches. So a one-file bump costs one file, and -f
 * (opt->force) is what makes a present file move again.
 *
 * The manifest a directory gets lists every registry file in that directory,
 * in registry order, in sha256sum format: the pinned digest for each file
 * this call verified, the line the old manifest had for each file it did
 * not touch, and the pinned digest for files not on disk at all -- so the
 * manifest never claims a file is current that this call did not confirm,
 * and a directory that is an exact copy of an upstream directory ends up
 * byte-identical to upstream's own SHA256SUMS. No upstream manifest is
 * downloaded: the digests the binary carries are the whole trust chain.
 */
int yame_assets_fetch_files(const yame_fetch_cfg_t *cfg, const char *root,
                            const yame_asset_file_t *const *want, size_t n_want,
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
/* `<tool> fetch ...` -- the whole subcommand: browser, -l, named targets,
 * the single-file form, mirrors, -t/-k. argv[0] is ignored. */
int yame_fetch_main(const yame_fetch_cfg_t *cfg, int argc, char *argv[]);

/* The registry the running tool considers its own, for library code that
 * needs one but is not handed one -- `summary -b` picking masks through the
 * browser is the case today. The executable registers it at startup; a
 * library object must never reference a symbol only the executable defines,
 * because a static archive links whole objects and every downstream that
 * pulls summary.o for summarize1() would then fail to link. (v1.43's first
 * cut did exactly that.) NULL until registered: a downstream tool that has
 * not registered one gets "no registry in this build", not a link error. */
void yame_set_default_fetch_cfg(const yame_fetch_cfg_t *cfg);
const yame_fetch_cfg_t *yame_default_fetch_cfg(void);

/* The mask picker `summary -b` uses, over the same registry. */
size_t yame_browse_pick(const yame_fetch_cfg_t *cfg, const char *open_unit,
                        char ***paths);

#endif /* _YAME_ASSETS_H */
