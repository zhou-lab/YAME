// SPDX-License-Identifier: LicenseRef-CHOP-Academic-BSD-2-Clause
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
 * `yame fetch` -- populate the shared asset store.
 *
 * A driver over src/assets.c, so the store can be filled without any
 * particular downstream tool installed, and so the engine they all link has a
 * way to be exercised on its own.
 *
 * Two forms:
 *   yame fetch <name>          a store directory, or one file in it
 *   yame fetch -u URL -s SHA -o DEST     one file, digest given by hand
 *
 * The second is the low-level form: a per-file digest supplied by the caller,
 * which is how a store gets a file this build's registry has never heard of.
 */

/* strcasestr is a GNU extension: -g matches case-insensitively, and this is
 * the one call that needs the feature macro. */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fnmatch.h>
#include <unistd.h>

#include "assets.h"
/* No registry.h here. fetch runs over whatever registry the caller hands it
 * through yame_fetch_cfg_t; the two macros below make the ~40 sites that read
 * the catalogue read the caller's. The pointer is file-scope because every
 * one of those sites is a static helper in this file, reached only from
 * yame_fetch_main() or yame_browse_pick(), which set it first. */
static const yame_fetch_cfg_t *cfg_;
#define TOOL          (cfg_->tool)
#include "yame_ui.h"

#include "wzmisc.h"   /* wzmalloc/wzstrdup: allocation that cannot return NULL */

/* ---- the registry, read as directories ----
 *
 * The compiled registry is a flat list of files. The browser and every
 * command here think in store DIRECTORIES -- hg38, EPIC/KYCG, hg38/models --
 * because that is what a person fetches, verifies and is told about. So the
 * list is grouped once, at entry, into this view: one unit per distinct
 * directory, its files in table order, and the one source@tag every file in
 * it carries (the generator refuses a directory that mixes two). Nothing is
 * compiled per directory; this is derived, and it dies with the process. */
typedef struct {
  char   dir[256];                 /* the store directory; the browser path */
  char   source[160];              /* the repo, out of the keys */
  char   tag[64];
  const yame_asset_file_t **files; /* into cfg_->files, table order */
  size_t n_files;
} unit_t;
static unit_t *UNITS;
static size_t  N_UNITS;
#define YAME_ASSETS   (UNITS)
#define YAME_ASSETS_N (N_UNITS)

static void build_units(void) {
  /* Rebuilt on every entry: a downstream tool may hand fetch a different
   * registry than last time (yame_browse_pick after yame_fetch_main). */
  for (size_t i = 0; i < N_UNITS; ++i) free(UNITS[i].files);
  free(UNITS); UNITS = NULL; N_UNITS = 0;
  UNITS = wzcalloc(cfg_->n_files ? cfg_->n_files : 1, sizeof(unit_t));
  for (size_t j = 0; j < cfg_->n_files; ++j) {
    const yame_asset_file_t *f = &cfg_->files[j];
    size_t dl = yame_file_dirlen(f);
    size_t u = 0;
    for (; u < N_UNITS; ++u)
      if (strlen(UNITS[u].dir) == dl && strncmp(UNITS[u].dir, f->store_path, dl) == 0) break;
    if (u == N_UNITS) {
      unit_t *nu = &UNITS[N_UNITS++];
      if (dl >= sizeof nu->dir) dl = sizeof nu->dir - 1;
      memcpy(nu->dir, f->store_path, dl); nu->dir[dl] = '\0';
      yame_key_source(f->key, nu->source, sizeof nu->source);
      yame_key_tag(f->key, nu->tag, sizeof nu->tag);
      nu->files = wzcalloc(cfg_->n_files, sizeof(*nu->files));
    }
    UNITS[u].files[UNITS[u].n_files++] = f;
  }
}

static int usage(void) {
  char root[4096];
  yame_assets_root(NULL, cfg_->tool_env, root, sizeof(root));

  { char head[128];
    snprintf(head, sizeof head, "%s fetch                              browse the catalogue", TOOL);
    yame_usage_head(head); }
  { char l1[128], l2[128];
    snprintf(l1, sizeof l1, "%s fetch [options] <name> ...", TOOL);
    snprintf(l2, sizeof l2, "%s fetch [options] -u <url> -s <sha256> -o <dest>", TOOL);
    yame_usage_text(l1); yame_usage_text(l2); }

  yame_usage_sec("Naming:");
  yame_usage_text("A name is a store directory, as the browser shows it: hg38,");
  yame_usage_text("hg38/KYCG, hg38/data, EPIC. It takes that directory's own files --");
  yame_usage_text("`hg38` is the genome annotation, not the knowledgebase and models");
  yame_usage_text("beneath it; name those explicitly. Narrow within one with -g.");
  yame_usage_text("A file resolves too, best written out: `hg38/data/test.cg`. The");
  yame_usage_text("bare name works when only one directory publishes it.");
  yame_usage_text("Several names may be given, separated by spaces or by commas --");
  yame_usage_text("commas so a list fits an option that takes one argument. A");
  yame_usage_text("directory named twice is taken once, and naming it whole absorbs a");
  yame_usage_text("file picked out of it.");

  yame_usage_sec("Browsing:");
  yame_usage_text("With no target on a terminal, opens a tree browser: species, then");
  yame_usage_text("platform or genome build, then its knowledgebase and files. Arrows");
  yame_usage_text("move, right/left open and close a row, space or x selects (a folder");
  yame_usage_text("takes everything under it), f fetches what is selected, h lists every");
  yame_usage_text("key, q leaves. `/` filters the tree -- by name, source, collection");
  yame_usage_text("or title -- and enter keeps the filter so you can then select.");
  yame_usage_text("What is already in the store shows as present and");
  yame_usage_text("cannot be selected.");
  yame_usage_text("Piped or redirected it dumps the registry as TSV instead (same as");
  yame_usage_text("-l), so a script never blocks on a keystroke.");

  yame_usage_sec("Purpose:");
  yame_usage_text("Download reference assets into the shared store that every tool in the");
  yame_usage_text("suite reads, verifying each file against a digest this build pins.");
  yame_usage_text("-c puts them in the current directory instead, for a one-off or a");
  yame_usage_text("demo: no SHA256SUMS is written beside them, but each file is checked");
  yame_usage_text("against the same digest.");

  yame_usage_sec("Store:");
  /* a tool with its own store variable lists it first: that is what its
   * users export, and it is read ahead of YAME_DATA_HOME. The suite's other
   * variables come last, so a store moved with any one tool is found by all. */
  {
    char l[256];
    size_t k = 0, i;
    k += (size_t)snprintf(l + k, sizeof l - k, "Resolved in order: -d,");
    if (cfg_->tool_env)
      k += (size_t)snprintf(l + k, sizeof l - k, " $%s,", cfg_->tool_env);
    k += (size_t)snprintf(l + k, sizeof l - k, " $YAME_DATA_HOME,");
    for (i = 0; yame_assets_suite_env[i] && k < sizeof l; ++i) {
      if (cfg_->tool_env && strcmp(cfg_->tool_env, yame_assets_suite_env[i]) == 0)
        continue;
      k += (size_t)snprintf(l + k, sizeof l - k, " $%s,", yame_assets_suite_env[i]);
    }
    yame_usage_text(l);
  }
  yame_usage_text("${XDG_DATA_HOME:-~/.local/share}/yame");
  {
    const char *var = yame_assets_root_env(cfg_->tool_env);
    fprintf(stderr, "  %s%s: %s%s\n", yame_ui_green(),
            var ? var : (cfg_->tool_env ? cfg_->tool_env : "YAME_DATA_HOME"),
            root, yame_ui_reset());
  }
  yame_usage_sec("Mirror:");
  yame_usage_text("$YAME_ASSETS_MIRROR=<scheme://host[:port]> downloads from a site that");
  yame_usage_text("mirrors the public repositories, keeping each URL's path: the file at");
  yame_usage_text("https://raw.githubusercontent.com/zhou-lab/X/v1/f is fetched from <mirror>/zhou-lab/X/v1/f.");
  yame_usage_text("Every byte is still checked against the compiled-in digest.");

  yame_usage_sec("Options:");
  yame_usage_opt("-d <dir>", "Store root, overriding the environment.");
  yame_usage_opt("-f", "Re-download what is present, and replace a file the store's");
  yame_usage_cont("manifest records at a different digest than this build pins.");
  yame_usage_opt("-R", "A directory name also takes every directory beneath it:");
  yame_usage_cont("`-R hg38` is hg38, hg38/KYCG, hg38/data and hg38/models.");
  yame_usage_cont("Without it a name is one directory's own files. Works");
  yame_usage_cont("with -l and -n, so `-n -R hg38` says what that reaches.");
  yame_usage_opt("-l", "Dump the registry as TSV and exit: one row per file, with");
  yame_usage_cont("its size, digest, description and whether the store has it.");
  yame_usage_cont("Takes the same <name> and -g a fetch does, so `-l -g");
  yame_usage_cont("methscope hg38` is the dry run for fetching exactly that.");
  yame_usage_opt("-g <a,b>", "Only files matching every term: name, source, collection,");
  yame_usage_cont("title or upstream database. `-g chromatin` inside a");
  yame_usage_cont("knowledgebase, `-g celltype` across a whole genome.");
  yame_usage_opt("-n", "Say what would be fetched and stop, successfully. The same");
  yame_usage_cont("plan a fetch prints before asking, but it exits 0, so a");
  yame_usage_cont("script can check first. -l gives the same set as TSV.");
  yame_usage_opt("-y", "Fetch a whole folder without asking. A folder is confirmed");
  yame_usage_cont("first, since a short name can reach a lot -- `hg38` is 3.5");
  yame_usage_cont("GB. On a terminal the browser is the confirmation: it opens");
  yame_usage_cont("on the folder with exactly those files checked, f fetches,");
  yame_usage_cont("q leaves. Off one, a folder needs -y. A name that picks out");
  yame_usage_cont("ONE file needs neither: naming the file is the confirmation,");
  yame_usage_cont("so a documented fetch line runs in a script as it stands.");
  yame_usage_opt("-q", "No progress output.");
  yame_usage_opt("-u <url>", "Single-file form: what to download.");
  yame_usage_opt("-s <sha256>", "Single-file form: the digest it must have.");
  yame_usage_opt("-o <dest>", "Single-file form: where it goes (a path, not a dir).");
  yame_usage_opt("-c", "Into the current directory rather than the store.");
  yame_usage_opt("-h", "This help.");

  yame_usage_sec("Notes:");
  yame_usage_text("* Each file is verified against the digest this build carries for it;");
  yame_usage_text("  the SHA256SUMS a directory keeps records what was verified there.");
  yame_usage_text("  A file recorded at another digest is stale and needs -f to replace.");
  yame_usage_text("* `shasum -a 256 -c SHA256SUMS` in any store directory re-verifies it");
  yame_usage_text("  by hand, with none of this code involved.");
  /**
   * Say which of the two this build is, POSITIVELY.
   *
   * Only the bad case used to print, so a package test could assert the good
   * one only by grepping for the absence of a warning -- which passes on any
   * output at all, including a yame that failed to run. A line that is there
   * when downloading works is a thing a test can require.
   */
  if (yame_assets_have_curl())
    yame_usage_text("* Built with libcurl: fetch available.");
  else
    yame_usage_text("* THIS BUILD HAS NO LIBCURL: fetching is unavailable.");
  fprintf(stderr, "\n");

  return 1;
}

static const unit_t *find_unit(const char *dir) {
  for (size_t i = 0; i < YAME_ASSETS_N; ++i)
    if (strcmp(YAME_ASSETS[i].dir, dir) == 0) return &YAME_ASSETS[i];
  return NULL;
}


/* ------------------------------------------------------------ the browser
 *
 * The tree widget from src/ui.c, over the registry.
 *
 * Arranged the way someone thinks about the data -- species, then platform or
 * build, then its knowledgebase -- rather than the way it is published. The
 * registry is keyed on the upstream repo (InfiniumAnnotation/EPIC,
 * InfiniumAnnotation/EPIC/KYCG, KYCGKB/hg38, genomes/hg38), which is right for
 * fetching and wrong for choosing: two repos supply one genome, nobody thinks
 * of a knowledgebase as living somewhere other than the platform it is indexed
 * against, and nobody working on mouse wants to read past six human arrays.
 *
 * The store follows this tree rather than the registry's keying, so what you
 * see under a name is what lands under that name on disk: hg38/data in the
 * browser is <root>/hg38/data. Which repo filled a directory is in `source`.
 *
 *   human                       group
 *   +- hg38                     unit: a genome build or an array platform
 *   |  +- cpg_nocontig.cr       the unit index -- comes with anything else
 *   |  +- KYCG                  its knowledgebase
 *   |  |  +- CGI.20220904.cm    a set, with its .idx folded in
 *   |  +- seqinfo.tsv.gz
 *   +- MSA
 */

typedef struct {
  size_t asset;                 /* index into YAME_ASSETS */
  char   name[256];             /* the file within it */
} pick_t;

/**
 * Room for a selection.
 *
 * `a` takes the whole catalogue, which is 311 files today and grows with it,
 * so these are sized well past that rather than at it. Whatever the size,
 * hitting it is COUNTED and reported: a picker that quietly fetched the first
 * N of what you asked for would be worse than one that refused, because the
 * result looks like success.
 */
enum { PICK_MAX = 4096, CHOSEN_MAX = 256 };

typedef struct {
  char   root[4096];            /* the store */
  pick_t pick[PICK_MAX];
  size_t n_pick;
  size_t n_dropped;             /* wanted, but past PICK_MAX */
  int    force;

  /* What the user picked to USE, as opposed to everything that has to be
   * fetched alongside it -- an .idx and a unit index are part of the download
   * but are not what anyone asked to summarize against. */
  size_t chosen_asset[CHOSEN_MAX];
  char   chosen_name[CHOSEN_MAX][256];
  size_t n_chosen;
  size_t n_chosen_dropped;

  /* The widget holds these arrays and re-reads them after a commit, so a
   * fetch has to rewrite the entries in place -- otherwise a unit whose files
   * all just arrived keeps saying it has none of them. The arrays must not
   * move; the strings in them may. */
  char          **roots;
  unsigned char  *styles;
  size_t          n_roots;

  /* The widget's title buffer, which names the store: `d` rewrites it in
   * place when the store changes, since the widget holds the pointer. */
  char           *title;
  size_t          title_sz;

  /* A caller's narrowing (yame_browse_pick_opt): only these units are
   * offered, and these names arrive checked. NULL for the whole catalogue
   * and nothing checked, which is what the fetch browser always has. */
  const char *const *units;
  size_t              n_units;
  const char         *preselect;
  const char         *offer;      /* file globs a picker offers; NULL: all */

  /* `yame fetch <dir>` on a terminal: exactly the files the plan would move,
   * as "<unit index>|<file name>" -- the same key a file row carries. */
  char   **pre_keys;
  size_t   n_pre_keys;
} browse_t;

/* Defined below, next to the code that builds these rows in the first place:
 * a fetch has to put the counts back. */
static void refresh_roots(browse_t *b);

static const yame_asset_file_t *file_of(const unit_t *a,
                                        const char *name) {
  for (size_t i = 0; i < a->n_files; ++i)
    if (strcmp(yame_file_name(a->files[i]), name) == 0) return a->files[i];
  return NULL;
}

/* One store directory's own SHA256SUMS, cached a directory at a time. Every
 * caller walks a unit's files in a row, so a single slot turns what would be
 * a parse per file into a parse per directory. */
static struct {
  char             dir[4096];
  yame_sums_ent_t *ents;
  size_t           n;
  int              loaded;
} g_dsums;

/* A fetch rewrites the manifest under us, so the next question has to re-read
 * it -- otherwise a file that just arrived keeps answering from the tag it
 * replaced. */
static void dir_sums_forget(void) {
  free(g_dsums.ents);
  g_dsums.ents = NULL;
  g_dsums.n = 0;
  g_dsums.loaded = 0;
  g_dsums.dir[0] = '\0';
}

static const yame_sums_ent_t *dir_sums(const char *dir, size_t *n) {
  if (!g_dsums.loaded || strcmp(g_dsums.dir, dir) != 0) {
    char p[4096];
    /* The slot outlives every caller, so hand it back at exit rather than
     * leave the one reachable block in an otherwise clean leak report --
     * a known-benign block is exactly what makes the next real one easy to
     * miss. */
    static int at_exit_set = 0;
    if (!at_exit_set) { atexit(dir_sums_forget); at_exit_set = 1; }
    free(g_dsums.ents);
    g_dsums.ents = NULL;
    g_dsums.n = 0;
    if (yame_assets_join(p, sizeof(p), dir, YAME_ASSETS_SUMS_FILE) == 0)
      g_dsums.ents = yame_assets_sums_load_file(p, &g_dsums.n);
    if (!g_dsums.ents) g_dsums.n = 0;
    snprintf(g_dsums.dir, sizeof(g_dsums.dir), "%s", dir);
    g_dsums.loaded = 1;
  }
  *n = g_dsums.n;
  return g_dsums.ents;
}

/**
 * Does the directory's own manifest record this file at a DIFFERENT digest
 * than this build pins?
 *
 * Then the copy on disk is the previous tag's, not this one's -- the tag moved
 * under it. Calling that "present" is what let a superseded hg38_wg.updecx sit
 * in an upgraded store forever: the selection dropped it for having the right
 * name, so no digest was ever read and no re-download was ever proposed.
 *
 * A name the manifest does not list carries no such evidence, so it stays
 * present on the strength of existing -- an .idx that upstream does not
 * publish must not read as stale.
 */
static int file_superseded(const char *dir, const char *name,
                           const char *want) {
  if (!want || !*want) return 0;
  size_t n = 0;
  const yame_sums_ent_t *e = dir_sums(dir, &n);
  for (size_t i = 0; i < n; ++i)
    if (strcmp(e[i].name, name) == 0)
      return !yame_assets_digest_equal(e[i].sha, want);
  return 0;
}

static int file_present(const char *store_root, const unit_t *a,
                        const char *name) {
  const yame_asset_file_t *f = file_of(a, name);
  char dir[4096], path[4096];
  if (yame_assets_join(dir, sizeof(dir), store_root,
                       a->dir) != 0) return 0;
  if (yame_assets_join(path, sizeof(path), dir, name) != 0) return 0;
  if (!yame_assets_is_file(path)) return 0;
  return !file_superseded(dir, name, f ? f->sha256 : NULL);
}

/* One ladder for the whole suite: yame_human_size() in assets.c. */
static void human_size(uint64_t b, char *out, size_t n) {
  yame_human_size(b, out, n);
}

/* ---- how the catalogue is arranged ----
 *
 * Presentation only: the registry stays the truth. A unit this table does not
 * name still appears, under "other" -- a platform added upstream must show up
 * unannounced rather than vanish because nobody edited a list here.
 */
static const struct { const char *group, *unit; } UNIT_ORDER[] = {
  { "human",        "hg38"     },   /* the build first, then arrays newest */
  { "human",        "MSA"      },   /* first, since that is how anyone */
  { "human",        "EPICv2"   },   /* chooses one */
  { "human",        "EPIC"     },
  { "human",        "HM450"    },
  { "human",        "HM27"     },
  { "mouse",        "mm39"     },
  { "mouse",        "mm10"     },
  { "mouse",        "MM285"    },
  { "multispecies", "Mammal40" },
  { NULL, NULL }
};

/* ---- registry rows, read as units ---- */

/* Which unit a registry row belongs to, and where it sits inside it. */
static void unit_of(const unit_t *a, char *unit, size_t nu,
                    char *sub, size_t ns) {
  const char *slash = strchr(a->dir, '/');
  if (slash) {
    size_t l = (size_t)(slash - a->dir);
    if (l >= nu) l = nu - 1;
    memcpy(unit, a->dir, l);
    unit[l] = '\0';
    snprintf(sub, ns, "%s", slash + 1);
  } else {
    snprintf(unit, nu, "%s", a->dir);
    sub[0] = '\0';
  }
}

/**
 * Every row at or under a browser path.
 *
 * A name is a scope, not a key: "hg38" means the hg38 unit and everything in
 * it, the same as ticking that folder in the browser. That is also what makes
 * the name workable at all -- three targets (hg38, mm10, mm39) are published
 * by two sources each, and as a scope those simply select both instead of
 * being ambiguous.
 */
/* -R: a directory address also takes every directory beneath it. */
static int recursive_ = 0;

static size_t collect_scope(const char *path, const unit_t **out,
                            size_t cap) {
  /* A directory address is that directory's own files. `hg38` is the genome
   * annotation, not its knowledgebase and models beneath; `hg38/KYCG` is
   * addressed explicitly. (Decided 2026-09-19: a directory glob and a fetch
   * address select the same set, and sesame's "fetch the ordering" advice
   * must not pull 45 MB of sets under it.) With -R the name reaches down:
   * `hg38` is then hg38, hg38/KYCG, hg38/data and hg38/models, in table
   * order, so "everything for hg38" is one name rather than four to know. */
  size_t n = 0, pl = strlen(path);
  for (size_t i = 0; i < YAME_ASSETS_N && n < cap; ++i) {
    const char *d = YAME_ASSETS[i].dir;
    if (strcasecmp(d, path) == 0 ||
        (recursive_ && strncasecmp(d, path, pl) == 0 && d[pl] == '/'))
      out[n++] = &YAME_ASSETS[i];
  }
  return n;
}


/* An array platform, or a genome build? Decides what the row calls
 * itself. */
static int unit_is_array(const char *unit) {
  const unit_t *a = find_unit(unit);
  if (!a) return 0;
  for (size_t j = 0; j < a->n_files; ++j)
    if (strstr(yame_file_name(a->files[j]), ".ordering.tsv.gz")) return 1;
  return 0;
}

static const char *group_of(const char *unit) {
  for (size_t i = 0; UNIT_ORDER[i].unit; ++i)
    if (strcmp(UNIT_ORDER[i].unit, unit) == 0) return UNIT_ORDER[i].group;
  return "other";
}

static int unit_known(const char *unit) {
  for (size_t i = 0; i < YAME_ASSETS_N; ++i) {
    char u[256], s[256];
    unit_of(&YAME_ASSETS[i], u, sizeof(u), s, sizeof(s));
    if (strcmp(u, unit) == 0) return 1;
  }
  return 0;
}

/* The units of one group, in the table's order. */
static size_t group_units(const char *group, char units[][256], size_t cap) {
  size_t n = 0;
  if (strcmp(group, "other") != 0) {
    for (size_t i = 0; UNIT_ORDER[i].unit && n < cap; ++i)
      if (strcmp(UNIT_ORDER[i].group, group) == 0 &&
          unit_known(UNIT_ORDER[i].unit))
        snprintf(units[n++], 256, "%s", UNIT_ORDER[i].unit);
    return n;
  }
  for (size_t i = 0; i < YAME_ASSETS_N && n < cap; ++i) {
    char u[256], s[256];
    unit_of(&YAME_ASSETS[i], u, sizeof(u), s, sizeof(s));
    if (strcmp(group_of(u), "other") != 0) continue;
    size_t k = 0;
    for (; k < n; ++k) if (strcmp(units[k], u) == 0) break;
    if (k == n) snprintf(units[n++], 256, "%s", u);
  }
  return n;
}

static size_t all_groups(char groups[][64], size_t cap) {
  size_t n = 0;
  for (size_t i = 0; UNIT_ORDER[i].group && n < cap; ++i) {
    size_t k = 0;
    for (; k < n; ++k) if (strcmp(groups[k], UNIT_ORDER[i].group) == 0) break;
    if (k < n) continue;
    char units[32][256];
    if (group_units(UNIT_ORDER[i].group, units, 32))
      snprintf(groups[n++], 64, "%s", UNIT_ORDER[i].group);
  }
  char units[32][256];
  if (n < cap && group_units("other", units, 32))
    snprintf(groups[n++], 64, "%s", "other");
  return n;
}

/**
 * The file a unit cannot be read without: the probe ordering for an array,
 * the CpG coordinate stream for a genome.
 *
 * Everything else in the unit is a bit vector indexed against it, so it comes
 * along with any fetch from that unit rather than being something to
 * remember. It is shown at the top of the unit even when it is published
 * inside the knowledgebase repo -- where it came from is not the useful fact
 * about it.
 */
static int is_unit_index(const char *name) {
  return strstr(name, ".ordering.tsv.gz") != NULL ||
         strcmp(name, "cpg_nocontig.cr") == 0;
}

static int has_file(const unit_t *a, const char *name) {
  for (size_t i = 0; i < a->n_files; ++i)
    if (strcmp(yame_file_name(a->files[i]), name) == 0) return 1;
  return 0;
}

/**
 * An index rides with the file it indexes rather than being offered beside
 * it: it is unreadable alone and its data file is unusable without it, so
 * listing the two apart only invited picking one.
 *
 * Two spellings, because two ecosystems: YAME writes <x>.idx, tabix writes
 * <x>.tbi, and a genes.bed.gz without its .tbi is exactly as much use as a
 * .cm without its .idx.
 */
static const char *const COMPANION_SFX[] = { ".idx", ".tbi", NULL };

/* Is this file an index whose data file sits in the same directory? */
static int is_companion(const unit_t *a, const char *name) {
  size_t l = strlen(name);
  for (size_t k = 0; COMPANION_SFX[k]; ++k) {
    size_t sl = strlen(COMPANION_SFX[k]);
    char base[256];
    if (l <= sl || strcmp(name + l - sl, COMPANION_SFX[k]) != 0) continue;
    if (l - sl >= sizeof(base)) continue;
    memcpy(base, name, l - sl);
    base[l - sl] = '\0';
    if (has_file(a, base)) return 1;
  }
  return 0;
}

/* The index belonging to `name`, if the directory publishes one. Returns the
 * suffix that matched, so a row can say which kind it carries. */
static const char *companion_of(const unit_t *a, const char *name,
                                char *out, size_t n) {
  for (size_t k = 0; COMPANION_SFX[k]; ++k) {
    snprintf(out, n, "%s%s", name, COMPANION_SFX[k]);
    if (has_file(a, out)) return COMPANION_SFX[k];
  }
  out[0] = '\0';
  return NULL;
}

static uint64_t companion_size(const unit_t *a, const char *name) {
  char idxname[256];
  if (!companion_of(a, name, idxname, sizeof(idxname))) return 0;
  for (size_t i = 0; i < a->n_files; ++i)
    if (strcmp(yame_file_name(a->files[i]), idxname) == 0) return a->files[i]->size;
  return 0;
}

/* One thing the catalogue offers: a file, plus its index when it has one. */
typedef struct {
  const unit_t  *a;
  const yame_asset_file_t *f;
  const char              *paired;   /* ".idx" / ".tbi", or NULL */
  int                      required;
} ent_t;

/**
 * The entries one path publishes.
 *
 * `sub` is "" for the unit itself or "KYCG" for the node inside it; with
 * `recursive` it is ignored and the whole unit is walked, which is how a
 * folded row counts what is underneath it. Rows and counts come from this one
 * walk, so they cannot disagree.
 */
static size_t unit_entries(const char *unit, const char *sub, int recursive,
                           ent_t *out, size_t cap) {
  size_t n = 0;
  for (size_t i = 0; i < YAME_ASSETS_N && n < cap; ++i) {
    const unit_t *a = &YAME_ASSETS[i];
    char u[256], s[256];
    unit_of(a, u, sizeof(u), s, sizeof(s));
    if (strcmp(u, unit) != 0) continue;

    for (size_t j = 0; j < a->n_files && n < cap; ++j) {
      const yame_asset_file_t *f = a->files[j];
      if (is_companion(a, yame_file_name(f))) continue;

      int required = is_unit_index(yame_file_name(f));
      /* The row sits at its unit's store path. `required` governs ordering
       * and auto-inclusion, not placement. */
      const char *fsub = a->dir;
      const char *slash = strchr(fsub, '/');
      const char *at = slash ? slash + 1 : "";
      (void)s;
      if (!recursive && strcmp(at, sub) != 0) continue;

      char idxname[256];
      out[n].a = a;
      out[n].f = f;
      out[n].paired = companion_of(a, yame_file_name(f), idxname, sizeof(idxname));
      (void)idxname;
      out[n].required = required;
      ++n;
    }
  }
  return n;
}

/* The subdirectories directly below a unit -- in practice KYCG, but derived
 * rather than assumed. */
static size_t unit_subs(const char *unit, char subs[][256], size_t cap) {
  size_t n = 0;
  for (size_t i = 0; i < YAME_ASSETS_N && n < cap; ++i) {
    char u[256], s[256];
    unit_of(&YAME_ASSETS[i], u, sizeof(u), s, sizeof(s));
    if (strcmp(u, unit) != 0 || !s[0]) continue;

    size_t k = 0;
    for (; k < n; ++k) if (strcmp(subs[k], s) == 0) break;
    if (k == n) snprintf(subs[n++], 256, "%s", s);
  }
  return n;
}

static int ent_present(const char *store_root, const ent_t *e) {
  if (!file_present(store_root, e->a, yame_file_name(e->f))) return 0;
  if (e->paired) {
    char idxname[256];
    snprintf(idxname, sizeof(idxname), "%s%s", yame_file_name(e->f), e->paired);
    if (!file_present(store_root, e->a, idxname)) return 0;
  }
  return 1;
}

static void unit_counts(const char *store_root, const char *unit,
                        const char *sub, int recursive,
                        size_t *total, size_t *have) {
  static ent_t ents[1024];
  size_t n = unit_entries(unit, sub, recursive, ents, 1024);
  *total = n;
  *have = 0;
  for (size_t i = 0; i < n; ++i) if (ent_present(store_root, &ents[i])) ++*have;
}

/**
 * Does any directory behind this unit hold a tag this build does not pin?
 *
 * Presence is a file test, and a file test cannot see this: the three small
 * files of an older genomes tag are byte-identical to the new one, so the unit
 * looks half-fetched rather than blocked. Selecting the missing file then
 * fails inside the guard in yame_assets_fetch_dir(), which is a bad way to
 * learn that the whole directory needs -f. Ask the pin, and say so on the row.
 *
 * A directory merely one of OUR tags behind is not that: it fetches on its own
 * now, and the files whose digest moved already read as missing, so the row's
 * ordinary gauge tells the truth and no marker is wanted.
 */
/* How many of a unit's own files are stale: on disk, recorded at a digest
 * this build does not pin. */
static size_t unit_stale_count(const char *store_root, const char *unit) {
  const unit_t *a = find_unit(unit);
  size_t n = 0;
  if (!a) return 0;
  for (size_t j = 0; j < a->n_files; ++j)
    if (yame_file_state(store_root, a->files[j]) == YAME_STORE_STALE) ++n;
  return n;
}

static void group_counts(const char *store_root, const char *group,
                         size_t *total, size_t *have) {
  char units[32][256];
  size_t nu = group_units(group, units, 32);
  *total = 0; *have = 0;
  for (size_t i = 0; i < nu; ++i) {
    size_t t, h;
    unit_counts(store_root, units[i], "", 1, &t, &h);
    *total += t;
    *have += h;
  }
}

static unsigned char count_style(size_t total, size_t have) {
  if (total && have == total) return YAME_ROW_HAVE;
  if (!have) return YAME_ROW_MISSING;
  return YAME_ROW_PARTIAL;    /* partly here: neither claim is true */
}

/**
 * A five-cell gauge of how much of something is already in the store.
 *
 * Reading "12/36" takes a moment; seeing a half-filled bar does not, and the
 * question the browser mostly answers is "what am I still missing". Any
 * progress at all lights the first cell, so "a few of these" never looks
 * identical to "none of these".
 */
static void have_bar(size_t total, size_t have, char *out, size_t n) {
  const int cells = 5;
  int uni = yame_ui_unicode();
  const char *full = uni ? "▰" : "#";
  const char *empty = uni ? "▱" : ".";

  int on = 0;
  if (total && have) {
    on = (int)((have * (size_t)cells) / total);
    if (!on) on = 1;
    if (have == total) on = cells;
  }

  size_t o = 0;
  out[0] = '\0';
  for (int i = 0; i < cells; ++i) {
    const char *g = (i < on) ? full : empty;
    size_t l = strlen(g);
    if (o + l + 1 >= n) break;
    memcpy(out + o, g, l);
    o += l;
    out[o] = '\0';
  }
}

/**
 * The gauge, at a fixed width.
 *
 * The widget right-aligns this, which aligns its right edge -- so a ratio
 * that is three characters on one row and six on the next walks the bar left
 * and right down the column. Padding the ratio to a constant width pins both
 * ends of the field.
 */
static void counts_note(size_t total, size_t have, char *out, size_t n) {
  char bar[32], ratio[32];
  have_bar(total, have, bar, sizeof(bar));
  snprintf(ratio, sizeof(ratio), "%zu/%zu", have, total);
  snprintf(out, n, "%s %7s", bar, ratio);
}

/**
 * What a path is pinned at, and where it comes from.
 *
 * Both can be ambiguous, and for the same reason: a unit may draw on more
 * than one repo -- hg38 takes its knowledgebase from KYCGKB and its
 * annotation from zhou-lab/genomes, which are separately tagged. Naming one
 * of them for a row that spans two would be worse than saying nothing, so a
 * disagreement returns "mixed" for the tag and NULL for the upstream, and
 * the caller leaves the line out.
 *
 * `sub` NULL means the whole unit; "" means the unit's own directory; "KYCG"
 * means that node, which is where an answer is usually unambiguous even when
 * the unit's is not.
 */
static const char *path_tag(const char *unit, const char *sub) {
  const char *tag = NULL;
  for (size_t i = 0; i < YAME_ASSETS_N; ++i) {
    char u[256], s[256];
    unit_of(&YAME_ASSETS[i], u, sizeof(u), s, sizeof(s));
    if (strcmp(u, unit) != 0) continue;
    if (sub && strcmp(s, sub) != 0) continue;
    if (!tag) tag = YAME_ASSETS[i].tag;
    else if (strcmp(tag, YAME_ASSETS[i].tag) != 0) return "mixed";
  }
  return tag ? tag : "-";
}

static const char *path_upstream(const char *unit, const char *sub) {
  const char *url = NULL;
  for (size_t i = 0; i < YAME_ASSETS_N; ++i) {
    char u[256], s[256];
    unit_of(&YAME_ASSETS[i], u, sizeof(u), s, sizeof(s));
    if (strcmp(u, unit) != 0) continue;
    if (sub && strcmp(s, sub) != 0) continue;
    if (!url) url = YAME_ASSETS[i].source;
    else if (strcmp(url, YAME_ASSETS[i].source) != 0) return NULL;
  }
  return url;
}

static const char *unit_tag(const char *unit) { return path_tag(unit, NULL); }

/* ---- tree paths ---- */

/**
 * <unit>[/<sub>] -- the chain of keys the widget hands back.
 *
 * A species is a heading rather than a component: it groups the units under
 * it without being a level of its own. Nesting ten platforms under three
 * labels cost a keystroke to reach any of them and hid the rest until you
 * pressed it.
 */
typedef struct {
  char unit[256];
  char sub[256];
} bpath_t;

static void bpath_parse(const char *path, bpath_t *p) {
  p->unit[0] = p->sub[0] = '\0';

  const char *b = strchr(path, '/');
  size_t l = b ? (size_t)(b - path) : strlen(path);
  if (l >= sizeof(p->unit)) l = sizeof(p->unit) - 1;
  memcpy(p->unit, path, l);
  p->unit[l] = '\0';
  if (b) snprintf(p->sub, sizeof(p->sub), "%s", b + 1);
}

/* A heading row carries a rule in front of the species name, which is also
 * how one is told from a unit when it comes back as a path. */
static const char *group_mark(void) {
  /* Two cells, not a full rule: the heading shares a column with the unit
   * names, so anything wider pads every row beneath it. */
  return yame_ui_unicode() ? "▪ " : "= ";
}


/* Defined below, beside the -g filter they serve; the listing wants the same
 * matching so that -l is the dry run for a fetch. */
static void file_facets(const unit_t *a, const char *name,
                        char *out, size_t cap);
static int fetch_names(const unit_t *a, const char *store_root,
                       const char *const *names, size_t n_names,
                       const yame_fetch_opt_t *opt, char **err);
static int facets_match(const char *facets, const char *terms);
static size_t collect_scope(const char *path, const unit_t **out,
                            size_t cap);

/* One registry entry with the files named out of it: empty list = the whole
 * directory. Declared here because `-l` prints exactly the selection a fetch
 * would make, so both paths take this. */
#define SEL_ONLY_MAX 64
typedef struct {
  const unit_t *a;
  const char *only[SEL_ONLY_MAX];  /* empty: the whole directory */
  size_t n_only;
} sel_t;
static int file_wanted(const unit_t *a, const char *name,
                       const char *filter,
                       const char *const *only, size_t n_only);

/**
 * Every file this build knows about, one per line, as TSV.
 *
 * The browser is the only other thing that can see individual files, and it
 * needs a terminal and a person. Anything else -- a script, a CI job, an
 * agent -- could previously enumerate the twenty fetchable specs and nothing
 * below them, so "is there a ChromHMM set for MSA?" had no non-interactive
 * answer. This is that answer: the registry, dumped.
 *
 * TSV rather than the aligned columns the summary uses, because the consumer
 * here is cut(1), not an eye. Descriptions come from the registry row: each
 * file carries its own.
 */
static int dump_registry(const char *dopt, const sel_t *sel, size_t n_sel,
                         const char *filter) {
  char root[4096];
  yame_assets_root(dopt, cfg_->tool_env, root, sizeof(root));

  printf("target\tsource\ttag\tstore_path\tfile\tbytes\tsha256\tdir_state\t"
         "local\tdescription\n");

  for (size_t i = 0; i < n_sel; ++i) {
    const unit_t *a = sel[i].a;
    char dir[4096];
    yame_assets_join(dir, sizeof(dir), root, a->dir);

    (void)dir;
    for (size_t j = 0; j < a->n_files; ++j) {
      const yame_asset_file_t *f = a->files[j];
      const char *state = "-";
      switch (yame_file_state(root, f)) {
      case YAME_STORE_CURRENT: state = "current"; break;
      case YAME_STORE_ABSENT:  state = "absent";  break;
      case YAME_STORE_STALE:   state = "stale";   break;
      default: break;
      }

      /* A name that picked out one file lists that file, not its directory:
       * `-l <file>` is the dry run for fetching that file, which is what the
       * usage text has always promised. */
      if (!file_wanted(a, yame_file_name(f), filter, sel[i].only, sel[i].n_only))
        continue;

      /* A tab or newline inside a title would shift every later column, so
       * fold whitespace to single spaces on the way out. */
      char title[512];
      const char *src = f->title && *f->title ? f->title : "-";
      size_t n = 0;
      for (; *src && n + 1 < sizeof(title); ++src) {
        char c = (*src == '\t' || *src == '\n' || *src == '\r') ? ' ' : *src;
        if (c == ' ' && n && title[n - 1] == ' ') continue;
        title[n++] = c;
      }
      title[n] = '\0';

      printf("%s\t%s\t%s\t%s\t%s\t%" PRIu64 "\t%s\t%s\t%s\t%s\n",
             a->dir, a->source, a->tag, a->dir, yame_file_name(f),
             f->size, f->sha256 ? f->sha256 : "-", state,
             file_present(root, a, yame_file_name(f)) ? "yes" : "no", title);
    }
  }
  return 0;
}

/**
 * Everything a file can be matched on, gathered into one string.
 *
 * The filename and the row it came from, plus what the catalogue already
 * knows about that kind of file: which collections it appears in, its short
 * title, and the upstream database it came from. Deliberately NOT the biology
 * or processing prose -- those run to hundreds of characters and mention
 * "CpG" or "methylation" in nearly every row, so matching them would make
 * every common word select the whole catalogue.
 */
static void file_facets(const unit_t *a, const char *name,
                        char *out, size_t cap) {
  const yame_asset_file_t *f = file_of(a, name);
  char u[256], sb[256];
  unit_of(a, u, sizeof(u), sb, sizeof(sb));
  snprintf(out, cap, "%s %s %s %s %s %s", name, a->source, u, sb,
           f && f->title  ? f->title  : "",
           f && f->source ? f->source : "");
}

/* All terms must appear: narrowing is the point, so a second term that
 * widened the result would be a surprise. */
static int facets_match(const char *facets, const char *terms) {
  char buf[512];
  snprintf(buf, sizeof(buf), "%s", terms);
  for (char *save = NULL, *t = strtok_r(buf, ",", &save);
       t; t = strtok_r(NULL, ",", &save)) {
    while (*t == ' ') ++t;
    if (!*t) continue;
    if (!strcasestr(facets, t)) return 0;
  }
  return 1;
}

/* ---- expanding ----
 *
 * Every row is "<content>\t<trailing field>", and the widget right-aligns
 * that trailing field at the margin. So a unit's gauge, a folder's gauge and
 * a file's size all land in one column however deep the row sits -- which
 * padding to a fixed width cannot do, since each level of nesting spends four
 * more cells on the indent.
 */

/* Lines its label up under a unit row's KIND column. */
#define SUB_ROW_FMT "%-11s %-7s"

/**
 * The trailing field, which the widget right-aligns: a tag, then whatever the
 * row measures -- a size for a file, a gauge for anything holding files.
 *
 * One shape for every kind of row, because they share the column. Putting a
 * file's tag beside its size while a folder's sat in a column forty cells to
 * the left meant the same fact appeared in two places depending on what you
 * were looking at, which is exactly when a column stops being readable.
 * Fixed widths: the field is right-aligned as a whole, so anything variable
 * inside it walks the tag back and forth down the list.
 */
#define TAIL_FMT "%-5s %13s"

/* Keys are "<asset index>|<filename>" for a file and the plain component name
 * for a directory, which is what makes a tree path parseable back. */
/**
 * One file's row: its name, the tag it comes from, and its size.
 *
 * Every row carries its tag, not only those under a unit reporting "mixed".
 * Showing it just where it disambiguates sounds tidier and is worse in use:
 * it makes the absence of a tag meaningful, so reading a row correctly means
 * first noticing what its folder said several lines above. A column that is
 * always there is read without that step.
 */
static void emit_entry(browse_t *b, yame_ui_kids_t *out, const ent_t *e) {
  size_t idx = (size_t)(e->a - YAME_ASSETS);
  int here = ent_present(b->root, e);

  char sz[24], name[256], line[352], key[288];
  /* A pair is one row, so it is one size too: the .idx is part of what a
   * fetch of this row will cost. */
  human_size(e->f->size + (e->paired ? companion_size(e->a, yame_file_name(e->f)) : 0),
             sz, sizeof(sz));
  snprintf(name, sizeof(name), "%s%s%s", yame_file_name(e->f), e->paired ? " +" : "",
           e->paired ? e->paired + 1 : "");
  /* Both fields fixed-width: the tail is right-aligned as a whole, so a size
   * that varies in width would walk the tag column left and right. */
  /* A stale file -- on disk, but from a release this build does not pin --
   * says so where the eye already is, and is offered again: green would
   * mean "nothing to do here", which is the one thing it is not. */
  int stale = yame_file_state(b->root, e->f) == YAME_STORE_STALE;
  char tail[32];
  snprintf(tail, sizeof(tail), "%s%s", stale ? "stale " : "", sz);
  snprintf(line, sizeof(line), "%s\t" TAIL_FMT, name, e->a->tag, tail);
  snprintf(key, sizeof(key), "%zu|%s", idx, yame_file_name(e->f));

  out->rows[out->n]   = wzstrdup(line);
  out->keys[out->n]   = wzstrdup(key);
  out->styles[out->n] = (unsigned char)(here && !stale ? YAME_ROW_HAVE
                                             : e->required ? YAME_ROW_REQUIRED
                                                           : YAME_ROW_MISSING);
  ++out->n;
}

/**
 * One level: the units of a group, or a unit's index, subdirectories and
 * files.
 *
 * Files are offered per file rather than per directory because a
 * knowledgebase holds dozens of sets and most callers want a few.
 */
/* Does a picker offer this file? With `offer` set, only files whose store
 * path matches one of its comma-separated globs (fnmatch, so `*` also crosses
 * `/`): `*.cm` for any mask, `* /KYCG/ *.cm` for knowledgebase sets alone. A
 * picker whose caller can only use sets must not list the row list or a
 * model, or the choice is refused after it is made. */
static int offered_file(const browse_t *b, const yame_asset_file_t *f) {
  if (!b->offer) return 1;
  char buf[512];
  snprintf(buf, sizeof buf, "%s", b->offer);
  for (char *save = NULL, *g = strtok_r(buf, ",", &save); g; g = strtok_r(NULL, ",", &save)) {
    while (*g == ' ') ++g;
    if (*g && fnmatch(g, f->store_path, 0) == 0) return 1;
  }
  return 0;
}

static void bx_expand(void *ctx, const char *path, yame_ui_kids_t *out) {
  browse_t *b = ctx;
  bpath_t p;
  bpath_parse(path, &p);

  enum { CAP = 512 };
  out->rows   = wzcalloc(CAP, sizeof(char *));
  out->keys   = wzcalloc(CAP, sizeof(char *));
  out->styles = wzcalloc(CAP, 1);
  out->branch = wzcalloc(CAP, 1);
  if (!out->rows || !out->keys || !out->styles || !out->branch) return;

  static ent_t ents[CAP];
  size_t n_ents = unit_entries(p.unit, p.sub, 0, ents, CAP);

  /* The index first: it is what everything below it is addressed against. */
  for (size_t i = 0; i < n_ents && out->n < CAP; ++i)
    if (ents[i].required && offered_file(b, ents[i].f)) emit_entry(b, out, &ents[i]);

  /* Then the subdirectories, only at the top of a unit. */
  if (!p.sub[0]) {
    char subs[8][256];
    size_t n_subs = unit_subs(p.unit, subs, 8);
    for (size_t i = 0; i < n_subs && out->n < CAP; ++i) {
      size_t total, have;
      unit_counts(b->root, p.unit, subs[i], 0, &total, &have);
      if (!total) continue;
      if (b->offer) {                    /* nothing in it to offer: no row */
        static ent_t sub_ents[CAP];
        size_t ns = unit_entries(p.unit, subs[i], 0, sub_ents, CAP), k = 0;
        while (k < ns && !offered_file(b, sub_ents[k].f)) ++k;
        if (k == ns) continue;
      }

      char note[64], line[256], subdir[512];
      counts_note(total, have, note, sizeof(note));
      snprintf(subdir, sizeof(subdir), "%s/%s", p.unit, subs[i]);
      size_t stale = unit_stale_count(b->root, subdir);
      if (stale) snprintf(note, sizeof(note), "%zu stale: -f", stale);
      snprintf(line, sizeof(line), SUB_ROW_FMT "\t" TAIL_FMT, subs[i], "sets",
               path_tag(p.unit, subs[i]), note);

      out->rows[out->n]   = wzstrdup(line);
      out->keys[out->n]   = wzstrdup(subs[i]);
      out->styles[out->n] = count_style(total, have);
      out->branch[out->n] = 1;
      ++out->n;
    }
  }

  for (size_t i = 0; i < n_ents && out->n < CAP; ++i)
    if (!ents[i].required && offered_file(b, ents[i].f)) emit_entry(b, out, &ents[i]);
}

/* ---- the info pane ----
 *
 * A detail callback rather than a modal panel, so the arrow keys keep working
 * while it is open and the text follows the cursor. It runs on every redraw,
 * which is affordable because the lookup is a scan of a compiled-in table and
 * the wrapping is a few hundred bytes of formatting.
 *
 * Ported from kycg, which is also where the descriptions started.
 */

#define INFO_MAX_LINES 64

typedef struct {
  char *line[INFO_MAX_LINES];
  int   n;
  int   cols;
} info_lay_t;

static void lay_push(info_lay_t *L, const char *s) {
  if (L->n >= INFO_MAX_LINES) return;
  L->line[L->n] = wzstrdup(s ? s : "");
  if (L->line[L->n]) ++L->n;
}

static void lay_free(info_lay_t *L) {
  for (int i = 0; i < L->n; ++i) free(L->line[i]);
  L->n = 0;
}

/* Name in bold, what it is in cyan beside it, then a blank line. */
static void lay_head(info_lay_t *L, const char *name, const char *title) {
  char buf[1024];
  snprintf(buf, sizeof(buf), "  %s%s%s   %s%s%s", yame_ui_bold(), name,
           yame_ui_reset(), yame_ui_cyan(), title ? title : "",
           yame_ui_reset());
  lay_push(L, buf);
  lay_push(L, "");
}

/**
 * Wrap `text` into the pane, under an optional dim label.
 *
 * NULL `label` means running prose at the left margin; otherwise the label is
 * printed once in a fixed-width gutter and continuation lines align under the
 * text rather than under the label, so a three-line citation reads as one
 * field.
 */
static void lay_wrap(info_lay_t *L, const char *label, const char *text) {
  if (!text || !*text) return;

  const int gutter = label ? 14 : 2;   /* "  processing  " is the widest */
  int avail = (L->cols > 0 ? L->cols : yame_ui_cols()) - gutter - 2;
  if (avail < 20) avail = 20;

  const char *p = text;
  int first = 1;
  while (*p) {
    while (*p == ' ') ++p;
    if (!*p) break;

    /* Longest prefix that fits, broken at the last space; a single word
     * longer than the line is emitted whole and allowed to be truncated,
     * which beats hyphenating a DOI. */
    size_t rest = strlen(p), take = rest;
    if (rest > (size_t)avail) {
      size_t brk = 0;
      for (size_t i = 0; i < (size_t)avail; ++i) if (p[i] == ' ') brk = i;
      take = brk ? brk : (size_t)avail;
    }

    char buf[1024], head[64];
    if (label && first)
      snprintf(head, sizeof(head), "  %s%-*s%s", yame_ui_dim(), gutter - 4,
               label, yame_ui_reset());
    else
      snprintf(head, sizeof(head), "%*s", gutter, "");

    snprintf(buf, sizeof(buf), "%s%s%.*s", head, label && first ? "  " : "",
             (int)take, p);
    lay_push(L, buf);

    p += take;
    first = 0;
  }
  lay_push(L, "");
}

/* Where a file comes from and where it lands. Facts from the registry, so
 * every file says something even when the table describes none of them. */
static void lay_provenance(info_lay_t *L, const unit_t *a,
                           const yame_asset_file_t *f, const char *paired) {
  char buf[1024];

  snprintf(buf, sizeof(buf), "%s @ %s", a->source, a->tag);
  lay_wrap(L, "upstream", buf);

  snprintf(buf, sizeof(buf), "%s/%s", a->dir, f ? yame_file_name(f) : "");
  lay_wrap(L, "store", buf);

  if (f) {
    snprintf(buf, sizeof(buf), "%.16s...%s%s", f->sha256,
             paired ? "   ships with its " : "", paired ? paired : "");
    lay_wrap(L, "sha256", buf);
  }
}

static void bx_detail(void *ctx, const char *path, const char *key, int cols,
                      yame_ui_detail_t *out) {
  browse_t *b = ctx;
  out->rows = NULL;
  out->n = 0;

  info_lay_t L = {0};
  L.cols = cols > 20 ? cols : 20;

  bpath_t p;
  bpath_parse(path, &p);

  const char *bar = key ? strchr(key, '|') : NULL;

  /* No branch for a species heading: the cursor never lands on one (a
   * heading never opens and is stepped over), so the pane is never asked
   * about it. */
  if (!bar) {
    /* A unit, or the knowledgebase inside one. */
    size_t total, have;
    unit_counts(b->root, p.unit, p.sub, p.sub[0] ? 0 : 1, &total, &have);

    char buf[1024];
    snprintf(buf, sizeof(buf), "%s%s%s", p.unit, p.sub[0] ? "/" : "", p.sub);
    lay_head(&L, buf, p.sub[0] ? "knowledgebase"
                               : unit_is_array(p.unit) ? "array platform"
                                                       : "genome build");

    if (p.sub[0])
      lay_wrap(&L, NULL, "Annotations projected onto this unit's row space, "
                         "one .cm per set, each with its .idx. What `kycg "
                         "test` enriches against.");
    else
      lay_wrap(&L, NULL, "Everything here is indexed against one row space. "
                         "The index at the top of the list is fetched with "
                         "anything else taken from this unit.");

    snprintf(buf, sizeof(buf), "%zu of %zu files already in the store",
             have, total);
    lay_wrap(&L, "in store", buf);

    /* Only when it is one answer: a unit spanning two repos gets neither
     * line rather than one repo's name standing for both. */
    const char *sub = p.sub[0] ? p.sub : NULL;
    const char *up = path_upstream(p.unit, sub);
    const char *tg = path_tag(p.unit, sub);
    if (up) {
      snprintf(buf, sizeof(buf), "%s @ %s", up, tg);
      lay_wrap(&L, "upstream", buf);
    } else {
      lay_wrap(&L, "upstream", "several repos -- open a row to see which");
    }
    lay_wrap(&L, "tag", tg);
  } else {
    size_t idx = (size_t)strtoul(key, NULL, 10);
    const char *name = bar + 1;
    if (idx >= YAME_ASSETS_N) { lay_free(&L); return; }
    const unit_t *a = &YAME_ASSETS[idx];

    const yame_asset_file_t *f = NULL;
    for (size_t i = 0; i < a->n_files; ++i)
      if (strcmp(yame_file_name(a->files[i]), name) == 0) { f = a->files[i]; break; }

    char idxname[256];
    const char *paired = companion_of(a, name, idxname, sizeof(idxname));

    if (f && f->title && *f->title) {
      lay_head(&L, name, f->title);
      lay_wrap(&L, NULL, f->description);
      lay_wrap(&L, "source", f->source);
      lay_wrap(&L, "citation", f->citation);
    } else {
      char buf[512];
      snprintf(buf, sizeof(buf), "  %s%s%s   %s(nothing recorded about this "
               "file)%s", yame_ui_bold(), name, yame_ui_reset(),
               yame_ui_dim(), yame_ui_reset());
      lay_push(&L, buf);
      lay_push(&L, "");
    }
    lay_provenance(&L, a, f, paired);
    if (f) {
      char adv[1024];
      if (yame_store_state(cfg_, f->store_path, adv, sizeof adv) == YAME_STORE_STALE)
        lay_wrap(&L, "in store", adv);
    }
  }

  out->rows = wzmalloc((size_t)(L.n ? L.n : 1) * sizeof(char *));
  if (!out->rows) { lay_free(&L); return; }
  for (int i = 0; i < L.n; ++i) out->rows[i] = L.line[i];
  out->n = (size_t)L.n;
  /* Ownership moves to the widget; do not free the strings here. */
  L.n = 0;
  lay_free(&L);
}

/**
 * What a row can be found by, beyond the filename it displays.
 *
 * Same fields the -g filter uses, for the same reason: the tree shows a
 * filename, so without this `/methscope` or `/chromatin` would find nothing,
 * though both are true of rows on the screen. A folder answers for itself --
 * its unit and source -- so filtering for a source keeps its folders as well
 * as its files.
 */
static const char *bx_facets(void *ctx, const char *path, const char *key) {
  static char buf[1024];
  (void)ctx;
  bpath_t bp;
  bpath_parse(path, &bp);
  if (!bp.unit[0]) return NULL;

  const char *bar = key ? strchr(key, '|') : NULL;
  if (bar) {                                  /* a file: ask the registry */
    size_t idx = (size_t)atoi(key);
    if (idx < YAME_ASSETS_N) {
      file_facets(&YAME_ASSETS[idx], bar + 1, buf, sizeof(buf));
      return buf;
    }
  }
  /* a unit or folder: its own name, and every source that fills it */
  size_t n = (size_t)snprintf(buf, sizeof(buf), "%s %s", bp.unit, bp.sub);
  for (size_t i = 0; i < YAME_ASSETS_N && n + 32 < sizeof(buf); ++i) {
    char u[256], sb[256];
    unit_of(&YAME_ASSETS[i], u, sizeof(u), sb, sizeof(sb));
    if (strcmp(u, bp.unit) != 0) continue;
    if (bp.sub[0] && strcmp(sb, bp.sub) != 0) continue;
    n += (size_t)snprintf(buf + n, sizeof(buf) - n, " %s", YAME_ASSETS[i].source);
  }
  return buf;
}

/* ---- choosing and fetching ---- */

static void pick_add(browse_t *b, size_t asset, const char *name) {
  if (b->n_pick >= PICK_MAX) { ++b->n_dropped; return; }
  for (size_t i = 0; i < b->n_pick; ++i)
    if (b->pick[i].asset == asset && strcmp(b->pick[i].name, name) == 0) return;
  b->pick[b->n_pick].asset = asset;
  snprintf(b->pick[b->n_pick].name, sizeof(b->pick[0].name), "%s", name);
  ++b->n_pick;
}

static void bx_accept(void *ctx, const char *path, const char *key) {
  browse_t *b = ctx;
  (void)path;
  const char *bar = key ? strchr(key, '|') : NULL;
  if (!bar) return;

  size_t idx = (size_t)strtoul(key, NULL, 10);
  if (idx >= YAME_ASSETS_N) return;
  const unit_t *a = &YAME_ASSETS[idx];
  const char *name = bar + 1;
  pick_add(b, idx, name);

  /* The pieces that are not choices. A .cm without its .idx is unusable, and
   * everything in a unit is addressed against that unit's index -- so both
   * come along rather than being something to remember. */
  char idxname[256];
  if (companion_of(a, name, idxname, sizeof(idxname)))
    pick_add(b, idx, idxname);

  char unit[256], sub[256];
  unit_of(a, unit, sizeof(unit), sub, sizeof(sub));
  static ent_t ents[512];
  size_t n = unit_entries(unit, "", 0, ents, 512);
  for (size_t i = 0; i < n; ++i)
    if (ents[i].required)
      pick_add(b, (size_t)(ents[i].a - YAME_ASSETS), yame_file_name(ents[i].f));
}

/**
 * One panel line, laid out the same way every time: a glyph in its own
 * column, a label, and a value at a fixed column.
 *
 * The panel used to print whatever each site had to say -- "ok" alone on a
 * line, a bare count, a path run together with a number -- which reads as
 * debug output rather than as the tool reporting. One shape makes it a small
 * table, and a small table is read at a glance.
 */
#define PANEL_LABELW 42

static void panel_row(int line, const char *color, const char *glyph,
                      const char *label, const char *value) {
  char lbl[700];
  snprintf(lbl, sizeof(lbl), "%s", label ? label : "");

  int pad = PANEL_LABELW - (int)strlen(lbl);
  if (pad < 1) pad = 1;
  yame_ui_panel_line(line, "  %s%s%s  %s%*s%s%s%s",
                     color ? color : "", glyph ? glyph : " ", yame_ui_reset(),
                     lbl, pad, "",
                     yame_ui_dim(), value ? value : "", yame_ui_reset());
}

/* Bytes moving, painted onto the row of the file they belong to. The panel
 * at the foot says what is happening overall; this says which of the rows you
 * are looking at is the one being fetched. */
typedef struct {
  size_t asset;
  char   key[288];
} fprog_t;

/* An .idx has no row of its own -- it rides with its .cm -- so its progress
 * is painted onto the row that asked for it. */
static void fprog_key(fprog_t *p, const char *name) {
  snprintf(p->key, sizeof(p->key), "%zu|%s", p->asset, name);
  if (yame_ui_tree_progress(p->key, 0, 1)) return;

  size_t l = strlen(name);
  if (l > 4 && strcmp(name + l - 4, ".idx") == 0) {
    char base[256];
    if (l - 4 < sizeof(base)) {
      memcpy(base, name, l - 4);
      base[l - 4] = '\0';
      snprintf(p->key, sizeof(p->key), "%zu|%s", p->asset, base);
    }
  }
}

static void fp_begin(void *ud, const char *name, uint64_t total) {
  (void)total;
  fprog_key((fprog_t *)ud, name);
}

static void fp_progress(void *ud, uint64_t now, uint64_t total) {
  fprog_t *p = ud;
  if (p->key[0]) yame_ui_tree_progress(p->key, now, total ? total : now + 1);
}

static void fp_done(void *ud, const char *name, uint64_t bytes, int ok) {
  fprog_t *p = ud;
  (void)bytes;
  if (!p->key[0]) return;

  /**
   * Green as soon as THIS file lands, not when the batch does.
   *
   * A row standing for a pair waits for the .idx: the .cm arriving is not the
   * row arriving, and a green row whose index is still in flight is a claim
   * the store cannot honour. The .idx is the second of the two, so settling
   * on it settles the pair.
   */
  int done = ok;
  if (ok && p->asset < YAME_ASSETS_N) {
    char idxname[256];
    if (companion_of(&YAME_ASSETS[p->asset], name, idxname, sizeof(idxname)))
      done = 0;                    /* its .idx is still to come */
  }
  yame_ui_tree_settle(p->key, done);
  p->key[0] = '\0';
}

/* What `f` is about to do, before it does it. A selection is easy to build up
 * without noticing -- space on a directory takes everything under it -- and a
 * knowledgebase runs to tens of megabytes, so the size is worth seeing while
 * it can still be reconsidered. */
static int confirm_plan(browse_t *b) {
  size_t n_files = 0, n_dirs = 0, unknown = 0;
  uint64_t bytes = 0;
  size_t seen[512], n_seen = 0;
  /* Whether any pick for that asset lands in the asset's OWN directory: only
   * those write its manifest, so only those can carry its tag forward. */

  for (size_t i = 0; i < b->n_pick; ++i) {
    size_t idx = b->pick[i].asset;
    if (idx >= YAME_ASSETS_N) continue;
    const unit_t *a = &YAME_ASSETS[idx];
    /* Count what the rows show. A .cm and its .idx are one row and one
     * "x of y" everywhere else in the browser, so the header counts pairs
     * too; the .idx is still a real download, so its bytes still count. */
    if (!is_companion(a, b->pick[i].name)) ++n_files;

    uint64_t sz = 0;
    for (size_t j = 0; j < a->n_files; ++j)
      if (strcmp(yame_file_name(a->files[j]), b->pick[i].name) == 0) {
        sz = a->files[j]->size;
        break;
      }
    if (sz) bytes += sz; else ++unknown;

    size_t k = 0;
    for (; k < n_seen; ++k) if (seen[k] == idx) break;
    if (k == n_seen && n_seen < 512) seen[n_seen++] = idx;
  }
  n_dirs = n_seen;

  char sz[32], label[256], value[128];
  human_size(bytes, sz, sizeof(sz));

  snprintf(label, sizeof(label), "Fetch %zu file%s from %zu director%s",
           n_files, n_files == 1 ? "" : "s", n_dirs, n_dirs == 1 ? "y" : "ies");
  if (!bytes)          snprintf(value, sizeof(value), "size not published");
  else if (unknown)    snprintf(value, sizeof(value), "at least %s", sz);
  else                 snprintf(value, sizeof(value), "%s", sz);

  int row = 1;
  yame_ui_panel_open(4);
  panel_row(0, yame_ui_bold(), yame_ui_unicode() ? "⤓" : ">", label, value);
  if (b->n_dropped) {
    char note[128];
    snprintf(note, sizeof(note), "%zu more did not fit and will NOT be fetched",
             b->n_dropped);
    panel_row(row, yame_ui_red(), yame_ui_cross(), "selection truncated", note);
  } else if (unknown) {
    panel_row(row, NULL, " ", "", "some sizes are not published upstream");
  } else {
    yame_ui_panel_line(row, "");
  }
  return yame_ui_panel_confirm(row + 1, "   Proceed?", 1);
}

/* Fetch everything picked, grouped by directory so each group costs a single
 * manifest request. `in_widget` reports through the panel; otherwise this is
 * an ordinary command-line fetch. */
static void fetch_picks(browse_t *b, int in_widget, int *ok_out,
                        int *bad_out, uint64_t *bytes_out);

/* The "use these" action: same collection as a fetch, plus a record of what
 * was actually asked for. A .cm's .idx and a unit's index come along with the
 * download but are not what anyone wants to summarize against. */
static void bx_choose(void *ctx, const char *path, const char *key) {
  browse_t *b = ctx;
  bx_accept(ctx, path, key);              /* everything the fetch needs */

  const char *bar = key ? strchr(key, '|') : NULL;
  if (!bar) return;
  if (b->n_chosen >= CHOSEN_MAX) { ++b->n_chosen_dropped; return; }

  b->chosen_asset[b->n_chosen] = (size_t)strtoul(key, NULL, 10);
  snprintf(b->chosen_name[b->n_chosen], sizeof(b->chosen_name[0]), "%s",
           bar + 1);
  ++b->n_chosen;
}

/* Runs inside the widget: the tree stays on screen, this draws in a panel and
 * onto the rows being fetched, and the tree resumes with its rows reloaded so
 * the new state shows. */
static int bx_commit(void *ctx) {
  browse_t *b = ctx;
  int ok = 0, bad = 0;

  if (!b->n_pick) return 0;

  if (!confirm_plan(b)) {
    yame_ui_panel_close();
    b->n_pick = 0;               /* the checks stay; the fetch does not run */
    b->n_dropped = 0;
    return 0;
  }

  uint64_t moved = 0;
  yame_ui_panel_open(4);
  fetch_picks(b, 1, &ok, &bad, &moved);

  char label[128], value[64];
  snprintf(label, sizeof(label), "%d file%s fetched", ok, ok == 1 ? "" : "s");
  human_size(moved, value, sizeof(value));
  panel_row(0, yame_ui_green(), yame_ui_check(), label,
            moved ? value : "");

  if (bad) {
    snprintf(label, sizeof(label), "%d failed", bad);
    panel_row(1, yame_ui_red(), yame_ui_cross(), label,
              "nothing was written for these");
  } else {
    yame_ui_panel_line(1, "");
  }
  yame_ui_panel_pause(2, "   press any key to return to the catalogue");
  yame_ui_panel_close();
  refresh_roots(b);
  b->n_pick = 0;
  b->n_dropped = 0;
  return 1;
}

/**
 * `d` in either browser: point it at a different store.
 *
 * The store is the one thing the browser shows that cannot otherwise be
 * changed from inside it -- every other question is about its contents -- so
 * without this a wrong $YAME_DATA_HOME meant quitting and starting again with
 * -d. The path is edited in place (a leading ~ is the home directory); a path
 * that does not exist yet is a new, empty store, as -d would make it; one
 * that exists and is not a directory is refused. Returning 1 makes the tree
 * reload everything, folds and checks included: they described the old store.
 */
static int bx_on_key(void *ctx, char key, const char *path, const char *node_key) {
  browse_t *b = ctx;
  (void) path; (void) node_key;           /* about the store, not a row */
  if (key != 'd') return 0;

  char buf[4096];
  snprintf(buf, sizeof buf, "%s", b->root);
  yame_ui_panel_open(3);
  yame_ui_panel_line(0, "  %s%s%s", yame_ui_dim(),
                     "store directory (enter to accept, esc to cancel)", yame_ui_reset());
  int ok = yame_ui_panel_ask(1, "store:", buf, sizeof buf);
  if (!ok || !buf[0] || strcmp(buf, b->root) == 0) { yame_ui_panel_close(); return 0; }

  char want[4096];
  const char *home = getenv("HOME");
  if (buf[0] == '~' && (buf[1] == '/' || !buf[1]) && home)
    snprintf(want, sizeof want, "%s%s", home, buf + 1);
  else
    snprintf(want, sizeof want, "%s", buf);
  if (yame_assets_is_file(want)) {
    char msg[4200];
    snprintf(msg, sizeof msg, "  %s is not a directory; the store is unchanged", want);
    yame_ui_panel_pause(2, msg);
    yame_ui_panel_close();
    return 0;
  }
  yame_ui_panel_close();

  snprintf(b->root, sizeof b->root, "%s", want);
  if (b->title) snprintf(b->title, b->title_sz, "store (d): %s", b->root);
  b->n_pick = 0; b->n_dropped = 0;
  refresh_roots(b);
  return 1;
}

static void fetch_picks(browse_t *b, int in_widget, int *ok_out, int *bad_out,
                        uint64_t *bytes_out) {
  int ok = 0, bad = 0;
  uint64_t moved = 0;

  for (size_t i = 0; i < b->n_pick; ++i) {
    size_t idx = b->pick[i].asset;
    /* Already fetched as part of an earlier entry's group. Not a failure --
     * counting it as one is what made three successful files report "3
     * fetched, 2 failed". */
    if (idx == (size_t)-1) continue;
    if (idx >= YAME_ASSETS_N) { ++bad; continue; }

    /* Collect every file picked from this same directory, then mark them
     * consumed so the group is only fetched once. */
    const char *names[PICK_MAX];
    size_t n_names = 0;
    for (size_t j = i; j < b->n_pick && n_names < PICK_MAX; ++j) {
      if (b->pick[j].asset != idx || !b->pick[j].name[0]) continue;
      names[n_names++] = b->pick[j].name;
      if (j != i) b->pick[j].asset = (size_t)-1;   /* consumed */
    }
    if (!n_names) continue;

    const unit_t *a = &YAME_ASSETS[idx];

    fprog_t fp;
    memset(&fp, 0, sizeof(fp));
    fp.asset = idx;

    yame_fetch_opt_t opt = {0};
    opt.force = b->force;
    if (in_widget) {
      char what[600], howmany[64];
      snprintf(what, sizeof(what), "%s/%s", a->source, a->dir);
      snprintf(howmany, sizeof(howmany), "%zu file%s", n_names,
               n_names == 1 ? "" : "s");
      panel_row(0, yame_ui_cyan(), yame_ui_unicode() ? "⤓" : ">", what,
                howmany);
      yame_ui_panel_line(1, "");
      opt.quiet = 1;                /* the panel and the rows are the output */
      opt.on_begin = fp_begin;
      opt.on_progress = fp_progress;
      opt.on_done = fp_done;
      opt.ud = &fp;
    } else {
      fprintf(stderr, "%s/%s  -  %zu file%s\n", a->source, a->dir,
              n_names, n_names == 1 ? "" : "s");
    }

    char *err = NULL;
    if (fetch_names(a, b->root, names, n_names, &opt, &err) == 0) {
      ok += (int)n_names;
      for (size_t j = 0; j < n_names; ++j)
        for (size_t k = 0; k < a->n_files; ++k)
          if (strcmp(yame_file_name(a->files[k]), names[j]) == 0) moved += a->files[k]->size;
    } else {
      bad += (int)n_names;
      if (in_widget)
        panel_row(1, yame_ui_red(), yame_ui_cross(), "failed",
                  err ? err : "(no detail)");
      else fprintf(stderr, "  failed: %s\n", err ? err : "(no detail)");
    }
    free(err);
  }

  if (ok_out) *ok_out = ok;
  if (bad_out) *bad_out = bad;
  if (bytes_out) *bytes_out = moved;
}

/* The header's last field spans the same two things every tail holds, laid
 * out the same way, so it sits over the columns it names. */
static const char *browse_header(void) {
  static char h[128];
  char tail[64];
  snprintf(tail, sizeof(tail), TAIL_FMT, "TAG", "IN STORE");
  snprintf(h, sizeof(h), "TARGET\tKIND\t%10s\t%s", "ROWS", tail);
  return h;
}

/* Digits grouped in threes, 866553 -> "866,553": a row count is compared by
 * eye against `yame info`, and ungrouped it is easy to misread by a digit. */
static void commify(uint64_t v, char *out, size_t n) {
  char raw[32];
  int len = snprintf(raw, sizeof raw, "%" PRIu64, v);
  size_t o = 0;
  for (int i = 0; i < len && o + 2 < n; ++i) {
    if (i && (len - i) % 3 == 0) out[o++] = ',';
    out[o++] = raw[i];
  }
  out[o] = '\0';
}

/* Fill in one root row: a species heading, or a unit. Every count in it is a
 * function of the store, so this is also the refresh after a fetch. */
static void root_row(browse_t *b, size_t i, const char *group,
                     const char *unit) {
  size_t total, have;
  char note[64], line[256];

  if (unit) {
    unit_counts(b->root, unit, "", 1, &total, &have);
    counts_note(total, have, note, sizeof(note));
    /* A conflicted directory cannot be written into at all, so the gauge
     * would be describing a fetch that is not on offer. Say what is actually
     * wrong, in the column the eye is already in. */
    size_t stale = unit_stale_count(b->root, unit);
    if (stale) snprintf(note, sizeof(note), "%zu stale: -f", stale);
    /* The row space's size: a knowledgebase is only usable against a query
     * in the same row space, and this is the number to compare with `yame
     * info` on the query before fetching anything. */
    char rows[32] = "";
    uint64_t nr = yame_ref_rows_by_name(unit);
    if (nr) {                       /* right-aligned, so sizes compare by eye */
      char c[24];
      commify(nr, c, sizeof c);
      snprintf(rows, sizeof rows, "%10s", c);
    }
    snprintf(line, sizeof(line), "%s\t%s\t%s\t" TAIL_FMT, unit,
             unit_is_array(unit) ? "array" : "genome", rows, unit_tag(unit), note);
  } else {
    char upper[64];
    size_t k = 0;
    for (; group[k] && k + 1 < sizeof(upper); ++k)
      upper[k] = (group[k] >= 'a' && group[k] <= 'z')
                   ? (char)(group[k] - 'a' + 'A') : group[k];
    upper[k] = '\0';

    /* No counts on a heading: it is a label for the rows under it, each of
     * which reports its own, and a third number in the same column only
     * invited adding them up. */
    group_counts(b->root, group, &total, &have);
    snprintf(line, sizeof(line), "%s%s\t\t\t", group_mark(), upper);
  }

  free(b->roots[i]);
  b->roots[i] = wzstrdup(line);
  b->styles[i] = count_style(total, have);
}

/* Is this unit offered? Everything is, unless a caller narrowed the list. */
static int unit_offered(const browse_t *b, const char *unit) {
  if (!b->units) return 1;
  for (size_t i = 0; i < b->n_units; ++i)
    if (b->units[i] && strcmp(b->units[i], unit) == 0) return 1;
  return 0;
}

/* The offered units of one group, and so whether its heading is shown: a
 * heading over nothing would be a row that opens onto nothing. */
static size_t offered_units(const browse_t *b, const char *group,
                            char units[][256], size_t cap) {
  char all[32][256];
  size_t n = group_units(group, all, 32 < cap ? 32 : cap), k = 0;
  for (size_t j = 0; j < n; ++j)
    if (unit_offered(b, all[j])) memcpy(units[k++], all[j], 256);
  return k;
}

/* Lay out the root rows -- a heading per group, then its units -- into
 * `branch` (headings never open) and return how many. Building and
 * refreshing both walk this, so a narrowed list stays aligned with the
 * arrays the widget holds. */
static size_t layout_roots(const browse_t *b, unsigned char *branch, size_t cap) {
  char groups[16][64];
  size_t ng = all_groups(groups, 16), n = 0;
  for (size_t g = 0; g < ng && n < cap; ++g) {
    char units[32][256];
    size_t nu = offered_units(b, groups[g], units, 32);
    if (!nu) continue;
    branch[n++] = 0;
    for (size_t j = 0; j < nu && n < cap; ++j) branch[n++] = 1;
  }
  return n;
}

/* Re-read every root row from the store. */
static void refresh_roots(browse_t *b) {
  char groups[16][64];
  size_t ng = all_groups(groups, 16), i = 0;
  for (size_t g = 0; g < ng && i < b->n_roots; ++g) {
    char units[32][256];
    size_t nu = offered_units(b, groups[g], units, 32);
    if (!nu) continue;
    root_row(b, i++, groups[g], NULL);
    for (size_t j = 0; j < nu && i < b->n_roots; ++j)
      root_row(b, i++, groups[g], units[j]);
  }
}

/* Does a file answer a preselect list? Each comma-separated name matches a
 * file named in full or its set name -- the part before the first dot --
 * ignoring case: the same rule `-m ChromHMM` resolves by, so a name means
 * the same thing checked in the browser as typed on the command line. */
static int preselect_match(const char *fname, const char *names) {
  const char *dot = strchr(fname, '.');
  size_t setlen = dot ? (size_t)(dot - fname) : strlen(fname);
  for (const char *p = names; p && *p; ) {
    const char *comma = strchr(p, ',');
    size_t len = comma ? (size_t)(comma - p) : strlen(p);
    while (len && *p == ' ') { ++p; --len; }
    if (len && ((strlen(fname) == len && strncasecmp(fname, p, len) == 0) ||
                (setlen == len && strncasecmp(fname, p, len) == 0)))
      return 1;
    p = comma ? comma + 1 : NULL;
  }
  return 0;
}

/* The tree asks this of every selectable file under the opened unit. A file
 * row's key is "<registry index>|<file name>". */
static int bx_preselect(void *ctx, const char *path, const char *key) {
  const browse_t *b = ctx;
  (void) path;
  if (!key) return 0;
  for (size_t i = 0; i < b->n_pre_keys; ++i)
    if (strcmp(b->pre_keys[i], key) == 0) return 1;
  const char *bar = strchr(key, '|');
  return bar && b->preselect && preselect_match(bar + 1, b->preselect);
}

/* Before the browser opens over a store holding stale files: say which,
 * directory by directory, and offer to replace them right here. Asked
 * once, on the terminal the person is looking at -- a report printed on
 * stderr a moment before the first frame was cleared by it, and one after
 * the browser closed described a store they had already left. `n` opens
 * the browser with the rows marked stale; `y` fetches exactly those files
 * with -f, then opens it current. Returns 0 to go on to the browser, -1 if
 * the replacement failed (the browser still opens; the rows say so). */
static int stale_dialog(const char *dopt, const yame_fetch_opt_t *opt) {
  const yame_asset_file_t *want[512];
  size_t n = yame_store_stale(cfg_, dopt, want, 512);
  if (!n) return 0;
  if (n > 512) n = 512;

  size_t n_dirs = 0;
  uint64_t bytes = 0;
  for (size_t i = 0; i < n; ++i) {
    bytes += want[i]->size;
    if (i == 0 || yame_file_dirlen(want[i]) != yame_file_dirlen(want[i-1]) ||
        strncmp(want[i]->store_path, want[i-1]->store_path, yame_file_dirlen(want[i])) != 0)
      ++n_dirs;
  }
  fprintf(stderr, "%s%s fetch: %zu director%s hold%s files from an earlier release "
          "than this build pins%s\n", yame_ui_bold(), TOOL, n_dirs,
          n_dirs == 1 ? "y" : "ies", n_dirs == 1 ? "s" : "", yame_ui_reset());
  for (size_t i = 0; i < n; ++i) {
    size_t dl = yame_file_dirlen(want[i]);
    if (i == 0 || dl != yame_file_dirlen(want[i-1]) ||
        strncmp(want[i]->store_path, want[i-1]->store_path, dl) != 0) {
      char dir[4096], src[160];
      memcpy(dir, want[i]->store_path, dl); dir[dl] = '\0';
      yame_key_source(want[i]->key, src, sizeof src);
      fprintf(stderr, "  %-12s %s(%s)%s\n", dir, yame_ui_dim(),
              strncmp(src, "hf:", 3) == 0 ? src + 3 : src, yame_ui_reset());
    }
    fprintf(stderr, "      %s\n", yame_file_name(want[i]));
  }
  char sz[32], q[96];
  human_size(bytes, sz, sizeof sz);
  snprintf(q, sizeof q, "Replace them now? (%s)", sz);
  if (!yame_ui_confirm(q, 0)) { fputc('\n', stderr); return 0; }

  yame_fetch_opt_t o = *opt;
  o.force = 1;
  char root[4096], *err = NULL;
  yame_assets_root(dopt, cfg_->tool_env, root, sizeof root);
  int rc = yame_assets_fetch_files(cfg_, root, want, n, &o, &err);
  if (rc != 0) fprintf(stderr, "%s fetch: %s\n", TOOL, err ? err : "replacement failed");
  free(err);
  fputc('\n', stderr);
  return rc == 0 ? 0 : -1;
}

/* The catalogue as a browsable tree: species, then platform or build, then
 * what each publishes. Returns 0 when it ran, -1 when the terminal cannot
 * host it and the caller should print the list instead. */
static int browse_catalog(const char *dopt, int force, const char *open_unit,
                          char **pre_keys, size_t n_pre_keys) {
  enum { MAXROOT = 64 };
  static char *roots[MAXROOT];
  static unsigned char styles[MAXROOT], branch[MAXROOT];
  size_t n_roots = 0;
  static browse_t b;

  memset(&b, 0, sizeof(b));
  b.force = force;
  b.pre_keys = pre_keys;
  b.n_pre_keys = n_pre_keys;
  yame_assets_root(dopt, cfg_->tool_env, b.root, sizeof(b.root));

  /* Every unit at the top level, with its species as a heading above it: a
   * label to read past, not a level to open. */
  memset(roots, 0, sizeof(roots));
  b.roots = roots;
  b.styles = styles;

  n_roots = layout_roots(&b, branch, MAXROOT);   /* a heading never opens */
  b.n_roots = n_roots;
  refresh_roots(&b);

  /* The title says where the store is and what decided that: the tool's own
   * variable, a sibling tool's (inherited), -d, or the default. "yame fetch"
   * in front of it said nothing the person did not already know. */
  static char title[4200];
  {
    const char *var = dopt ? NULL : yame_assets_root_env(cfg_->tool_env);
    const char *own = cfg_->tool_env ? cfg_->tool_env : "YAME_DATA_HOME";
    if (dopt)          snprintf(title, sizeof(title), "-d: %s", b.root);
    else if (!var)     snprintf(title, sizeof(title), "%s (default): %s", own, b.root);
    else if (strcmp(var, own) == 0)
                       snprintf(title, sizeof(title), "%s: %s", var, b.root);
    else               snprintf(title, sizeof(title), "%s (inherited): %s", var, b.root);
  }

  yame_ui_tree_t spec;
  memset(&spec, 0, sizeof(spec));
  spec.title       = title;
  spec.header      = browse_header();
  spec.roots       = roots;
  spec.root_styles = styles;
  spec.root_branch = branch;
  spec.n_roots     = n_roots;
  spec.expand      = bx_expand;
  spec.detail      = bx_detail;
  spec.detail_key  = 'i';
  spec.detail_verb = "info";
  spec.facets      = bx_facets;
  spec.actions[0].key    = 'f';
  spec.actions[0].verb   = "fetch";
  spec.actions[0].accept = bx_accept;
  spec.actions[0].commit = bx_commit;   /* non-NULL: the tree stays open */
  spec.n_actions   = 1;
  /* nothing to ask for if it is here -- unless -f, which re-fetches it */
  spec.have_selectable = force;
  if (open_unit) {
    spec.open_root = open_unit;
    if (n_pre_keys) spec.preselect = bx_preselect;
  }
  spec.on_key      = bx_on_key;
  spec.hint        = "d        change the store";
  b.title = title; b.title_sz = sizeof title;
  spec.ctx         = &b;

  return yame_ui_tree(&spec) < 0 ? -1 : 0;
}

/**
 * Browse the catalogue and hand back what was chosen, fetching it if needed.
 *
 * The same screen as `yame fetch`, opened at the row space the caller cares
 * about, with one extra verb: `u` ends the session and returns the selection.
 * So a command that needs a knowledgebase can offer "show me what there is"
 * instead of requiring a path to a file the user has not downloaded yet.
 *
 * Returns the number of paths (malloc'd, caller frees both the strings and
 * the array), or 0 if nothing was chosen or the terminal cannot host a tree.
 */
size_t yame_browse_pick(const yame_fetch_cfg_t *cfg, const char *open_unit, char ***out_paths) {
  yame_pick_opt_t o;
  memset(&o, 0, sizeof o);
  o.open_unit = open_unit;
  return yame_browse_pick_opt(cfg, &o, out_paths);
}

size_t yame_browse_pick_opt(const yame_fetch_cfg_t *cfg, const yame_pick_opt_t *opt,
                            char ***out_paths) {
  cfg_ = cfg;
  build_units();
  enum { MAXROOT = 64 };
  static char *roots[MAXROOT];
  static unsigned char styles[MAXROOT], branch[MAXROOT];
  size_t n_roots = 0;
  static browse_t b;

  *out_paths = NULL;
  memset(&b, 0, sizeof(b));
  yame_assets_root(opt->store, cfg_->tool_env, b.root, sizeof(b.root));
  b.units = opt->units;
  b.n_units = opt->n_units;
  b.preselect = (opt->preselect && *opt->preselect) ? opt->preselect : NULL;
  b.offer = (opt->offer && *opt->offer) ? opt->offer : NULL;

  memset(roots, 0, sizeof(roots));
  b.roots = roots;
  b.styles = styles;
  n_roots = layout_roots(&b, branch, MAXROOT);
  if (!n_roots) return 0;                       /* the narrowing left nothing */
  b.n_roots = n_roots;
  refresh_roots(&b);

  char vkey = opt->verb_key ? opt->verb_key : 'u';
  static char title[4200];
  if (opt->title) snprintf(title, sizeof(title), "%s", opt->title);
  else snprintf(title, sizeof(title), "%s  %s   choose, then %c",
                cfg_->tool ? cfg_->tool : "yame", yame_ui_bullet(), vkey);

  yame_ui_tree_t spec;
  memset(&spec, 0, sizeof(spec));
  spec.title       = title;
  spec.header      = browse_header();
  spec.roots       = roots;
  spec.root_styles = styles;
  spec.root_branch = branch;
  spec.n_roots     = n_roots;
  spec.expand      = bx_expand;
  spec.detail      = bx_detail;
  spec.detail_key  = 'i';
  spec.detail_verb = "info";
  spec.facets      = bx_facets;
  spec.open_root   = opt->open_unit;    /* start where the caller is working */
  if (opt->open_unit && b.preselect) spec.preselect = bx_preselect;
  spec.actions[0].key    = 'f';
  spec.actions[0].verb   = "fetch";
  spec.actions[0].accept = bx_accept;
  spec.actions[0].commit = bx_commit;   /* stays open */
  spec.actions[1].key    = vkey;
  spec.actions[1].verb   = opt->verb ? opt->verb : "use";
  spec.actions[1].accept = bx_choose;
  spec.actions[1].commit = NULL;        /* ends the session, returns the pick */
  spec.n_actions   = 2;
  /* Something already in the store is exactly what a caller wants to use, so
   * unlike a fetch it stays selectable. */
  spec.have_selectable = 1;
  spec.on_key      = bx_on_key;
  spec.hint        = "d        change the store";
  b.title = title; b.title_sz = sizeof title;
  spec.ctx         = &b;

  if (yame_ui_tree(&spec) != 2 || !b.n_chosen) return 0;

  if (b.n_chosen_dropped)
    fprintf(stderr, "yame: %zu chosen file%s did not fit and %s ignored.\n",
            b.n_chosen_dropped, b.n_chosen_dropped == 1 ? "" : "s",
            b.n_chosen_dropped == 1 ? "was" : "were");

  /* Anything chosen but not downloaded is fetched now, on the normal screen:
   * the widget is gone, so this is an ordinary transfer with ordinary output. */
  int ok = 0, bad = 0;
  fetch_picks(&b, 0, &ok, &bad, NULL);

  char **paths = calloc(b.n_chosen, sizeof(char *));
  if (!paths) return 0;
  size_t n = 0;
  for (size_t i = 0; i < b.n_chosen; ++i) {
    size_t idx = b.chosen_asset[i];
    if (idx >= YAME_ASSETS_N) continue;
    char dir[4096], path[4096];
    if (yame_assets_join(dir, sizeof(dir), b.root, YAME_ASSETS[idx].dir) != 0)
      continue;
    if (yame_assets_join(path, sizeof(path), dir, b.chosen_name[i]) != 0) continue;
    if (!yame_assets_is_file(path)) continue;   /* the fetch did not land it */
    paths[n++] = wzstrdup(path);
  }
  if (!n) { free(paths); return 0; }
  *out_paths = paths;
  return n;
}

/* Progress for the one-shot command line.
 *
 * One line, repainted in place: which file of how many, a bar, and the bytes
 * so far. Only on a terminal -- the '\r' that makes it an indicator turns
 * into hundreds of duplicated lines in a log or a CI transcript, so there it
 * falls back to one finished line per file. */
static struct {
  char     name[56];
  uint64_t total, moved;      /* moved: bytes actually transferred so far */
  size_t   idx, n, failed;
  int      tty;
} PROG;

static void prog_bar(uint64_t now, uint64_t total) {
  const int cells = 22;
  int uni = yame_ui_unicode();
  int on = (total && now <= total) ? (int)((now * cells) / total) : 0;
  int pct = total ? (int)((now * 100) / total) : 0;

  /* Filled cells in colour, the remainder dimmed, so the boundary reads at a
   * glance instead of having to be counted. */
  char bar[256]; size_t o = 0;
  o += (size_t)snprintf(bar + o, sizeof(bar) - o, "%s", yame_ui_cyan());
  for (int i = 0; i < cells && o + 24 < sizeof(bar); ++i) {
    if (i == on) o += (size_t)snprintf(bar + o, sizeof(bar) - o, "%s%s",
                                       yame_ui_reset(), yame_ui_dim());
    o += (size_t)snprintf(bar + o, sizeof(bar) - o, "%s",
                          i < on ? (uni ? "\u2501" : "=") : (uni ? "\u2501" : "-"));
  }
  snprintf(bar + o, sizeof(bar) - o, "%s", yame_ui_reset());

  char a[32], b[32];
  human_size(now, a, sizeof(a));
  human_size(total, b, sizeof(b));

  /* \033[K clears to end of line: a shorter repaint would otherwise leave the
   * tail of the longer one behind it. */
  fprintf(stderr, "\r\033[K  %s[%zu/%zu]%s %-26.26s %s %s%3d%%%s %8s/%s",
          yame_ui_dim(), PROG.idx, PROG.n, yame_ui_reset(), PROG.name, bar,
          yame_ui_bold(), pct, yame_ui_reset(), a, total ? b : "?");
  fflush(stderr);
}

static void prog_begin(void *ud, const char *name, uint64_t total) {
  (void)ud;
  snprintf(PROG.name, sizeof(PROG.name), "%s", name);
  PROG.total = total;
  if (PROG.idx < PROG.n) PROG.idx++;
  if (PROG.tty && total) prog_bar(0, total);   /* size unknown yet: wait for the first byte */
}

static void prog_progress(void *ud, uint64_t now, uint64_t total) {
  (void)ud;
  uint64_t t = total ? total : PROG.total;
  /* Until the response headers arrive the size is unknown, and a bar with no
   * denominator is worse than no bar -- it reads as stalled at 0%. */
  if (PROG.tty && t) prog_bar(now, t);
}

static void prog_done(void *ud, const char *name, uint64_t bytes, int ok) {
  (void)ud;
  if (ok) PROG.moved += bytes;
  else    PROG.failed++;

  /* On a terminal the one line is reused: the next file paints over it, and
   * the summary replaces it at the end. Scrolling a name per file says
   * nothing the summary will not, and buries the line that is still moving.
   *
   * Piped, there is nothing to repaint onto, and a log wants the record --
   * so there each file settles on its own line. */
  if (PROG.tty) return;

  char hs[32];
  human_size(bytes, hs, sizeof(hs));
  fprintf(stderr, "  [%zu/%zu] %-40.40s %9s\n", PROG.idx, PROG.n, name,
          ok ? hs : "failed");
  fflush(stderr);
}

/* Fetch one catalogued entry: the files that entry declares, which for almost
 * every row is the whole upstream directory.
 *
 * Not fetch_subtree(), because an entry and an upstream directory are not
 * always the same set. One published directory can back more than one entry --
 * the methscope bundles are one flat repo holding both hg38 and mm10 models,
 * split here so each lands under the genome it belongs to. Taking the file
 * list from the entry rather than the manifest is what keeps that honest;
 * where the two coincide this is exactly what it did before. */
/**
 * Does one file survive the selection?
 *
 * Two independent narrowings share this: -g, which matches facets, and a spec
 * that named a file rather than a directory. Keeping them in one predicate is
 * what stops the plan, the listing and the transfer from disagreeing about
 * what is about to move -- they each used to inline the -g test.
 */
/* Where a file would land, and whether it is already there: the store's
 * layout, or the current directory under -c. */
static int sel_present(const char *root, const unit_t *a,
                       const char *name, int here) {
  if (!here) return file_present(root, a, name);
  /* In a working directory a name collision is ordinary: your own
   * human_hg38_test.cg is not the catalogue's, and the store's habit of
   * treating any file with the right name as present would quietly leave a
   * demo running on it. Here, present means the right bytes. */
  if (!yame_assets_is_file(name)) return 0;
  for (size_t i = 0; i < a->n_files; ++i)
    if (strcmp(yame_file_name(a->files[i]), name) == 0) {
      char got[65];
      if (yame_assets_sha256_file(name, got) != 0) return 0;
      return yame_assets_digest_equal(got, a->files[i]->sha256);
    }
  return 1;
}

/* One named file against one catalogue name, with its index riding along. */
static int one_file_wanted(const char *name, const char *only_file) {
  size_t l = strlen(only_file);
  if (strncmp(name, only_file, l) != 0) return 0;
  /* An index rides with its data file, the same way the browser folds it
   * into that file's row. Handing back a .cg without its .idx gives you
   * something `subset` and `split` cannot open, which is not what naming
   * the file meant. */
  if (name[l]) {
    const char *sfx = yame_assets_index_suffix(name);
    if (!sfx || strcmp(name + l, sfx) != 0) return 0;
  }
  return 1;
}

/* `only` is a LIST, because one directory can be named file by file:
 * `fetch hg38/models/a.clfx hg38/models/b.updecx` resolves to one entry with
 * two names. It held a single name until v1.44, and the second name was
 * silently dropped as a duplicate of the first. Empty list: the whole
 * directory. */
static int file_wanted(const unit_t *a, const char *name,
                       const char *filter,
                       const char *const *only, size_t n_only) {
  if (n_only) {
    int hit = 0;
    for (size_t i = 0; i < n_only && !hit; ++i)
      if (only[i] && one_file_wanted(name, only[i])) hit = 1;
    if (!hit) return 0;
  }
  if (filter) {
    char f[1024];
    file_facets(a, name, f, sizeof(f));
    if (!facets_match(f, filter)) return 0;
  }
  return 1;
}

/**
 * -c: the files themselves, in the current directory.
 *
 * A demo wants ./human_hg38_test.cg, not a store path -- and the store's
 * fetch writes SHA256SUMS verbatim beside what it takes, which is right for a
 * directory that has to stay checkable and wrong for someone's working
 * directory, where it would be litter and would collide between units.
 *
 * Verification is not weakened by skipping the manifest: the per-file digest
 * compiled into the registry is the one the manifest would have supplied, and
 * it is what the anchor made trustworthy at build time.
 */
static int fetch_entry_here(const unit_t *a, const char *filter,
                            const char *const *only, size_t n_only,
                            const yame_fetch_opt_t *opt, char **err) {
  for (size_t i = 0; i < a->n_files; ++i) {
    const char *name = yame_file_name(a->files[i]);
    if (!file_wanted(a, name, filter, only, n_only)) continue;
    int got = 0;
    if (yame_assets_download_verify(a->files[i]->url, a->files[i]->sha256, name,
                                    opt, &got, err) != 0)
      return -1;
  }
  return 0;
}

/* Fetch the named files of one directory, each verified against the digest
 * the registry compiles in, and rewrite the directory's manifest. */
static int fetch_names(const unit_t *a, const char *store_root,
                       const char *const *names, size_t n_names,
                       const yame_fetch_opt_t *opt, char **err) {
  if (!n_names) return 0;
  const yame_asset_file_t **want = wzmalloc(n_names * sizeof(*want));
  size_t n = 0;
  for (size_t i = 0; i < n_names; ++i) {
    const yame_asset_file_t *f = file_of(a, names[i]);
    if (f) want[n++] = f;
  }
  int rc = yame_assets_fetch_files(cfg_, store_root, want, n, opt, err);
  free(want);
  /* The manifest just changed on disk; the next presence question has to read
   * the new one. */
  dir_sums_forget();
  return rc;
}

static int fetch_entry(const unit_t *a, const char *store_root,
                       const char *filter,
                       const char *const *only, size_t n_only,
                       const yame_fetch_opt_t *opt, char **err) {
  const char **names = malloc((a->n_files ? a->n_files : 1) * sizeof(*names));
  if (!names) return -1;
  size_t n_names = 0;
  for (size_t i = 0; i < a->n_files; ++i) {
    const char *name = yame_file_name(a->files[i]);
    if (!file_wanted(a, name, filter, only, n_only)) continue;
    names[n_names++] = name;
  }
  if (!n_names) { free(names); return 0; }
  int rc = fetch_names(a, store_root, names, n_names, opt, err);
  free(names);
  return rc;
}

/* One selected directory, and what of it: the whole unit when `only` is
 * empty, the files a name picked out otherwise. Per entry rather than per
 * command, because `fetch a/x.cg b/y.cg` restricts each of the two
 * differently -- and each may carry its own @tag. */

/**
 * Resolve one name onto selections, appending to `out`.
 *
 * Naming the same directory twice is not an error, and a whole-unit selection
 * absorbs a file-level one: `fetch hg38/data hg38/data/x.cg` takes the
 * directory once, not the directory and then the file again.
 */
static int resolve_spec(const char *arg,
                        sel_t *out, size_t *n_out, size_t cap, int quiet) {
  char spec[512];
  if (snprintf(spec, sizeof(spec), "%s", arg) >= (int)sizeof(spec)) {
    fprintf(stderr, "%s fetch: " "target name too long.\n", TOOL);
    return 1;
  }

  const unit_t *hits[64];
  const char *only = NULL;

  /* Whatever the browser showed you is what you can type: a row is named
   * <unit>[/<folder>] -- hg38/data, EPIC/KYCG, hg38 -- and that is the
   * store path, so the browser can always tell you the command for the
   * thing you are looking at. (The registry's <source>/<target> was once a
   * second spelling; it was accepted through 1.51 and is gone.)
   *
   * A name is one store directory's own files -- "hg38" the annotation,
   * "hg38/data" the datasets -- or, with -R, that directory and every one
   * beneath it. */
  size_t n_sel = collect_scope(spec, hits, 64);

  /* A file name, bare or with a scope in front. Someone copying a command
   * out of the documentation types the file it names; making them work out
   * that human_hg38_test.cg lives in hg38/data is a lookup this can do. A
   * name claimed by several directories (cpg_nocontig.cr is in three) is
   * reported rather than guessed at. */
  if (!n_sel) {
    const char *cut = strrchr(spec, '/');
    const char *fname = cut ? cut + 1 : spec;
    char head[256] = "";
    if (cut && (size_t)(cut - spec) < sizeof(head))
      memcpy(head, spec, (size_t)(cut - spec));

    const unit_t *hit[16];
    char hitpath[16][544];   /* a 264-byte browser path plus a file name */
    size_t n_hit = 0, n_claim = 0;   /* shown, and how many there really are */
    for (size_t i = 0; i < YAME_ASSETS_N; ++i) {
      for (size_t j = 0; j < YAME_ASSETS[i].n_files; ++j) {
        if (strcmp(yame_file_name(YAME_ASSETS[i].files[j]), fname) != 0) continue;
        const char *dir = YAME_ASSETS[i].dir;
        if (head[0] && strcasecmp(dir, head) != 0) break;
        ++n_claim;
        if (n_hit < 16) {
          hit[n_hit] = &YAME_ASSETS[i];
          snprintf(hitpath[n_hit], sizeof(hitpath[0]), "%.263s/%.270s",
                   dir, fname);
          ++n_hit;
        }
        break;
      }
    }
    if (n_claim == 1) {
      hits[0] = hit[0];
      n_sel = 1;
      /* The registry string outlives this resolver. Pointing back into `arg`
       * retained an @tag suffix, while pointing into local `spec` would
       * dangle on return. */
      only = yame_file_name(file_of(hit[0], fname));
    }
    else if (n_claim > 1) {
      fprintf(stderr, "%s fetch: " "%zu directories publish a file called "
                      "%s. Name one:\n", TOOL, n_claim, fname);
      for (size_t i = 0; i < n_hit; ++i)
        fprintf(stderr, "  %s\n", hitpath[i]);
      if (n_claim > n_hit)
        fprintf(stderr, "  %s and %zu more\n",
                yame_ui_unicode() ? "\u2026" : "...", n_claim - n_hit);
      return 1;
    }
  }

  if (!n_sel) {
    if (quiet) return 1;              /* the caller has another reading to try */
    fprintf(stderr,
            "%s fetch: nothing in the catalogue is called %s.\n"
            "  Name it the way the browser shows it: hg38, hg38/KYCG,\n"
            "  hg38/data, EPIC. A name is one store directory's own files;\n"
            "  name a subdirectory to add it: %s fetch -y EPICv2 EPICv2/KYCG.\n"
            "  `%s fetch -l` lists what there is.\n", TOOL, spec, TOOL, TOOL);
    return 1;
  }
  for (size_t i = 0; i < n_sel; ++i) {
    int dup = 0;
    for (size_t k = 0; k < *n_out; ++k)
      if (out[k].a == hits[i]) {
        dup = 1;
        /* Naming the directory itself absorbs any file already picked out of
         * it, and stays absorbing: a whole-directory selection carries an
         * empty list and a later file name does not narrow it back down. */
        if (!only) out[k].n_only = 0;
        else if (out[k].n_only) {
          /* A SECOND file from the SAME directory. Until v1.44 this was taken
           * for a duplicate of the first and dropped, so `fetch dir/a dir/b`
           * silently fetched only `a`. */
          int seen = 0;
          for (size_t m = 0; m < out[k].n_only; ++m)
            if (strcmp(out[k].only[m], only) == 0) { seen = 1; break; }
          if (!seen) {
            if (out[k].n_only >= SEL_ONLY_MAX) {
              fprintf(stderr, "%s fetch: " "more than %d files named out of "
                      "one directory; name the directory instead.\n",
                      TOOL, SEL_ONLY_MAX);
              return 1;
            }
            out[k].only[out[k].n_only++] = only;
          }
        }
        break;
      }
    if (dup) continue;
    if (*n_out >= cap) { fprintf(stderr, "%s fetch: " "too many names.\n", TOOL); return 1; }
    out[*n_out].a = hits[i];
    out[*n_out].n_only = 0;
    if (only) out[*n_out].only[out[*n_out].n_only++] = only;
    ++*n_out;
  }
  return 0;
}

/* Every name on the command line onto selections. `-l` and a real fetch take
 * the same path, so `-l <names>` prints exactly what fetching them would
 * take -- including a file-level name, which used to list nothing at all.
 * `shown` collects the names for the plan line and the errors; pass NULL when
 * there is no plan to print. */
/* The whole catalogue as a selection: every entry, no file narrowed. */
static size_t sel_all(sel_t *out, size_t cap) {
  size_t n = 0;
  for (size_t i = 0; i < YAME_ASSETS_N && n < cap; ++i) {
    out[n].a = &YAME_ASSETS[i];
    out[n].n_only = 0;
    ++n;
  }
  return n;
}

static int resolve_args(int argc, char *argv[], int first,
                        sel_t *sel, size_t *n_sel, size_t SELCAP,
                        char *shown, size_t shown_cap) {
  for (int ai = first; ai < argc; ++ai) {
    /* Commas separate names, the same way -m takes CGI,ChromHMM. Split
     * first: parsing @tag before the comma would let `a@v3,b` read the tag
     * as "v3,b" and swallow the second name.
     *
     * The whole string is still tried if any piece fails, so the syntax does
     * not forbid a comma inside a catalogue name -- such a name would fail
     * as pieces and resolve as itself. */
    int rc = 1;
    size_t before = *n_sel;
    const char *bad = NULL;           /* the piece that did not resolve */
    if (strchr(argv[ai], ',')) {
      /* Deliberately never freed: a selection's `only` points into this copy
       * and outlives the loop. One small allocation per argument. */
      char *work = strdup(argv[ai]);
      if (!work) return 1;
      rc = 0;
      for (char *save = NULL, *tok = strtok_r(work, ",", &save);
           tok && rc == 0; tok = strtok_r(NULL, ",", &save)) {
        while (*tok == ' ') ++tok;        /* "a, b" reads like -m does */
        if (!*tok) continue;              /* a trailing or doubled comma */
        rc = resolve_spec(tok, sel, n_sel, SELCAP, 1);
        if (rc != 0) bad = tok;
      }
      if (rc != 0) *n_sel = before;
    }
    if (rc != 0) {
      /* The whole string may be a name in its own right. If it is not, the
       * piece that failed is the more useful thing to name -- "nothing is
       * called nope" beats quoting the entire list back. */
      if (resolve_spec(argv[ai], sel, n_sel, SELCAP, bad != NULL) != 0) {
        if (bad) resolve_spec(bad, sel, n_sel, SELCAP, 0);
        return 1;
      }
    }
    if (shown) {
      size_t l = strlen(shown);
      snprintf(shown + l, shown_cap - l, "%s%s", l ? " " : "", argv[ai]);
    }
  }
  return 0;
}

/*
 * Move the option arguments in front of the names, so `fetch hg38/KYCG -g CGI`
 * means what `fetch -g CGI hg38/KYCG` means.
 *
 * GNU getopt does this permutation itself, which is why the form the docs
 * promise -- "options may come before or after the names" -- worked on Linux
 * and failed on macOS, where BSD getopt stops at the first non-option and the
 * rest arrive as names. `-g` then read as a catalogue name nobody has. Our own
 * macOS CI leg never caught it because no test passed an option after a name;
 * POSIXLY_CORRECT=1 reproduces the BSD behaviour on Linux, and one now does.
 *
 * `optstring` says which options take an argument, so `-g CGI` travels as a
 * pair. Everything after a bare `--` is a name, whatever it looks like.
 */
static void permute_opts(int argc, char **argv, const char *optstring) {
  char **opts = wzmalloc((size_t)argc * sizeof(char *));
  char **names = wzmalloc((size_t)argc * sizeof(char *));
  int n_opt = 0, n_name = 0, i = 1, done = 0;

  for (; i < argc; ++i) {
    char *a = argv[i];
    if (!done && a[0] == '-' && a[1] == '-' && a[2] == '\0') {
      done = 1;                     /* keep the -- itself, where it was */
      opts[n_opt++] = a;
      continue;
    }
    if (done || a[0] != '-' || a[1] == '\0') { names[n_name++] = a; continue; }

    opts[n_opt++] = a;
    /* Does the LAST letter of this cluster take an argument? -qg CGI does. */
    char last = a[strlen(a) - 1];
    const char *p = strchr(optstring, last);
    if (p && p[1] == ':' && strlen(a) == 2 && i + 1 < argc)
      opts[n_opt++] = argv[++i];
  }

  int k = 1;
  for (i = 0; i < n_opt; ++i)  argv[k++] = opts[i];
  for (i = 0; i < n_name; ++i) argv[k++] = names[i];
  free(opts); free(names);
}

int yame_fetch_main(const yame_fetch_cfg_t *cfg, int argc, char *argv[]) {
  cfg_ = cfg;
  build_units();
  optind = 1;                       /* a library entry: never trust the caller's */
  recursive_ = 0;
  const char *dopt = NULL;
  const char *url = NULL, *sha = NULL, *dest = NULL;
  int force = 0, quiet = 0, list = 0, assume_yes = 0;
  int dry_run = 0, here = 0;
  const char *filter = NULL;
  int c;

  permute_opts(argc, argv, "cd:flqu:s:o:yng:hR");
  while ((c = getopt(argc, argv, "cd:flqu:s:o:yng:hR")) >= 0) {
    switch (c) {
    case 'c': here = 1; break;
    case 'd': dopt = optarg; break;
    case 'f': force = 1; break;
    case 'R': recursive_ = 1; break;
    case 'l': list = 1; break;
    case 'g': filter = optarg; break;
    case 'n': dry_run = 1; break;
    case 'y': assume_yes = 1; break;
    case 'q': quiet = 1; break;
    case 'u': url = optarg; break;
    case 's': sha = optarg; break;
    case 'o': dest = optarg; break;
    case 'h': return usage();
    default: return usage();
    }
  }

  /* After the loop, not inside it, so `-l -d <dir>` reports presence against
   * the store that was actually asked for. Returning from the case label read
   * fine until -d could change the answer. */
  /* A store that is behind this binary says so, once, on stderr -- before
   * the listing; the browser asks instead (stale_dialog) -- and only when
   * nothing was named, since a named fetch is already the fix. One block
   * per directory, naming the files, what they came from and the command
   * that repairs it; a script sees them and does not have to know to read
   * a dir_state column. The store state comes from the same helper a
   * downstream tool calls at load time, so yame and sesame describe the
   * situation in the same words. */
  int browsing = optind >= argc && !list && !url && !sha && !dest &&
                 isatty(STDIN_FILENO) && isatty(STDOUT_FILENO);
  if (optind >= argc && !quiet && !browsing) {
    if (yame_store_report(cfg_, dopt, stderr) > 0 && !list)
      fputc('\n', stderr);
  }

  if (list) {
    sel_t lsel[64];
    size_t n_l = 0;
    if (optind >= argc) n_l = sel_all(lsel, 64);
    else if (resolve_args(argc, argv, optind, lsel, &n_l,
                          sizeof(lsel)/sizeof(lsel[0]), NULL, 0) != 0)
      return 1;
    return dump_registry(dopt, lsel, n_l, filter);
  }

  /* About to actually use the store, so this is where the one-time notice
   * about a pre-consolidation cache belongs -- not in every path that merely
   * asks where the store is. */
  {
    char nroot[4096];
    yame_assets_root(dopt, cfg_->tool_env, nroot, sizeof(nroot));
    if (!quiet) yame_assets_legacy_notice(nroot);
  }

  yame_fetch_opt_t opt = {0};
  opt.force = force;
  opt.quiet = quiet;
  if (!quiet) { opt.on_begin = prog_begin; opt.on_done = prog_done; }

  char *err = NULL;

  /* ---- single-file form ---- */
  if (url || sha || dest) {
    if (!url || !sha || !dest) {
      fprintf(stderr, "%s fetch: " "-u, -s and -o go together.\n", TOOL);
      return 1;
    }
    if (yame_assets_download_verify(url, sha, dest, &opt, NULL, &err) != 0) {
      fprintf(stderr, "%s fetch: " "%s\n", TOOL, err ? err : "failed");
      free(err);
      return 1;
    }
    free(err);
    return 0;
  }

  /* ---- no target: browse ----
   *
   * Only on a terminal. Piped or redirected, this falls back to the flat
   * listing, so a script that runs `yame fetch` never blocks on a widget
   * waiting for a keystroke that will not come. */
  if (optind >= argc) {
    sel_t lsel[64];
    size_t n_l = sel_all(lsel, sizeof(lsel)/sizeof(lsel[0]));
    if (!isatty(STDIN_FILENO) || !isatty(STDOUT_FILENO))
      return dump_registry(dopt, lsel, n_l, filter);
    if (!quiet) stale_dialog(dopt, &opt);
    if (browse_catalog(dopt, force, NULL, NULL, 0) == 0) return 0;
    return dump_registry(dopt, lsel, n_l, filter); /* no terminal for the widget */
  }

  /* ---- catalogued form: <name>[@tag] ... ---- */
  sel_t sel[64];
  size_t n_sel = 0;
  char shown[512] = "";           /* the names, for the plan and the errors */
  const size_t SELCAP = sizeof(sel)/sizeof(sel[0]);
  if (resolve_args(argc, argv, optind, sel, &n_sel, SELCAP,
                   shown, sizeof(shown)) != 0)
    return 1;


  char root[4096];
  yame_assets_root(dopt, cfg_->tool_env, root, sizeof(root));
  /* -c never touches the store, so a read-only one must not stop it. */
  if (!here && !yame_assets_root_writable(root)) {
    fprintf(stderr,
            "%s fetch: %s is not writable. A read-only shared store is fine "
            "to read from, but nothing can be fetched into it; set -d or "
            "%s to somewhere you own.\n", TOOL, root,
            cfg_->tool_env ? cfg_->tool_env : "YAME_DATA_HOME");
    return 1;
  }

  /* One plan, whether or not a filter narrowed it: how many files, from how
   * many directories, and how large. Quoting the unfiltered total next to a
   * filtered list would name a number that is not going to be transferred. */
  size_t n_files = 0, n_dirs = 0, n_have = 0;
  uint64_t total = 0;
  for (size_t i = 0; i < n_sel; ++i) {
    size_t n_here = 0;
    for (size_t j = 0; j < sel[i].a->n_files; ++j) {
      if (!file_wanted(sel[i].a, yame_file_name(sel[i].a->files[j]), filter,
                       sel[i].only, sel[i].n_only))
        continue;
      ++n_here;
      /* Already-present files are skipped unless -f, so counting their bytes
       * in the total would quote a transfer that is not going to happen. */
      if (!force && sel_present(root, sel[i].a, yame_file_name(sel[i].a->files[j]), here))
        ++n_have;
      else total += sel[i].a->files[j]->size;
    }
    if (n_here) { ++n_dirs; n_files += n_here; }
  }

  if (!n_files) {
    if (filter)
      fprintf(stderr, "%s fetch: " "nothing under %s matches -g %s.\n"
                      "  Terms are ANDed and match the file name, its source, "
                      "collection, title or upstream database.\n", TOOL, shown, filter);
    else
      fprintf(stderr, "%s fetch: " "%s holds no files.\n", TOOL, shown);
    return 1;
  }

  /* Nothing to move is a finished job, not an empty plan: say so plainly and
   * stop, rather than printing a size of nothing and asking to confirm it. */
  if (n_have == n_files) {
    /* The filter is a qualifier on the count, not part of the name: "hg38/KYCG
     * -g TFBS:" read as a store path the page had just finished defining. */
    fprintf(stderr, "%s: all %zu file%s%s%s already %s.\n", shown, n_files,
            n_files == 1 ? "" : "s", filter ? " matching -g " : "",
            filter ? filter : "", here ? "here" : "in the store");
    return 0;
  }

  char hs[32];
  human_size(total, hs, sizeof(hs));

  /* Say what is about to happen when the answer is not obvious: a filter, or
   * a name that reached more than one directory. A folder can be very large
   * -- hg38 reaches the 2.9 GB whole-genome decoder -- and the browser has
   * always confirmed before transferring, so the shorter name should not be
   * the more dangerous one. */
  /* Name the files, biggest first, rather than only counting them: the whole
   * question before a fetch is what is about to land, and one 2.8 GB decoder
   * among six small files is the thing worth seeing. Capped, because -g array
   * matches 225 files and a prompt nobody reads is not a confirmation. */
  {
    fprintf(stderr, "%s: %zu file%s%s%s in %zu director%s", shown,
            n_files, n_files == 1 ? "" : "s",
            filter ? " matching -g " : "", filter ? filter : "",
            n_dirs, n_dirs == 1 ? "y" : "ies");
    if (n_have)
      fprintf(stderr, " -- %zu already %s, %zu to fetch",
              n_have, here ? "here" : "in the store", n_files - n_have);
    fprintf(stderr, ", %s\n", hs);

    struct { const char *name; uint64_t size; } *v =
        wzmalloc(n_files * sizeof(*v));
    if (v) {
      size_t k = 0;
      for (size_t i = 0; i < n_sel; ++i)
        for (size_t j = 0; j < sel[i].a->n_files && k < n_files; ++j) {
          if (!file_wanted(sel[i].a, yame_file_name(sel[i].a->files[j]), filter,
                       sel[i].only, sel[i].n_only))
            continue;
          if (!force && sel_present(root, sel[i].a, yame_file_name(sel[i].a->files[j]), here))
            continue;                 /* listing what will move, not what is */
          v[k].name = yame_file_name(sel[i].a->files[j]);
          v[k].size = sel[i].a->files[j]->size;
          ++k;
        }
      /* Insertion sort: k is at most a few hundred, and this keeps the
       * comparator next to the thing it orders. */
      for (size_t a = 1; a < k; ++a) {
        typeof(v[0]) t = v[a];
        size_t b = a;
        while (b && v[b - 1].size < t.size) { v[b] = v[b - 1]; --b; }
        v[b] = t;
      }
      /* A fetch is about to show a live progress line per file, so its plan
       * only needs to say what kind of thing is coming -- three names on one
       * line. -n has nothing following it, so there the list is the answer
       * and it gets room. */
      const size_t SHOW = dry_run ? 10 : 3;
      if (dry_run) {
        for (size_t i = 0; i < k && i < SHOW; ++i) {
          char one[32];
          human_size(v[i].size, one, sizeof(one));
          fprintf(stderr, "  %-44.44s %9s\n", v[i].name, one);
        }
        if (k > SHOW)
          fprintf(stderr, "  %s and %zu more\n",
                  yame_ui_unicode() ? "\u2026" : "...", k - SHOW);
      } else {
        fprintf(stderr, "  %s", yame_ui_dim());
        for (size_t i = 0; i < k && i < SHOW; ++i)
          fprintf(stderr, "%s%s", i ? ", " : "", v[i].name);
        if (k > SHOW)
          fprintf(stderr, " %s and %zu more",
                  yame_ui_unicode() ? "\u2026" : "...", k - SHOW);
        fprintf(stderr, "%s\n", yame_ui_reset());
      }
      free(v);
    }
  }

  PROG.n = n_files - n_have;
  PROG.idx = 0;
  PROG.tty = isatty(STDERR_FILENO);
  opt.on_progress = prog_progress;

  /* -n stops here, successfully: the plan above is the whole answer. Without
   * it, refusing non-interactively exits 1, which is right for a fetch that
   * did not happen but wrong for a question that was answered. */
  if (dry_run) return 0;

  /* An explicitly named file is its own confirmation.
   *
   * The guard exists so a short name cannot pull a directory silently --
   * `hg38` reaches 3.5 GB -- but a name that resolved to ONE file is already
   * the specific request the prompt would be asking for, and its size is in
   * the plan above either way. Without this, a documented `yame fetch
   * <file>` line is copy-pasteable only into an interactive terminal: it
   * fails in a script, in CI and for an agent, which is how every fetch line
   * in the methscope docs came to be unrunnable. Selecting a directory --
   * with or without -g -- still asks. */
  {
    int all_named = (n_sel > 0);
    for (size_t i = 0; i < n_sel; ++i)
      if (!sel[i].n_only) { all_named = 0; break; }
    if (all_named) assume_yes = 1;
  }

  /* Otherwise always ask. A fetch writes to a shared store and can be very
   * large, and the size is only knowable from the plan just printed -- so the
   * plan and the question belong together rather than the question being
   * reserved for cases someone guessed would be big. */
  if (!assume_yes) {
    if (!isatty(STDIN_FILENO)) {
      fprintf(stderr,
              "Refusing to fetch %zu file%s (%s) without confirmation. "
              "Re-run with -y.\n", n_files - n_have,
              (n_files - n_have) == 1 ? "" : "s", hs);
      return 1;
    }
    /* On a terminal the browser is the confirmation: it opens on the unit
     * with exactly the files this plan would move already checked, so the
     * list can be read, narrowed, and fetched with f -- or left with q. The
     * tree opens one unit, so names spanning several keep the prompt. -c
     * writes outside the store the browser shows, so it keeps it too. */
    if (!here && isatty(STDERR_FILENO) && yame_ui_fancy()) {
      char top[256] = "", u[256], sb[256];
      int one_unit = 1;
      for (size_t i = 0; i < n_sel && one_unit; ++i) {
        unit_of(sel[i].a, u, sizeof u, sb, sizeof sb);
        if (!top[0]) snprintf(top, sizeof top, "%s", u);
        else if (strcmp(top, u)) one_unit = 0;
      }
      if (one_unit && top[0]) {
        char **keys = wzcalloc(n_files ? n_files : 1, sizeof(char *));
        size_t nk = 0;
        for (size_t i = 0; i < n_sel; ++i)
          for (size_t j = 0; j < sel[i].a->n_files && nk < n_files; ++j) {
            const char *fn = yame_file_name(sel[i].a->files[j]);
            if (!file_wanted(sel[i].a, fn, filter, sel[i].only, sel[i].n_only)) continue;
            char k[512];
            snprintf(k, sizeof k, "%zu|%s", (size_t) (sel[i].a - YAME_ASSETS), fn);
            keys[nk++] = wzstrdup(k);
          }
        int rc = browse_catalog(dopt, force, top, keys, nk);
        for (size_t i = 0; i < nk; ++i) free(keys[i]);
        free(keys);
        if (rc == 0) return 0;
        /* no widget after all (a terminal too small, say): ask as before */
      }
    }
    fprintf(stderr, "Proceed? [y/N] ");
    int c = getchar();
    if (c != 'y' && c != 'Y') { fprintf(stderr, "nothing fetched.\n"); return 1; }
    /* The answer has served its purpose. Move back over the echoed newline
     * and the question, so the progress line takes their place instead of
     * leaving a dead prompt above it. */
    if (isatty(STDERR_FILENO)) fprintf(stderr, "\033[A\r\033[K");
  }

  for (size_t i = 0; i < n_sel; ++i) {
    const unit_t *a = sel[i].a;
    int rc = here
      ? fetch_entry_here(a, filter, sel[i].only, sel[i].n_only, &opt, &err)
      : fetch_entry(a, root, filter, sel[i].only, sel[i].n_only, &opt, &err);
    if (rc != 0) {
      fprintf(stderr, "%s fetch: " "%s\n", TOOL, err ? err : "failed");
      free(err);
      return 1;
    }
    free(err);
    err = NULL;
  }

  /* Replace the progress line rather than leaving it at 100%: one summary of
   * what moved, where the live line was. */
  if (PROG.tty) fprintf(stderr, "\r\033[K");
  if (!quiet) {
    char hs[32];
    human_size(PROG.moved, hs, sizeof(hs));
    fprintf(stderr, "  %s%s%s %zu file%s, %s%s\n",
            PROG.failed ? yame_ui_red() : yame_ui_green(),
            PROG.failed ? yame_ui_cross() : yame_ui_check(), yame_ui_reset(),
            PROG.idx - PROG.failed, (PROG.idx - PROG.failed) == 1 ? "" : "s", hs,
            PROG.failed ? " -- some failed" : "");
  }
  return 0;
}
