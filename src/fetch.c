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
 * `yame fetch` -- populate the shared asset store.
 *
 * A driver over src/assets.c, so the store can be filled without any
 * particular downstream tool installed, and so the engine they all link has a
 * way to be exercised on its own.
 *
 * Two forms:
 *   yame fetch <name>[@tag]   a directory this build pins
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
#include <unistd.h>

#include "assets.h"
#include "registry.h"
#include "assetinfo.h"
#include "yame_ui.h"

static int usage(void) {
  char root[4096];
  yame_assets_root(NULL, NULL, root, sizeof(root));

  yame_usage_head("yame fetch                              browse the catalogue");
  yame_usage_text("yame fetch [options] <name>[@tag] ...");
  yame_usage_text("yame fetch [options] -u <url> -s <sha256> -o <dest>");

  yame_usage_sec("Naming:");
  yame_usage_text("A name is what the browser shows: hg38, hg38/KYCG, hg38/data,");
  yame_usage_text("EPIC. It is a scope -- `hg38` takes the unit and everything under");
  yame_usage_text("it, `hg38/data` takes the one directory. Narrow within it with -g.");
  yame_usage_text("A file resolves too, best written out: `hg38/data/test.cg`. The");
  yame_usage_text("bare name works when only one directory publishes it.");
  yame_usage_text("Several names may be given, separated by spaces or by commas --");
  yame_usage_text("commas so a list fits an option that takes one argument. A");
  yame_usage_text("directory named twice is taken once, and naming it whole absorbs a");
  yame_usage_text("file picked out of it.");
  yame_usage_text("The registry's own <source>/<target> still resolves.");

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
  yame_usage_text("Resolved in order: -d, $YAME_DATA_HOME,");
  yame_usage_text("${XDG_DATA_HOME:-~/.local/share}/yame");
  fprintf(stderr, "  %sYAME_DATA_HOME: %s%s\n",
          yame_ui_green(), root, yame_ui_reset());

  yame_usage_sec("Options:");
  yame_usage_opt("-d <dir>", "Store root, overriding the environment.");
  yame_usage_opt("-t <tag>", "Upstream tag, overriding the pinned one. Without a digest");
  yame_usage_cont("for that tag nothing can be verified, so this needs -k.");
  yame_usage_opt("-k", "Accept an unpinned tag (no anchor check). Files are still");
  yame_usage_cont("verified against the manifest that tag publishes.");
  yame_usage_opt("-f", "Re-download what is present, and overwrite a store that");
  yame_usage_cont("was populated from a different tag.");
  yame_usage_opt("-l", "Dump the registry as TSV and exit: one row per file, with");
  yame_usage_cont("its size, digest, description and whether the store has it.");
  yame_usage_cont("Takes the same <name> and -g a fetch does, so `-l hg38 -g");
  yame_usage_cont("methscope` is the dry run for fetching exactly that.");
  yame_usage_opt("-g <a,b>", "Only files matching every term: name, source, collection,");
  yame_usage_cont("title or upstream database. `-g chromatin` inside a");
  yame_usage_cont("knowledgebase, `-g celltype` across a whole genome.");
  yame_usage_opt("-n", "Say what would be fetched and stop, successfully. The same");
  yame_usage_cont("plan a fetch prints before asking, but it exits 0, so a");
  yame_usage_cont("script can check first. -l gives the same set as TSV.");
  yame_usage_opt("-y", "Fetch a whole folder without asking. A name covering more");
  yame_usage_cont("than one directory is confirmed first, since a short name");
  yame_usage_cont("can reach a lot -- `hg38` is 3.5 GB.");
  yame_usage_opt("-q", "No progress output.");
  yame_usage_opt("-u <url>", "Single-file form: what to download.");
  yame_usage_opt("-s <sha256>", "Single-file form: the digest it must have.");
  yame_usage_opt("-o <dest>", "Single-file form: where it goes (a path, not a dir).");
  yame_usage_opt("-c", "Into the current directory rather than the store.");
  yame_usage_opt("-h", "This help.");

  yame_usage_sec("Notes:");
  yame_usage_text("* A directory records which upstream tag filled it, in the SHA256SUMS");
  yame_usage_text("  it keeps. A build pinned elsewhere refuses to overwrite it rather");
  yame_usage_text("  than start a re-download war; -f overrules that.");
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

static const yame_asset_reg_t *find_asset(const char *source, const char *target) {
  for (size_t i = 0; i < YAME_ASSETS_N; ++i)
    if (strcmp(YAME_ASSETS[i].source, source) == 0 &&
        strcmp(YAME_ASSETS[i].target, target) == 0)
      return &YAME_ASSETS[i];
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
} browse_t;

/* Defined below, next to the code that builds these rows in the first place:
 * a fetch has to put the counts back. */
static void refresh_roots(browse_t *b);

static int file_present(const char *store_root, const yame_asset_reg_t *a,
                        const char *name) {
  char dir[4096], path[4096];
  if (yame_assets_join(dir, sizeof(dir), store_root, a->store_sub) != 0) return 0;
  if (yame_assets_join(path, sizeof(path), dir, name) != 0) return 0;
  return yame_assets_is_file(path);
}

static void human_size(uint64_t b, char *out, size_t n) {
  static const char *u[] = { "B", "KB", "MB", "GB" };
  double v = (double)b; int i = 0;
  while (v >= 1024.0 && i < 3) { v /= 1024.0; ++i; }
  if (!b) snprintf(out, n, "%s", "");
  else if (i == 0) snprintf(out, n, "%.0f %s", v, u[i]);
  else snprintf(out, n, "%.1f %s", v, u[i]);
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
static void unit_of(const yame_asset_reg_t *a, char *unit, size_t nu,
                    char *sub, size_t ns) {
  /* A genome's knowledgebase is its own upstream repo, so the mapping cannot
   * come from the target alone. */
  if (strcmp(a->source, "KYCGKB") == 0) {
    snprintf(unit, nu, "%s", a->target);
    snprintf(sub, ns, "KYCG");
    return;
  }
  const char *slash = strchr(a->target, '/');
  if (slash) {
    size_t l = (size_t)(slash - a->target);
    if (l >= nu) l = nu - 1;
    memcpy(unit, a->target, l);
    unit[l] = '\0';
    snprintf(sub, ns, "%s", slash + 1);
  } else {
    snprintf(unit, nu, "%s", a->target);
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
static size_t collect_scope(const char *path, const yame_asset_reg_t **out,
                            size_t cap) {
  size_t n = 0, plen = strlen(path);
  for (size_t i = 0; i < YAME_ASSETS_N && n < cap; ++i) {
    char u[128], sb[128], full[264];
    unit_of(&YAME_ASSETS[i], u, sizeof(u), sb, sizeof(sb));
    if (sb[0]) snprintf(full, sizeof(full), "%s/%s", u, sb);
    else       snprintf(full, sizeof(full), "%s", u);
    if (strncasecmp(full, path, plen) == 0 &&
        (full[plen] == '\0' || full[plen] == '/'))
      out[n++] = &YAME_ASSETS[i];
  }
  return n;
}


/* An array platform, or a genome build? Decides which recommended list
 * applies and what the row calls itself. */
static int unit_is_array(const char *unit) {
  for (size_t i = 0; i < YAME_ASSETS_N; ++i) {
    char u[128], s[128];
    unit_of(&YAME_ASSETS[i], u, sizeof(u), s, sizeof(s));
    if (strcmp(u, unit) == 0)
      return strcmp(YAME_ASSETS[i].source, "InfiniumAnnotation") == 0;
  }
  return 0;
}

static const char *group_of(const char *unit) {
  for (size_t i = 0; UNIT_ORDER[i].unit; ++i)
    if (strcmp(UNIT_ORDER[i].unit, unit) == 0) return UNIT_ORDER[i].group;
  return "other";
}

static int unit_known(const char *unit) {
  for (size_t i = 0; i < YAME_ASSETS_N; ++i) {
    char u[128], s[128];
    unit_of(&YAME_ASSETS[i], u, sizeof(u), s, sizeof(s));
    if (strcmp(u, unit) == 0) return 1;
  }
  return 0;
}

/* The units of one group, in the table's order. */
static size_t group_units(const char *group, char units[][128], size_t cap) {
  size_t n = 0;
  if (strcmp(group, "other") != 0) {
    for (size_t i = 0; UNIT_ORDER[i].unit && n < cap; ++i)
      if (strcmp(UNIT_ORDER[i].group, group) == 0 &&
          unit_known(UNIT_ORDER[i].unit))
        snprintf(units[n++], 128, "%s", UNIT_ORDER[i].unit);
    return n;
  }
  for (size_t i = 0; i < YAME_ASSETS_N && n < cap; ++i) {
    char u[128], s[128];
    unit_of(&YAME_ASSETS[i], u, sizeof(u), s, sizeof(s));
    if (strcmp(group_of(u), "other") != 0) continue;
    size_t k = 0;
    for (; k < n; ++k) if (strcmp(units[k], u) == 0) break;
    if (k == n) snprintf(units[n++], 128, "%s", u);
  }
  return n;
}

static size_t all_groups(char groups[][64], size_t cap) {
  size_t n = 0;
  for (size_t i = 0; UNIT_ORDER[i].group && n < cap; ++i) {
    size_t k = 0;
    for (; k < n; ++k) if (strcmp(groups[k], UNIT_ORDER[i].group) == 0) break;
    if (k < n) continue;
    char units[32][128];
    if (group_units(UNIT_ORDER[i].group, units, 32))
      snprintf(groups[n++], 64, "%s", UNIT_ORDER[i].group);
  }
  char units[32][128];
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

static int has_file(const yame_asset_reg_t *a, const char *name) {
  for (size_t i = 0; i < a->n_files; ++i)
    if (strcmp(a->files[i].name, name) == 0) return 1;
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
static int is_companion(const yame_asset_reg_t *a, const char *name) {
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
static const char *companion_of(const yame_asset_reg_t *a, const char *name,
                                char *out, size_t n) {
  for (size_t k = 0; COMPANION_SFX[k]; ++k) {
    snprintf(out, n, "%s%s", name, COMPANION_SFX[k]);
    if (has_file(a, out)) return COMPANION_SFX[k];
  }
  out[0] = '\0';
  return NULL;
}

static uint64_t companion_size(const yame_asset_reg_t *a, const char *name) {
  char idxname[256];
  if (!companion_of(a, name, idxname, sizeof(idxname))) return 0;
  for (size_t i = 0; i < a->n_files; ++i)
    if (strcmp(a->files[i].name, idxname) == 0) return a->files[i].size;
  return 0;
}

/* One thing the catalogue offers: a file, plus its index when it has one. */
typedef struct {
  const yame_asset_reg_t  *a;
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
    const yame_asset_reg_t *a = &YAME_ASSETS[i];
    char u[128], s[128];
    unit_of(a, u, sizeof(u), s, sizeof(s));
    if (strcmp(u, unit) != 0) continue;

    for (size_t j = 0; j < a->n_files && n < cap; ++j) {
      const yame_asset_file_t *f = &a->files[j];
      if (is_companion(a, f->name)) continue;

      int required = is_unit_index(f->name);
      /* The index renders at the top of the unit wherever it is published. */
      const char *at = required ? "" : s;
      if (!recursive && strcmp(at, sub) != 0) continue;

      char idxname[256];
      out[n].a = a;
      out[n].f = f;
      out[n].paired = companion_of(a, f->name, idxname, sizeof(idxname));
      (void)idxname;
      out[n].required = required;
      ++n;
    }
  }
  return n;
}

/* The subdirectories directly below a unit -- in practice KYCG, but derived
 * rather than assumed. */
static size_t unit_subs(const char *unit, char subs[][128], size_t cap) {
  size_t n = 0;
  for (size_t i = 0; i < YAME_ASSETS_N && n < cap; ++i) {
    char u[128], s[128];
    unit_of(&YAME_ASSETS[i], u, sizeof(u), s, sizeof(s));
    if (strcmp(u, unit) != 0 || !s[0]) continue;

    size_t k = 0;
    for (; k < n; ++k) if (strcmp(subs[k], s) == 0) break;
    if (k == n) snprintf(subs[n++], 128, "%s", s);
  }
  return n;
}

static int ent_present(const char *store_root, const ent_t *e) {
  if (!file_present(store_root, e->a, e->f->name)) return 0;
  if (e->paired) {
    char idxname[256];
    snprintf(idxname, sizeof(idxname), "%s%s", e->f->name, e->paired);
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
 */
static int unit_pin_conflict(const char *store_root, const char *unit) {
  for (size_t i = 0; i < YAME_ASSETS_N; ++i) {
    const yame_asset_reg_t *a = &YAME_ASSETS[i];
    char u[128], s[128];
    unit_of(a, u, sizeof(u), s, sizeof(s));
    if (strcmp(u, unit) != 0) continue;

    char dir[4096];
    if (yame_assets_join(dir, sizeof(dir), store_root, a->store_sub) != 0)
      continue;
    if (yame_assets_pin_check(dir, a->anchor) == YAME_PIN_CONFLICT) return 1;
  }
  return 0;
}

static void group_counts(const char *store_root, const char *group,
                         size_t *total, size_t *have) {
  char units[32][128];
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
    char u[128], s[128];
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
    char u[128], s[128];
    unit_of(&YAME_ASSETS[i], u, sizeof(u), s, sizeof(s));
    if (strcmp(u, unit) != 0) continue;
    if (sub && strcmp(s, sub) != 0) continue;
    if (!url) url = YAME_ASSETS[i].base_url;
    else if (strcmp(url, YAME_ASSETS[i].base_url) != 0) return NULL;
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
  char unit[128];
  char sub[128];
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

static int is_group_row(const char *path) {
  const char *m = group_mark();
  return strncmp(path, m, strlen(m)) == 0;
}

/* The group a heading row names, lowercased back from its display form. */
static void group_of_row(const char *path, char *out, size_t n) {
  const char *p = path + strlen(group_mark());
  size_t i = 0;
  for (; p[i] && i + 1 < n; ++i)
    out[i] = (p[i] >= 'A' && p[i] <= 'Z') ? (char)(p[i] - 'A' + 'a') : p[i];
  out[i] = '\0';
}

/* ---- what data/assets.tsv calls a file ---- */

/**
 * The key a filename is described under.
 *
 * A knowledgebase set is named by the part before the first dot, so one row
 * covers every platform publishing it. Everything else is named by its role,
 * because "EPIC.hg38.mask.cm" and "MM285.mm39.mask.cm" are the same kind of
 * thing and describing them once is the point.
 */
static void info_key_of(const char *name, char *out, size_t n) {
  static const struct { const char *needle, *key; } roles[] = {
    { ".ordering.",    "ordering" },
    { ".coord.",       "coord" },
    { ".snp.",         "snp" },
    { ".typeI_ext.",   "typeI_ext" },
    { ".mask.",        "mask" },
    { "seqinfo.",      "seqinfo" },
    { "gaps.",         "gaps" },
    { "cytoband.",     "cytoband" },
    { "cpg_nocontig.", "cpg_nocontig" },
    { "genes.",        "genes" },
    { NULL, NULL }
  };
  for (size_t i = 0; roles[i].needle; ++i)
    if (strstr(name, roles[i].needle)) {
      snprintf(out, n, "%s", roles[i].key);
      return;
    }

  const char *dot = strchr(name, '.');
  size_t l = dot ? (size_t)(dot - name) : strlen(name);
  if (l >= n) l = n - 1;
  memcpy(out, name, l);
  out[l] = '\0';
}

/* Defined below, beside the -g filter they serve; the listing wants the same
 * matching so that -l is the dry run for a fetch. */
static void file_facets(const yame_asset_reg_t *a, const char *name,
                        char *out, size_t cap);
static int facets_match(const char *facets, const char *terms);
static size_t collect_scope(const char *path, const yame_asset_reg_t **out,
                            size_t cap);

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
 * here is cut(1), not an eye. Descriptions come from data/assets.tsv keyed by
 * role, so one row describes `mask` for every platform that publishes one.
 */
static int dump_registry(const char *dopt, const char *scope,
                         const char *filter) {
  char root[4096];
  yame_assets_root(dopt, NULL, root, sizeof(root));

  /* The same selection a fetch would make, printed instead of downloaded, so
   * `-l <name> -g <terms>` is the dry run for the command without -l. */
  const yame_asset_reg_t *sel[64];
  size_t n_sel = 0;
  if (scope) n_sel = collect_scope(scope, sel, 64);
  else for (size_t i = 0; i < YAME_ASSETS_N && n_sel < 64; ++i)
         sel[n_sel++] = &YAME_ASSETS[i];

  printf("target\tsource\ttag\tstore_path\tfile\tbytes\tsha256\tdir_state\t"
         "local\tdescription\n");

  for (size_t i = 0; i < n_sel; ++i) {
    const yame_asset_reg_t *a = sel[i];
    char dir[4096];
    yame_assets_join(dir, sizeof(dir), root, a->store_sub);

    const char *state = "-";
    switch (yame_assets_pin_check(dir, a->anchor)) {
    case YAME_PIN_MATCH:    state = "present";   break;
    case YAME_PIN_CONFLICT: state = "other_tag"; break;
    case YAME_PIN_UNKNOWN:  state = "unpinned";  break;
    default:                state = "-";         break;
    }

    for (size_t j = 0; j < a->n_files; ++j) {
      const yame_asset_file_t *f = &a->files[j];

      if (filter) {
        char facets[1024];
        file_facets(a, f->name, facets, sizeof(facets));
        if (!facets_match(facets, filter)) continue;
      }

      char path[4096];
      yame_assets_join(path, sizeof(path), dir, f->name);

      char ikey[256];
      info_key_of(f->name, ikey, sizeof(ikey));
      const yame_assetinfo_t *k = yame_assetinfo_find(ikey);

      /* A tab or newline inside a title would shift every later column, so
       * fold whitespace to single spaces on the way out. */
      char title[512];
      const char *src = k && k->title ? k->title : "-";
      size_t n = 0;
      for (; *src && n + 1 < sizeof(title); ++src) {
        char c = (*src == '\t' || *src == '\n' || *src == '\r') ? ' ' : *src;
        if (c == ' ' && n && title[n - 1] == ' ') continue;
        title[n++] = c;
      }
      title[n] = '\0';

      printf("%s\t%s\t%s\t%s\t%s\t%" PRIu64 "\t%s\t%s\t%s\t%s\n",
             a->target, a->source, a->tag, a->store_sub, f->name,
             f->size, f->sha256 ? f->sha256 : "-", state,
             yame_assets_is_file(path) ? "yes" : "no", title);
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
static void file_facets(const yame_asset_reg_t *a, const char *name,
                        char *out, size_t cap) {
  char ikey[256];
  info_key_of(name, ikey, sizeof(ikey));
  const yame_assetinfo_t *k = yame_assetinfo_find(ikey);
  char u[128], sb[128];
  unit_of(a, u, sizeof(u), sb, sizeof(sb));
  snprintf(out, cap, "%s %s %s %s %s %s %s", name, a->source, u, sb,
           k && k->collections ? k->collections : "",
           k && k->title       ? k->title       : "",
           k && k->source      ? k->source      : "");
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
  human_size(e->f->size + (e->paired ? companion_size(e->a, e->f->name) : 0),
             sz, sizeof(sz));
  snprintf(name, sizeof(name), "%s%s%s", e->f->name, e->paired ? " +" : "",
           e->paired ? e->paired + 1 : "");
  /* Both fields fixed-width: the tail is right-aligned as a whole, so a size
   * that varies in width would walk the tag column left and right. */
  snprintf(line, sizeof(line), "%s\t" TAIL_FMT, name, e->a->tag, sz);
  snprintf(key, sizeof(key), "%zu|%s", idx, e->f->name);

  out->rows[out->n]   = strdup(line);
  out->keys[out->n]   = strdup(key);
  out->styles[out->n] = (unsigned char)(here ? YAME_ROW_HAVE
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
static void bx_expand(void *ctx, const char *path, yame_ui_kids_t *out) {
  browse_t *b = ctx;
  bpath_t p;
  bpath_parse(path, &p);

  enum { CAP = 512 };
  out->rows   = calloc(CAP, sizeof(char *));
  out->keys   = calloc(CAP, sizeof(char *));
  out->styles = calloc(CAP, 1);
  out->branch = calloc(CAP, 1);
  if (!out->rows || !out->keys || !out->styles || !out->branch) return;

  static ent_t ents[CAP];
  size_t n_ents = unit_entries(p.unit, p.sub, 0, ents, CAP);

  /* The index first: it is what everything below it is addressed against. */
  for (size_t i = 0; i < n_ents && out->n < CAP; ++i)
    if (ents[i].required) emit_entry(b, out, &ents[i]);

  /* Then the subdirectories, only at the top of a unit. */
  if (!p.sub[0]) {
    char subs[8][128];
    size_t n_subs = unit_subs(p.unit, subs, 8);
    for (size_t i = 0; i < n_subs && out->n < CAP; ++i) {
      size_t total, have;
      unit_counts(b->root, p.unit, subs[i], 0, &total, &have);
      if (!total) continue;

      char note[64], line[256];
      counts_note(total, have, note, sizeof(note));
      snprintf(line, sizeof(line), SUB_ROW_FMT "\t" TAIL_FMT, subs[i], "sets",
               path_tag(p.unit, subs[i]), note);

      out->rows[out->n]   = strdup(line);
      out->keys[out->n]   = strdup(subs[i]);
      out->styles[out->n] = count_style(total, have);
      out->branch[out->n] = 1;
      ++out->n;
    }
  }

  for (size_t i = 0; i < n_ents && out->n < CAP; ++i)
    if (!ents[i].required) emit_entry(b, out, &ents[i]);
}

/* ---- the info pane ----
 *
 * A detail callback rather than a modal panel, so the arrow keys keep working
 * while it is open and the text follows the cursor. It runs on every redraw,
 * which is affordable because the lookup is a scan of a compiled-in table and
 * the wrapping is a few hundred bytes of formatting.
 *
 * Ported from kycg, which is also where data/assets.tsv started.
 */

#define INFO_MAX_LINES 64

typedef struct {
  char *line[INFO_MAX_LINES];
  int   n;
  int   cols;
} info_lay_t;

static void lay_push(info_lay_t *L, const char *s) {
  if (L->n >= INFO_MAX_LINES) return;
  L->line[L->n] = strdup(s ? s : "");
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
static void lay_provenance(info_lay_t *L, const yame_asset_reg_t *a,
                           const yame_asset_file_t *f, const char *paired) {
  char buf[1024];

  snprintf(buf, sizeof(buf), "%s @ %s", a->base_url, a->tag);
  lay_wrap(L, "upstream", buf);

  snprintf(buf, sizeof(buf), "%s/%s", a->store_sub, f ? f->name : "");
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

  if (is_group_row(path)) {
    /* A species heading. */
    char group[64];
    group_of_row(path, group, sizeof(group));

    char units[32][128], list[512] = "";
    size_t nu = group_units(group, units, 32);
    for (size_t i = 0; i < nu; ++i) {
      strncat(list, units[i], sizeof(list) - strlen(list) - 1);
      if (i + 1 < nu) strncat(list, ", ", sizeof(list) - strlen(list) - 1);
    }

    size_t total, have;
    group_counts(b->root, group, &total, &have);

    lay_head(&L, group, "genome builds and array platforms");
    lay_wrap(&L, NULL, list);
    char buf[128];
    snprintf(buf, sizeof(buf), "%zu of %zu files already in the store",
             have, total);
    lay_wrap(&L, "in store", buf);
  } else if (!bar) {
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
    const yame_asset_reg_t *a = &YAME_ASSETS[idx];

    const yame_asset_file_t *f = NULL;
    for (size_t i = 0; i < a->n_files; ++i)
      if (strcmp(a->files[i].name, name) == 0) { f = &a->files[i]; break; }

    char idxname[256];
    const char *paired = companion_of(a, name, idxname, sizeof(idxname));

    char ikey[256];
    info_key_of(name, ikey, sizeof(ikey));
    const yame_assetinfo_t *k = yame_assetinfo_find(ikey);

    if (k) {
      lay_head(&L, name, k->title);
      lay_wrap(&L, NULL, k->biology);
      lay_wrap(&L, "source", k->source);
      lay_wrap(&L, "citation", k->citation);
      lay_wrap(&L, "processing", k->processing);
    } else {
      char buf[512];
      snprintf(buf, sizeof(buf), "  %s%s%s   %s(nothing recorded about this "
               "file)%s", yame_ui_bold(), name, yame_ui_reset(),
               yame_ui_dim(), yame_ui_reset());
      lay_push(&L, buf);
      lay_push(&L, "");
    }
    lay_provenance(&L, a, f, paired);
  }

  out->rows = malloc((size_t)(L.n ? L.n : 1) * sizeof(char *));
  if (!out->rows) { lay_free(&L); return; }
  for (int i = 0; i < L.n; ++i) out->rows[i] = L.line[i];
  out->n = (size_t)L.n;
  /* Ownership moves to the widget; do not free the strings here. */
  L.n = 0;
  lay_free(&L);
}

/** Is this file part of its unit's recommended selection? */
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
    char u[128], sb[128];
    unit_of(&YAME_ASSETS[i], u, sizeof(u), sb, sizeof(sb));
    if (strcmp(u, bp.unit) != 0) continue;
    if (bp.sub[0] && strcmp(sb, bp.sub) != 0) continue;
    n += (size_t)snprintf(buf + n, sizeof(buf) - n, " %s", YAME_ASSETS[i].source);
  }
  return buf;
}

static int bx_recommend(void *ctx, const char *path, const char *key) {
  (void)ctx;
  const char *bar = key ? strchr(key, '|') : NULL;
  if (!bar) return 0;

  bpath_t p;
  bpath_parse(path, &p);
  if (!p.unit[0]) return 0;

  char ikey[256];
  info_key_of(bar + 1, ikey, sizeof(ikey));
  /* Every array platform shares one recommended list: the sets are the same
   * annotation projected onto different probe orderings. */
  return yame_assetinfo_recommended(ikey,
                                    unit_is_array(p.unit) ? "array" : p.unit);
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
  const yame_asset_reg_t *a = &YAME_ASSETS[idx];
  const char *name = bar + 1;
  pick_add(b, idx, name);

  /* The pieces that are not choices. A .cm without its .idx is unusable, and
   * everything in a unit is addressed against that unit's index -- so both
   * come along rather than being something to remember. */
  char idxname[256];
  if (companion_of(a, name, idxname, sizeof(idxname)))
    pick_add(b, idx, idxname);

  char unit[128], sub[128];
  unit_of(a, unit, sizeof(unit), sub, sizeof(sub));
  static ent_t ents[512];
  size_t n = unit_entries(unit, "", 0, ents, 512);
  for (size_t i = 0; i < n; ++i)
    if (ents[i].required)
      pick_add(b, (size_t)(ents[i].a - YAME_ASSETS), ents[i].f->name);
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
  char lbl[512];
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

  for (size_t i = 0; i < b->n_pick; ++i) {
    size_t idx = b->pick[i].asset;
    if (idx >= YAME_ASSETS_N) continue;
    ++n_files;

    const yame_asset_reg_t *a = &YAME_ASSETS[idx];
    uint64_t sz = 0;
    for (size_t j = 0; j < a->n_files; ++j)
      if (strcmp(a->files[j].name, b->pick[i].name) == 0) {
        sz = a->files[j].size;
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

  yame_ui_panel_open(4);
  panel_row(0, yame_ui_bold(), yame_ui_unicode() ? "⤓" : ">", label, value);
  if (b->n_dropped) {
    char note[128];
    snprintf(note, sizeof(note), "%zu more did not fit and will NOT be fetched",
             b->n_dropped);
    panel_row(1, yame_ui_red(), yame_ui_cross(), "selection truncated", note);
  } else if (unknown) {
    panel_row(1, NULL, " ", "", "some sizes are not published upstream");
  } else {
    yame_ui_panel_line(1, "");
  }
  return yame_ui_panel_confirm(2, "   Proceed?", 1);
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

    const yame_asset_reg_t *a = &YAME_ASSETS[idx];
    char store_sub[4096];
    if (yame_assets_join(store_sub, sizeof(store_sub), b->root, a->store_sub) != 0) {
      ++bad; continue;
    }

    fprog_t fp;
    memset(&fp, 0, sizeof(fp));
    fp.asset = idx;

    yame_fetch_opt_t opt = {0};
    opt.force = b->force;
    if (in_widget) {
      char what[256], howmany[64];
      snprintf(what, sizeof(what), "%s/%s", a->source, a->target);
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
      fprintf(stderr, "%s/%s  -  %zu file%s\n", a->source, a->target,
              n_names, n_names == 1 ? "" : "s");
    }

    char *err = NULL;
    if (yame_assets_fetch_subset(a->base_url, a->tag, a->remote_sub, store_sub,
                                 a->anchor, names, n_names, &opt, &err) == 0) {
      ok += (int)n_names;
      for (size_t j = 0; j < n_names; ++j)
        for (size_t k = 0; k < a->n_files; ++k)
          if (strcmp(a->files[k].name, names[j]) == 0) moved += a->files[k].size;
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
  snprintf(h, sizeof(h), "TARGET\tKIND\t%s", tail);
  return h;
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
    if (unit_pin_conflict(b->root, unit))
      snprintf(note, sizeof(note), "%s", "stale tag: -f");
    snprintf(line, sizeof(line), "%s\t%s\t" TAIL_FMT, unit,
             unit_is_array(unit) ? "array" : "genome", unit_tag(unit), note);
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
    snprintf(line, sizeof(line), "%s%s\t\t", group_mark(), upper);
  }

  free(b->roots[i]);
  b->roots[i] = strdup(line);
  b->styles[i] = count_style(total, have);
}

/* Re-read every root row from the store. */
static void refresh_roots(browse_t *b) {
  char groups[16][64];
  size_t ng = all_groups(groups, 16), i = 0;
  for (size_t g = 0; g < ng && i < b->n_roots; ++g) {
    root_row(b, i++, groups[g], NULL);
    char units[32][128];
    size_t nu = group_units(groups[g], units, 32);
    for (size_t j = 0; j < nu && i < b->n_roots; ++j)
      root_row(b, i++, groups[g], units[j]);
  }
}

/* The catalogue as a browsable tree: species, then platform or build, then
 * what each publishes. Returns 0 when it ran, -1 when the terminal cannot
 * host it and the caller should print the list instead. */
static int browse_catalog(const char *dopt, int force) {
  enum { MAXROOT = 64 };
  static char *roots[MAXROOT];
  static unsigned char styles[MAXROOT], branch[MAXROOT];
  size_t n_roots = 0;
  static browse_t b;

  memset(&b, 0, sizeof(b));
  b.force = force;
  yame_assets_root(dopt, NULL, b.root, sizeof(b.root));

  /* Every unit at the top level, with its species as a heading above it: a
   * label to read past, not a level to open. */
  memset(roots, 0, sizeof(roots));
  b.roots = roots;
  b.styles = styles;

  char groups[16][64];
  size_t ng = all_groups(groups, 16);
  for (size_t i = 0; i < ng && n_roots < MAXROOT; ++i) {
    branch[n_roots++] = 0;                  /* a heading never opens */
    char units[32][128];
    size_t nu = group_units(groups[i], units, 32);
    for (size_t j = 0; j < nu && n_roots < MAXROOT; ++j) branch[n_roots++] = 1;
  }
  b.n_roots = n_roots;
  refresh_roots(&b);

  static char title[4200];
  snprintf(title, sizeof(title), "yame fetch   %s   %s", yame_ui_bullet(),
           b.root);

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
  spec.recommend   = bx_recommend;
  spec.facets      = bx_facets;
  spec.actions[0].key    = 'f';
  spec.actions[0].verb   = "fetch";
  spec.actions[0].accept = bx_accept;
  spec.actions[0].commit = bx_commit;   /* non-NULL: the tree stays open */
  spec.n_actions   = 1;
  spec.have_selectable = 0;             /* nothing to ask for if it is here */
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
size_t yame_browse_pick(const char *open_unit, char ***out_paths) {
  enum { MAXROOT = 64 };
  static char *roots[MAXROOT];
  static unsigned char styles[MAXROOT], branch[MAXROOT];
  size_t n_roots = 0;
  static browse_t b;

  *out_paths = NULL;
  memset(&b, 0, sizeof(b));
  yame_assets_root(NULL, NULL, b.root, sizeof(b.root));

  memset(roots, 0, sizeof(roots));
  b.roots = roots;
  b.styles = styles;

  char groups[16][64];
  size_t ng = all_groups(groups, 16);
  for (size_t i = 0; i < ng && n_roots < MAXROOT; ++i) {
    branch[n_roots++] = 0;
    char units[32][128];
    size_t nu = group_units(groups[i], units, 32);
    for (size_t j = 0; j < nu && n_roots < MAXROOT; ++j) branch[n_roots++] = 1;
  }
  b.n_roots = n_roots;
  refresh_roots(&b);

  static char title[4200];
  snprintf(title, sizeof(title), "yame  %s   choose, then u", yame_ui_bullet());

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
  spec.recommend   = bx_recommend;
  spec.facets      = bx_facets;
  spec.open_root   = open_unit;         /* start where the caller is working */
  spec.actions[0].key    = 'f';
  spec.actions[0].verb   = "fetch";
  spec.actions[0].accept = bx_accept;
  spec.actions[0].commit = bx_commit;   /* stays open */
  spec.actions[1].key    = 'u';
  spec.actions[1].verb   = "use";
  spec.actions[1].accept = bx_choose;
  spec.actions[1].commit = NULL;        /* ends the session, returns the pick */
  spec.n_actions   = 2;
  /* Something already in the store is exactly what a caller wants to use, so
   * unlike a fetch it stays selectable. */
  spec.have_selectable = 1;
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
    if (yame_assets_join(dir, sizeof(dir), b.root, YAME_ASSETS[idx].store_sub) != 0)
      continue;
    if (yame_assets_join(path, sizeof(path), dir, b.chosen_name[i]) != 0) continue;
    if (!yame_assets_is_file(path)) continue;   /* the fetch did not land it */
    paths[n++] = strdup(path);
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
static int sel_present(const char *root, const yame_asset_reg_t *a,
                       const char *name, int here) {
  if (!here) return file_present(root, a, name);
  /* In a working directory a name collision is ordinary: your own
   * human_hg38_test.cg is not the catalogue's, and the store's habit of
   * treating any file with the right name as present would quietly leave a
   * demo running on it. Here, present means the right bytes. */
  if (!yame_assets_is_file(name)) return 0;
  for (size_t i = 0; i < a->n_files; ++i)
    if (strcmp(a->files[i].name, name) == 0) {
      char got[65];
      if (yame_assets_sha256_file(name, got) != 0) return 0;
      return yame_assets_digest_equal(got, a->files[i].sha256);
    }
  return 1;
}

static int file_wanted(const yame_asset_reg_t *a, const char *name,
                       const char *filter, const char *only_file) {
  if (only_file) {
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
static int fetch_entry_here(const yame_asset_reg_t *a, const char *tag,
                            const char *filter, const char *only_file,
                            const yame_fetch_opt_t *opt, char **err) {
  for (size_t i = 0; i < a->n_files; ++i) {
    const char *name = a->files[i].name;
    if (!file_wanted(a, name, filter, only_file)) continue;

    char url[4096];
    int n = (a->remote_sub && a->remote_sub[0])
      ? snprintf(url, sizeof(url), "%s/%s/%s/%s",
                 a->base_url, tag, a->remote_sub, name)
      : snprintf(url, sizeof(url), "%s/%s/%s", a->base_url, tag, name);
    if (n < 0 || (size_t)n >= sizeof(url)) return -1;

    int got = 0;
    if (yame_assets_download_verify(url, a->files[i].sha256, name,
                                    opt, &got, err) != 0)
      return -1;
  }
  return 0;
}

static int fetch_entry(const yame_asset_reg_t *a, const char *store_root,
                       const char *tag, const char *anchor, const char *filter,
                       const char *only_file,
                       const yame_fetch_opt_t *opt, char **err) {
  char store_sub[4096];
  if (yame_assets_join(store_sub, sizeof(store_sub), store_root, a->store_sub) != 0)
    return -1;

  const char **names = malloc(a->n_files * sizeof(*names));
  if (!names) return -1;
  size_t n_names = 0;
  for (size_t i = 0; i < a->n_files; ++i) {
    if (!file_wanted(a, a->files[i].name, filter, only_file)) continue;
    names[n_names++] = a->files[i].name;
  }
  if (!n_names) { free(names); return 0; }

  int rc = yame_assets_fetch_subset(a->base_url, tag, a->remote_sub, store_sub,
                                    anchor, names, n_names, opt, err);
  free(names);
  return rc;
}

/* One selected directory, and what of it: the whole unit when `only` is
 * NULL, one file when a name picked one out. Per entry rather than per
 * command, because `fetch a/x.cg b/y.cg` restricts each of the two
 * differently -- and each may carry its own @tag. */
typedef struct {
  const yame_asset_reg_t *a;
  const char *only;               /* NULL: the whole directory */
  const char *tag;                /* @tag on this name, else -t, else pinned */
} sel_t;

/**
 * Resolve one name onto selections, appending to `out`.
 *
 * Naming the same directory twice is not an error, and a whole-unit selection
 * absorbs a file-level one: `fetch hg38/data hg38/data/x.cg` takes the
 * directory once, not the directory and then the file again.
 */
static int resolve_spec(const char *arg, const char *tag_opt,
                        sel_t *out, size_t *n_out, size_t cap, int quiet) {
  char spec[512];
  if (snprintf(spec, sizeof(spec), "%s", arg) >= (int)sizeof(spec)) {
    fprintf(stderr, "yame fetch: target name too long.\n");
    return 1;
  }
  const char *tag = tag_opt;
  char *at = strchr(spec, '@');
  if (at) { *at = '\0'; tag = at + 1; }

  const yame_asset_reg_t *hits[64];
  const char *only = NULL;

  /* Whatever the browser showed you is what you can type. The tree names a
   * row as <unit>[/<folder>] -- hg38/data, EPIC/KYCG, hg38 -- while the
   * registry names it <source>/<target>, and until now only the latter was
   * accepted. That left the browser unable to tell you the command for the
   * thing you were looking at: the source appears nowhere in the tree, so
   * "hg38 > data" gave no hint that it is spelled methscope/hg38/data.
   *
   * Both spellings work. The browser path is tried first because it is the
   * one a reader can see; the registry spelling stays valid for anything that
   * already uses it, including the two error messages below and the second
   * column of `fetch -l`. The two cannot be confused -- no browser path
   * matches a source name -- and no browser path is claimed by two rows. */
  /* A name is a scope. "hg38" takes the unit and everything under it, the
   * same as ticking that folder in the browser; "hg38/data" takes one row.
   * The registry's own <source>/<target> still resolves, for anything already
   * written against it. */
  size_t n_sel = collect_scope(spec, hits, 64);

  if (!n_sel) {
    char *slash = strchr(spec, '/');
    if (slash) {
      *slash = '\0';
      const yame_asset_reg_t *a = find_asset(spec, slash + 1);
      if (a) { hits[0] = a; n_sel = 1; }
      else *slash = '/';
    }
  }

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

    const yame_asset_reg_t *hit[16];
    char hitpath[16][544];   /* a 264-byte browser path plus a file name */
    size_t n_hit = 0, n_claim = 0;   /* shown, and how many there really are */
    for (size_t i = 0; i < YAME_ASSETS_N; ++i) {
      char u[128], sb[128], full[264];
      unit_of(&YAME_ASSETS[i], u, sizeof(u), sb, sizeof(sb));
      if (sb[0]) snprintf(full, sizeof(full), "%s/%s", u, sb);
      else       snprintf(full, sizeof(full), "%s", u);
      /* A unit index renders at the top of its unit rather than in the
       * directory that publishes it, so what a reader sees is
       * hg38/cpg_nocontig.cr while the file lives in hg38/KYCG. Take the
       * address the browser shows as well as the true one -- the visible
       * spelling should not be the one that fails. Only the three genome
       * .cr files differ this way; an array's ordering is published in the
       * unit it renders under. */
      const char *at = is_unit_index(fname) ? u : full;
      if (head[0] && strcasecmp(full, head) != 0 &&
          strcasecmp(at, head) != 0) continue;
      for (size_t j = 0; j < YAME_ASSETS[i].n_files; ++j)
        if (strcmp(YAME_ASSETS[i].files[j].name, fname) == 0) {
          ++n_claim;
          if (n_hit < 16) {
            hit[n_hit] = &YAME_ASSETS[i];
            snprintf(hitpath[n_hit], sizeof(hitpath[0]), "%.263s/%.270s",
                     at, fname);
            ++n_hit;
          }
          break;
        }
    }
    if (n_claim == 1) {
      hits[0] = hit[0];
      n_sel = 1;
      /* Point into argv, not into `spec`: spec is this function's local and
       * dies on return, while the selection it feeds outlives it. As a
       * local in main_fetch this was accidentally safe; extracting the
       * resolver made it a dangling pointer. */
      only = arg + (fname - spec);
    }
    else if (n_claim > 1) {
      fprintf(stderr, "yame fetch: %zu directories publish a file called "
                      "%s. Name one:\n", n_claim, fname);
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
            "yame fetch: nothing in the catalogue is called %s.\n"
            "  Name it the way the browser shows it: hg38, hg38/KYCG,\n"
            "  hg38/data, EPIC. A name takes everything under it.\n"
            "  `yame fetch -l` lists what there is.\n", spec);
    return 1;
  }
  for (size_t i = 0; i < n_sel; ++i) {
    int dup = 0;
    for (size_t k = 0; k < *n_out; ++k)
      if (out[k].a == hits[i]) {
        if (!only) out[k].only = NULL;   /* the wider name wins */
        dup = 1; break;
      }
    if (dup) continue;
    if (*n_out >= cap) { fprintf(stderr, "yame fetch: too many names.\n"); return 1; }
    out[*n_out].a = hits[i];
    out[*n_out].only = only;
    out[*n_out].tag = tag;
    ++*n_out;
  }
  return 0;
}

int main_fetch(int argc, char *argv[]) {
  const char *dopt = NULL, *tag_override = NULL;
  const char *url = NULL, *sha = NULL, *dest = NULL;
  int force = 0, quiet = 0, unpinned_ok = 0, list = 0, assume_yes = 0;
  int dry_run = 0, here = 0;
  const char *filter = NULL;
  int c;

  while ((c = getopt(argc, argv, "cd:t:kflqu:s:o:yng:h")) >= 0) {
    switch (c) {
    case 'c': here = 1; break;
    case 'd': dopt = optarg; break;
    case 't': tag_override = optarg; break;
    case 'k': unpinned_ok = 1; break;
    case 'f': force = 1; break;
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
  if (list)
    return dump_registry(dopt, optind < argc ? argv[optind] : NULL, filter);

  /* About to actually use the store, so this is where the one-time notice
   * about a pre-consolidation cache belongs -- not in every path that merely
   * asks where the store is. */
  {
    char nroot[4096];
    yame_assets_root(dopt, NULL, nroot, sizeof(nroot));
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
      fprintf(stderr, "yame fetch: -u, -s and -o go together.\n");
      return 1;
    }
    if (yame_assets_download_verify(url, sha, dest, &opt, NULL, &err) != 0) {
      fprintf(stderr, "yame fetch: %s\n", err ? err : "failed");
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
    if (!isatty(STDIN_FILENO) || !isatty(STDOUT_FILENO))
      return dump_registry(dopt, NULL, filter);
    if (browse_catalog(dopt, force) == 0) return 0;
    return dump_registry(dopt, NULL, filter); /* terminal cannot host the widget */
  }

  /* ---- catalogued form: <name>[@tag] ... ---- */
  sel_t sel[64];
  size_t n_sel = 0;
  char shown[512] = "";           /* the names, for the plan and the errors */
  const size_t SELCAP = sizeof(sel)/sizeof(sel[0]);
  for (int ai = optind; ai < argc; ++ai) {
    /* Commas separate names, the same way -m takes CGI,ChromHMM. Split
     * first: parsing @tag before the comma would let `a@v3,b` read the tag
     * as "v3,b" and swallow the second name.
     *
     * The whole string is still tried if any piece fails, so the syntax does
     * not forbid a comma inside a catalogue name -- such a name would fail
     * as pieces and resolve as itself. */
    int rc = 1;
    size_t before = n_sel;
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
        rc = resolve_spec(tok, tag_override, sel, &n_sel, SELCAP, 1);
        if (rc != 0) bad = tok;
      }
      if (rc != 0) n_sel = before;
    }
    if (rc != 0) {
      /* The whole string may be a name in its own right. If it is not, the
       * piece that failed is the more useful thing to name -- "nothing is
       * called nope" beats quoting the entire list back. */
      if (resolve_spec(argv[ai], tag_override, sel, &n_sel, SELCAP,
                       bad != NULL) != 0) {
        if (bad) resolve_spec(bad, tag_override, sel, &n_sel, SELCAP, 0);
        return 1;
      }
    }
    size_t l = strlen(shown);
    snprintf(shown + l, sizeof(shown) - l, "%s%s", l ? " " : "", argv[ai]);
  }


  char root[4096];
  yame_assets_root(dopt, NULL, root, sizeof(root));
  /* -c never touches the store, so a read-only one must not stop it. */
  if (!here && !yame_assets_root_writable(root)) {
    fprintf(stderr,
            "yame fetch: %s is not writable. A read-only shared store is fine "
            "to read from, but nothing can be fetched into it; set -d or "
            "YAME_DATA_HOME to somewhere you own.\n", root);
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
      if (!file_wanted(sel[i].a, sel[i].a->files[j].name, filter, sel[i].only))
        continue;
      ++n_here;
      /* Already-present files are skipped unless -f, so counting their bytes
       * in the total would quote a transfer that is not going to happen. */
      if (!force && sel_present(root, sel[i].a, sel[i].a->files[j].name, here))
        ++n_have;
      else total += sel[i].a->files[j].size;
    }
    if (n_here) { ++n_dirs; n_files += n_here; }
  }

  if (!n_files) {
    if (filter)
      fprintf(stderr, "yame fetch: nothing under %s matches -g %s.\n"
                      "  Terms are ANDed and match the file name, its source, "
                      "collection, title or upstream database.\n", shown, filter);
    else
      fprintf(stderr, "yame fetch: %s holds no files.\n", shown);
    return 1;
  }

  /* Nothing to move is a finished job, not an empty plan: say so plainly and
   * stop, rather than printing a size of nothing and asking to confirm it. */
  if (n_have == n_files) {
    fprintf(stderr, "%s%s%s: all %zu file%s already %s.\n", shown,
            filter ? " -g " : "", filter ? filter : "", n_files,
            n_files == 1 ? "" : "s", here ? "here" : "in the store");
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
    fprintf(stderr, "%s%s%s: %zu file%s in %zu director%s", shown,
            filter ? " -g " : "", filter ? filter : "",
            n_files, n_files == 1 ? "" : "s",
            n_dirs, n_dirs == 1 ? "y" : "ies");
    if (n_have)
      fprintf(stderr, " -- %zu already %s, %zu to fetch",
              n_have, here ? "here" : "in the store", n_files - n_have);
    fprintf(stderr, ", %s\n", hs);

    struct { const char *name; uint64_t size; } *v =
        malloc(n_files * sizeof(*v));
    if (v) {
      size_t k = 0;
      for (size_t i = 0; i < n_sel; ++i)
        for (size_t j = 0; j < sel[i].a->n_files && k < n_files; ++j) {
          if (!file_wanted(sel[i].a, sel[i].a->files[j].name, filter, sel[i].only))
            continue;
          if (!force && sel_present(root, sel[i].a, sel[i].a->files[j].name, here))
            continue;                 /* listing what will move, not what is */
          v[k].name = sel[i].a->files[j].name;
          v[k].size = sel[i].a->files[j].size;
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

  /* Always ask. A fetch writes to a shared store and can be very large, and
   * the size is only knowable from the plan just printed -- so the plan and
   * the question belong together rather than the question being reserved for
   * cases someone guessed would be big. */
  if (!assume_yes) {
    if (!isatty(STDIN_FILENO)) {
      fprintf(stderr,
              "Refusing to fetch %zu file%s (%s) without confirmation. "
              "Re-run with -y.\n", n_files - n_have,
              (n_files - n_have) == 1 ? "" : "s", hs);
      return 1;
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
    const yame_asset_reg_t *a = sel[i].a;
    const char *tag = sel[i].tag ? sel[i].tag : a->tag;
    const char *anchor = a->anchor;

    /* A tag this build does not pin carries no digest for its manifest, so
     * the anchor check has to be dropped -- a decision the caller makes
     * explicitly, not something that happens quietly because they typed a
     * tag. */
    if (sel[i].tag && strcmp(sel[i].tag, a->tag) != 0) {
      if (!unpinned_ok) {
        fprintf(stderr,
                "yame fetch: this build pins %s at %s, so it holds no digest "
                "for %s and cannot verify the manifest published there. "
                "Re-run with -k to accept that, or regenerate the registry "
                "and rebuild.\n", a->target, a->tag, sel[i].tag);
        return 1;
      }
      anchor = NULL;
    }

    char store_sub[4096];
    if (yame_assets_join(store_sub, sizeof(store_sub), root, a->store_sub) != 0) {
      fprintf(stderr, "yame fetch: store path too long.\n");
      return 1;
    }

    int rc = here ? fetch_entry_here(a, tag, filter, sel[i].only, &opt, &err)
                  : fetch_entry(a, root, tag, anchor, filter, sel[i].only,
                                &opt, &err);
    if (rc != 0) {
      fprintf(stderr, "yame fetch: %s\n", err ? err : "failed");
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
