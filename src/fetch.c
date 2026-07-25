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
 *   yame fetch <source>/<target>[@tag]   a directory this build pins
 *   yame fetch -u URL -s SHA -o DEST     one file, digest given by hand
 *
 * The second is the low-level form: a per-file digest supplied by the caller,
 * which is how a store gets a file this build's registry has never heard of.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "assets.h"
#include "registry.h"
#include "yame_ui.h"

static int usage(void) {
  char root[4096];
  yame_assets_root(NULL, NULL, root, sizeof(root));

  fprintf(stderr, "\n");
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "  yame fetch                              browse the catalogue\n");
  fprintf(stderr, "  yame fetch [options] <source>/<target>[@tag]\n");
  fprintf(stderr, "  yame fetch [options] -u <url> -s <sha256> -o <dest>\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Browsing:\n");
  fprintf(stderr, "  With no target on a terminal, opens a tree browser: arrows move,\n");
  fprintf(stderr, "  right/left open and close a source, space checks a row, f fetches\n");
  fprintf(stderr, "  what is checked, q leaves. What is already in the store shows as\n");
  fprintf(stderr, "  present and cannot be checked.\n");
  fprintf(stderr, "  Piped or redirected it prints the plain listing instead, so a script\n");
  fprintf(stderr, "  never blocks on a keystroke.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Purpose:\n");
  fprintf(stderr, "  Download reference assets into the shared store that every tool in the\n");
  fprintf(stderr, "  suite reads, verifying each file against a digest this build pins.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Store:\n");
  fprintf(stderr, "  Resolved in order: -d, $YAME_DATA_HOME,\n");
  fprintf(stderr, "  ${XDG_DATA_HOME:-~/.local/share}/yame\n");
  fprintf(stderr, "  Currently: %s\n", root);
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "  -d <dir>      Store root, overriding the environment.\n");
  fprintf(stderr, "  -t <tag>      Upstream tag, overriding the pinned one. Without a digest\n");
  fprintf(stderr, "                for that tag nothing can be verified, so this needs -k.\n");
  fprintf(stderr, "  -k            Accept an unpinned tag (no anchor check). Files are still\n");
  fprintf(stderr, "                verified against the manifest that tag publishes.\n");
  fprintf(stderr, "  -f            Re-download what is present, and overwrite a store that\n");
  fprintf(stderr, "                was populated from a different tag.\n");
  fprintf(stderr, "  -l            List what this build knows how to fetch, and exit.\n");
  fprintf(stderr, "  -q            No progress output.\n");
  fprintf(stderr, "  -u <url>      Single-file form: what to download.\n");
  fprintf(stderr, "  -s <sha256>   Single-file form: the digest it must have.\n");
  fprintf(stderr, "  -o <dest>     Single-file form: where it goes (a path, not a dir).\n");
  fprintf(stderr, "  -h            This help.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Notes:\n");
  fprintf(stderr, "  * A directory records which upstream tag filled it, in the SHA256SUMS\n");
  fprintf(stderr, "    it keeps. A build pinned elsewhere refuses to overwrite it rather\n");
  fprintf(stderr, "    than start a re-download war; -f overrules that.\n");
  fprintf(stderr, "  * `shasum -a 256 -c SHA256SUMS` in any store directory re-verifies it\n");
  fprintf(stderr, "    by hand, with none of this code involved.\n");
  if (!yame_assets_have_curl())
    fprintf(stderr, "  * THIS BUILD HAS NO LIBCURL: fetching is unavailable.\n");
  fprintf(stderr, "\n");
  return 1;
}

static int list_catalog(void) {
  char root[4096];
  yame_assets_root(NULL, NULL, root, sizeof(root));

  printf("%-20s %-14s %-8s %s\n", "SOURCE", "TARGET", "TAG", "STORE");
  for (size_t i = 0; i < YAME_ASSETS_N; ++i) {
    const yame_asset_reg_t *a = &YAME_ASSETS[i];
    char dir[4096];
    yame_assets_join(dir, sizeof(dir), root, a->store_sub);

    const char *state = "-";
    switch (yame_assets_pin_check(dir, a->anchor)) {
    case YAME_PIN_MATCH:    state = "present";  break;
    case YAME_PIN_CONFLICT: state = "OTHER TAG"; break;
    case YAME_PIN_UNKNOWN:  state = "unpinned"; break;
    default:                state = "-";        break;
    }
    printf("%-20s %-14s %-8s %s [%s]\n",
           a->source, a->target, a->tag, a->store_sub, state);
  }
  return 0;
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
 * The tree widget from src/ui.c. Sources are the top level and their targets
 * unfold underneath; what is already in the store is shown as present and
 * cannot be checked, since there is nothing to ask for.
 */

typedef struct {
  char   root[4096];            /* the store */
  char   sel[64][256];          /* "source/target", from the checked rows */
  size_t n_sel;
  int    force;
} browse_t;

/* Is this entry's directory already filled, and by our tag? */
static int asset_state(const char *store_root, const yame_asset_reg_t *a) {
  char dir[4096];
  if (yame_assets_join(dir, sizeof(dir), store_root, a->store_sub) != 0)
    return YAME_PIN_ABSENT;
  return yame_assets_pin_check(dir, a->anchor);
}

static void bx_expand(void *ctx, const char *row, yame_ui_kids_t *out) {
  browse_t *b = ctx;
  size_t cap = 0;

  for (size_t i = 0; i < YAME_ASSETS_N; ++i)
    if (strcmp(YAME_ASSETS[i].source, row) == 0) ++cap;
  if (!cap) return;

  out->rows   = calloc(cap, sizeof(char *));
  out->keys   = calloc(cap, sizeof(char *));
  out->styles = calloc(cap, 1);
  if (!out->rows || !out->keys || !out->styles) return;

  for (size_t i = 0; i < YAME_ASSETS_N; ++i) {
    const yame_asset_reg_t *a = &YAME_ASSETS[i];
    if (strcmp(a->source, row) != 0) continue;

    int st = asset_state(b->root, a);
    const char *note = st == YAME_PIN_MATCH    ? "in the store"
                     : st == YAME_PIN_CONFLICT ? "ANOTHER TAG"
                     : st == YAME_PIN_UNKNOWN  ? "unverified"
                                               : "not fetched";
    char line[256], key[256];
    snprintf(line, sizeof(line), "%-16s %-6s %s", a->target, a->tag, note);
    snprintf(key, sizeof(key), "%s/%s", a->source, a->target);

    out->rows[out->n] = strdup(line);
    out->keys[out->n] = strdup(key);
    /* Present rows are not checkable; a conflicting one IS, since re-fetching
     * (with -f) is exactly how you resolve it. */
    out->styles[out->n] = (unsigned char)(st == YAME_PIN_MATCH ? YAME_ROW_HAVE
                                                               : YAME_ROW_MISSING);
    ++out->n;
  }
}

static void bx_accept(void *ctx, const char *root, const char *key) {
  browse_t *b = ctx;
  (void)root;
  if (!key || b->n_sel >= sizeof(b->sel)/sizeof(b->sel[0])) return;
  snprintf(b->sel[b->n_sel], sizeof(b->sel[0]), "%s", key);
  ++b->n_sel;
}

/* Fetch one catalogued entry into the store. Shared by the browser and the
 * one-shot command line so they cannot drift apart. */
static int fetch_entry(const yame_asset_reg_t *a, const char *store_root,
                       const char *tag, const char *anchor,
                       const yame_fetch_opt_t *opt, char **err) {
  char store_sub[4096];
  if (yame_assets_join(store_sub, sizeof(store_sub), store_root, a->store_sub) != 0)
    return -1;
  return yame_assets_fetch_subtree(a->base_url, tag, a->remote_sub, store_sub,
                                   anchor, opt, err);
}

/* Runs inside the widget: the tree suspends, this draws in a panel, and the
 * tree resumes with the rows reloaded so the new state is visible. */
static void bx_commit(void *ctx) {
  browse_t *b = ctx;
  int ok = 0, bad = 0;

  if (!b->n_sel) return;

  yame_ui_panel_open(4);
  for (size_t i = 0; i < b->n_sel; ++i) {
    char spec[256];
    snprintf(spec, sizeof(spec), "%s", b->sel[i]);
    char *slash = strchr(spec, '/');
    /* Two slashes are possible ("InfiniumAnnotation/MSA/KYCG"): the source is
     * everything before the FIRST one, the target everything after. */
    if (!slash) { ++bad; continue; }
    *slash = '\0';

    const yame_asset_reg_t *a = find_asset(spec, slash + 1);
    if (!a) { ++bad; continue; }

    yame_ui_panel_line(0, "fetching %s/%s  (%zu of %zu)",
                       a->source, a->target, i + 1, b->n_sel);

    yame_fetch_opt_t opt = {0};
    opt.force = b->force;
    opt.quiet = 1;                 /* the panel is the only output surface */

    char *err = NULL;
    if (fetch_entry(a, b->root, a->tag, a->anchor, &opt, &err) == 0) {
      ++ok;
      yame_ui_panel_line(1, "ok");
    } else {
      ++bad;
      yame_ui_panel_line(1, "failed: %s", err ? err : "(no detail)");
    }
    free(err);
  }

  yame_ui_panel_line(0, "%d fetched, %d failed", ok, bad);
  yame_ui_panel_line(1, "");
  yame_ui_panel_pause(2, "press any key to return to the catalogue");
  yame_ui_panel_close();
  b->n_sel = 0;
}

/* The catalogue as a browsable tree. Returns 0 when it ran, -1 when the
 * terminal cannot host it and the caller should fall back to the flat list. */
static int browse_catalog(const char *dopt, int force) {
  static char *roots[8];
  static unsigned char styles[8];
  size_t n_roots = 0;
  browse_t b;

  memset(&b, 0, sizeof(b));
  b.force = force;
  yame_assets_root(dopt, NULL, b.root, sizeof(b.root));

  /* One root per source, in catalogue order, deduplicated. */
  for (size_t i = 0; i < YAME_ASSETS_N && n_roots < 8; ++i) {
    int seen = 0;
    for (size_t k = 0; k < n_roots; ++k)
      if (strcmp(roots[k], YAME_ASSETS[i].source) == 0) { seen = 1; break; }
    if (!seen) {
      roots[n_roots] = strdup(YAME_ASSETS[i].source);
      styles[n_roots] = YAME_ROW_PLAIN;
      ++n_roots;
    }
  }

  char title[4200];
  snprintf(title, sizeof(title), "yame fetch  --  %s", b.root);

  yame_ui_tree_t spec;
  memset(&spec, 0, sizeof(spec));
  spec.title       = title;
  spec.header      = NULL;
  spec.roots       = roots;
  spec.root_styles = styles;
  spec.n_roots     = n_roots;
  spec.expand      = bx_expand;
  spec.actions[0].key    = 'f';
  spec.actions[0].verb   = "fetch";
  spec.actions[0].accept = bx_accept;
  spec.actions[0].commit = bx_commit;   /* non-NULL: the tree stays open */
  spec.n_actions   = 1;
  spec.ctx         = &b;

  return yame_ui_tree(&spec) < 0 ? -1 : 0;
}

/* Progress: one line per file, rewritten in place while bytes move. Nothing
 * fancy -- a tool with its own renderer passes its own hooks to the library
 * instead of inheriting this one. */
static void prog_begin(void *ud, const char *name, uint64_t total) {
  (void)ud; (void)total;
  /* Only on a terminal: the '\r' that makes this an in-place indicator turns
   * into a duplicated line in a log file or a CI transcript. */
  if (!isatty(STDERR_FILENO)) return;
  fprintf(stderr, "  %-48s ...\r", name);
  fflush(stderr);
}

static void prog_done(void *ud, const char *name, uint64_t bytes, int ok) {
  (void)ud;
  if (ok) {
    double mb = (double)bytes / (1024.0 * 1024.0);
    fprintf(stderr, "  %-48s %8.1f MB\n", name, mb);
  } else {
    fprintf(stderr, "  %-48s failed\n", name);
  }
}

int main_fetch(int argc, char *argv[]) {
  const char *dopt = NULL, *tag_override = NULL;
  const char *url = NULL, *sha = NULL, *dest = NULL;
  int force = 0, quiet = 0, unpinned_ok = 0;
  int c;

  while ((c = getopt(argc, argv, "d:t:kflqu:s:o:h")) >= 0) {
    switch (c) {
    case 'd': dopt = optarg; break;
    case 't': tag_override = optarg; break;
    case 'k': unpinned_ok = 1; break;
    case 'f': force = 1; break;
    case 'l': return list_catalog();
    case 'q': quiet = 1; break;
    case 'u': url = optarg; break;
    case 's': sha = optarg; break;
    case 'o': dest = optarg; break;
    case 'h': return usage();
    default: return usage();
    }
  }

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
    if (!isatty(STDIN_FILENO) || !isatty(STDOUT_FILENO)) return list_catalog();
    if (browse_catalog(dopt, force) == 0) return 0;
    return list_catalog();            /* terminal cannot host the widget */
  }

  /* ---- catalogued form: <source>/<target>[@tag] ---- */
  char spec[512];
  if (snprintf(spec, sizeof(spec), "%s", argv[optind]) >= (int)sizeof(spec)) {
    fprintf(stderr, "yame fetch: target name too long.\n");
    return 1;
  }

  char *at = strchr(spec, '@');
  if (at) { *at = '\0'; tag_override = at + 1; }

  char *slash = strchr(spec, '/');
  if (!slash) {
    fprintf(stderr,
            "yame fetch: expected <source>/<target>, e.g. "
            "InfiniumAnnotation/MSA. Run `yame fetch -l` for the list.\n");
    return 1;
  }
  *slash = '\0';

  const yame_asset_reg_t *a = find_asset(spec, slash + 1);
  if (!a) {
    fprintf(stderr,
            "yame fetch: this build knows nothing about %s/%s. "
            "Run `yame fetch -l` for the list.\n", spec, slash + 1);
    return 1;
  }

  const char *tag = tag_override ? tag_override : a->tag;
  const char *anchor = a->anchor;

  /* A tag this build does not pin carries no digest for its manifest, so the
   * anchor check has to be dropped -- which is a decision the caller makes
   * explicitly, not something that happens quietly because they typed a tag. */
  if (tag_override && strcmp(tag_override, a->tag) != 0) {
    if (!unpinned_ok) {
      fprintf(stderr,
              "yame fetch: this build pins %s at %s, so it holds no digest for "
              "%s and cannot verify the manifest published there. Re-run with "
              "-k to accept that, or regenerate the registry and rebuild.\n",
              a->target, a->tag, tag_override);
      return 1;
    }
    anchor = NULL;
  }

  char root[4096], store_sub[4096];
  yame_assets_root(dopt, NULL, root, sizeof(root));
  if (yame_assets_join(store_sub, sizeof(store_sub), root, a->store_sub) != 0) {
    fprintf(stderr, "yame fetch: store path too long.\n");
    return 1;
  }

  if (!yame_assets_root_writable(root)) {
    fprintf(stderr,
            "yame fetch: %s is not writable. A read-only shared store is fine "
            "to read from, but nothing can be fetched into it; set -d or "
            "YAME_DATA_HOME to somewhere you own.\n", root);
    return 1;
  }

  if (!quiet)
    fprintf(stderr, "%s/%s @ %s -> %s\n", a->source, a->target, tag, store_sub);

  if (fetch_entry(a, root, tag, anchor, &opt, &err) != 0) {
    fprintf(stderr, "yame fetch: %s\n", err ? err : "failed");
    free(err);
    return 1;
  }
  free(err);

  if (!quiet) fprintf(stderr, "done.\n");
  return 0;
}
