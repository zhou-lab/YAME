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
 * A thin driver over src/assets.c, deliberately not a browser: kycg's `fetch`
 * has an interactive catalogue with per-knowledgebase descriptions and that
 * stays where it is. This exists so the store can be filled without any
 * particular tool installed, and so the engine every tool links has a way to
 * be exercised on its own.
 *
 * Two forms:
 *   yame fetch <source>/<target>[@tag]   a directory this build pins
 *   yame fetch -u URL -s SHA -o DEST     one file, digest given by hand
 *
 * The second is the low-level form. It is what methscope-cli's per-file
 * compiled digests will map onto, and it is how a store gets a file this
 * build's registry has never heard of.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "assets.h"
#include "registry.h"

static int usage(void) {
  char root[4096];
  yame_assets_root(NULL, NULL, root, sizeof(root));

  fprintf(stderr, "\n");
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "  yame fetch [options] <source>/<target>[@tag]\n");
  fprintf(stderr, "  yame fetch [options] -u <url> -s <sha256> -o <dest>\n");
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

  if (optind >= argc) return usage();

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

  if (yame_assets_fetch_subtree(a->base_url, tag, a->remote_sub, store_sub,
                                anchor, &opt, &err) != 0) {
    fprintf(stderr, "yame fetch: %s\n", err ? err : "failed");
    free(err);
    return 1;
  }
  free(err);

  if (!quiet) fprintf(stderr, "done.\n");
  return 0;
}
