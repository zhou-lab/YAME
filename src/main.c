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

#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "cfile.h"
#include "yame_version.h"

const int unit_base[40] = {
  0,  1,  1,  4,  4,  8,  8,  0,
  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  0,  0,  0,  8,
  8,  0,  0,  0,  0,  0,  0,  0
};

int main_pack(int argc, char *argv[]);
int main_unpack(int argc, char *argv[]);
int main_hprint(int argc, char *argv[]);
int main_subset(int argc, char *argv[]);
int main_rowsub(int argc, char *argv[]);
int main_split(int argc, char *argv[]);
int main_pairwise(int argc, char *argv[]);
int main_info(int argc, char *argv[]);
int main_summary(int argc, char *argv[]);
int main_chunk(int argc, char *argv[]);
int main_chunkchar(int argc, char *argv[]);
int main_rowop(int argc, char *argv[]);
int main_index(int argc, char *argv[]);
int main_mask(int argc, char *argv[]);
int main_dsample(int argc, char *argv[]);
int main_binarize(int argc, char *argv[]);
int main_perturb(int argc, char *argv[]);
int main_fetch(int argc, char *argv[]);

#include "yame_ui.h"
#include "assets.h"

#define PACKAGE_VERSION YAME_VERSION

/* Help is rendered with colour when stderr is a terminal; yame_ui_* return
 * empty strings otherwise, so a redirected --help stays plain text. */
static void sec(const char *title) {
  fprintf(stderr, "%s%s%s\n", yame_ui_bold(), title, yame_ui_reset());
}

static void cmd(const char *name, const char *desc) {
  fprintf(stderr, "  %s%-12s%s %s\n",
          yame_ui_green(), name, yame_ui_reset(), desc);
}

/* A continuation line for a command whose description needs two. */
static void cmd2(const char *desc) {
  fprintf(stderr, "               %s%s%s\n",
          yame_ui_dim(), desc, yame_ui_reset());
}

static int usage(void)
{
  char store[4096];
  yame_assets_root(NULL, NULL, store, sizeof(store));

  fprintf(stderr, "\n");
  fprintf(stderr, "%syame%s (Yet Another Methylation Encoder)  %s%s%s\n",
          yame_ui_bold(), yame_ui_reset(),
          yame_ui_dim(), PACKAGE_VERSION, yame_ui_reset());
  fprintf(stderr, "%sWhole-genome DNA methylation data management using CX "
                  "formats.%s\n", yame_ui_dim(), yame_ui_reset());
  fprintf(stderr, "%sContact: Wanding Zhou <wanding.zhou@pennmedicine.upenn.edu>"
                  "%s\n", yame_ui_dim(), yame_ui_reset());
  fprintf(stderr, "\n");

  /* Where the reference data is, named by the variable that moves it -- not a
   * command, so not in the list of them, and worth seeing before them. */
  {
    const char *env = getenv("YAME_DATA_HOME");
    fprintf(stderr, "%sYAME_DATA_HOME: %s%s%s\n", yame_ui_blue(), store,
            (env && *env) ? "" : " (unset, default)", yame_ui_reset());
  }
  fprintf(stderr, "\n");

  sec("Usage:");
  fprintf(stderr, "  yame <command> [options] [args]\n");
  fprintf(stderr, "\n");

  sec("Core I/O:");
  cmd("pack",   "Pack text/bed-like inputs into a .cx stream");
  cmd("unpack", "Unpack a .cx stream back to text");
  cmd("hprint", "Horizontal printing (primarily format 6)");
  fprintf(stderr, "\n");

  sec("Indexing / file management:");
  cmd("index", "Create/refresh a sample index for a .cx file");
  cmd("split", "Split a multi-sample .cx into single-sample files");
  cmd("info",  "Show basic metadata/parameters of a .cx file");
  fprintf(stderr, "\n");

  sec("Subsetting / chunking:");
  cmd("subset",    "Subset samples from a .cx (or terms from format 2 with -s)");
  cmd("rowsub",    "Subset rows by index list / mask / coordinates / block range");
  cmd("chunk",     "Chunk binary CX into smaller fragments");
  cmd("chunkchar", "Chunk text data into smaller fragments");
  fprintf(stderr, "\n");

  sec("Summaries / comparisons:");
  cmd("summary",  "Summarize query features, optionally against masks");
  cmd("pairwise", "Call pairwise differential methylation (fmt3 -> fmt6)");
  fprintf(stderr, "\n");

  sec("Transforms / utilities:");
  cmd("binarize", "Convert fmt3 (M/U) to fmt6 (set+universe) by beta/M threshold");
  cmd("mask",     "Invalidate specific CpG sites using an external mask file");
  cmd2("(deterministic, file-driven: controls *which* sites are valid)");
  cmd("dsample",  "Randomly subsample N covered CpG sites, masking the rest");
  cmd2("(stochastic, rate-driven: reduces site count and coverage)");
  cmd("perturb",  "Randomly flip 0/1 bits in fmt0 or fmt6 (noise injection)");
  cmd("rowop",    "Row-wise operations (e.g., sum / combine binary tracks)");
  fprintf(stderr, "\n");

  sec("Reference data:");
  cmd("fetch", "Download reference assets into the shared store");
  cmd2("('yame fetch' to browse, 'yame fetch -l' to list)");
  fprintf(stderr, "\n");

  fprintf(stderr, "Run '%syame <command> -h%s' for command-specific options.\n",
          yame_ui_bold(), yame_ui_reset());
  fprintf(stderr, "\n");
  return 1;
}

int main(int argc, char *argv[]) {
  int ret;
  if (argc < 2) return usage();
  if (strcmp(argv[1], "pack") == 0) ret = main_pack(argc-1, argv+1);
  else if (strcmp(argv[1], "unpack") == 0) ret = main_unpack(argc-1, argv+1);
  else if (strcmp(argv[1], "hprint") == 0) ret = main_hprint(argc-1, argv+1);
  else if (strcmp(argv[1], "subset") == 0) ret = main_subset(argc-1, argv+1);
  else if (strcmp(argv[1], "rowsub") == 0) ret = main_rowsub(argc-1, argv+1);
  else if (strcmp(argv[1], "split") == 0) ret = main_split(argc-1, argv+1);
  else if (strcmp(argv[1], "pairwise") == 0) ret = main_pairwise(argc-1, argv+1);
  else if (strcmp(argv[1], "info") == 0) ret = main_info(argc-1, argv+1);
  else if (strcmp(argv[1], "summary") == 0) ret = main_summary(argc-1, argv+1);
  else if (strcmp(argv[1], "index") == 0) ret = main_index(argc-1, argv+1);
  else if (strcmp(argv[1], "chunk") == 0) ret = main_chunk(argc-1, argv+1);
  else if (strcmp(argv[1], "chunkchar") == 0) ret = main_chunkchar(argc-1, argv+1);
  else if (strcmp(argv[1], "rowop") == 0) ret = main_rowop(argc-1, argv+1);
  else if (strcmp(argv[1], "mask") == 0) ret = main_mask(argc-1, argv+1);
  else if (strcmp(argv[1], "binarize") == 0) ret = main_binarize(argc-1, argv+1);
  else if (strcmp(argv[1], "dsample") == 0) ret = main_dsample(argc-1, argv+1);
  else if (strcmp(argv[1], "perturb") == 0) ret = main_perturb(argc-1, argv+1);
  else if (strcmp(argv[1], "fetch") == 0) ret = main_fetch(argc-1, argv+1);
  else {
    fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
    return 1;
  }

  fflush(stdout);             /* not enough for remote file systems */
  fclose(stdout);

  return ret;
}
