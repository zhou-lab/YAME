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

#define PACKAGE_VERSION "v1.8"

static int usage(void)
{
  fprintf(stderr, "\n");
  fprintf(stderr, "yame (Yet Another Methylation Encoder)\n");
  fprintf(stderr, "Whole-genome DNA methylation data management using CX formats.\n");
  fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
  fprintf(stderr, "Contact: Wanding Zhou <wanding.zhou@pennmedicine.upenn.edu>\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "  yame <command> [options] [args]\n");
  fprintf(stderr, "\n");

  fprintf(stderr, "Core I/O:\n");
  fprintf(stderr, "  pack         Pack text/bed-like inputs into a .cx stream\n");
  fprintf(stderr, "  unpack       Unpack a .cx stream back to text\n");
  fprintf(stderr, "  hprint       Horizontal printing (primarily format 6)\n");
  fprintf(stderr, "\n");

  fprintf(stderr, "Indexing / file management:\n");
  fprintf(stderr, "  index        Create/refresh a sample index for a .cx file\n");
  fprintf(stderr, "  split        Split a multi-sample .cx into single-sample files\n");
  fprintf(stderr, "  info         Show basic metadata/parameters of a .cx file\n");
  fprintf(stderr, "\n");

  fprintf(stderr, "Subsetting / chunking:\n");
  fprintf(stderr, "  subset       Subset samples from a .cx (or terms from format 2 with -s)\n");
  fprintf(stderr, "  rowsub       Subset rows by index list / mask / coordinates / block range\n");
  fprintf(stderr, "  chunk        Chunk binary CX into smaller fragments\n");
  fprintf(stderr, "  chunkchar    Chunk text data into smaller fragments\n");
  fprintf(stderr, "\n");

  fprintf(stderr, "Summaries / comparisons:\n");
  fprintf(stderr, "  summary      Summarize query features, optionally against masks\n");
  fprintf(stderr, "  pairwise     Call pairwise differential methylation (fmt3 -> fmt6)\n");
  fprintf(stderr, "\n");

  fprintf(stderr, "Transforms / utilities:\n");
  fprintf(stderr, "  binarize     Convert fmt3 (M/U) to fmt6 (set+universe) by beta/M threshold\n");
  fprintf(stderr, "  mask         Mask methylation data (e.g., set M=U=0 for masked sites)\n");
  fprintf(stderr, "  dsample      Downsample methylation data (fmt3 or fmt6)\n");
  fprintf(stderr, "  rowop        Row-wise operations (e.g., sum / combine binary tracks)\n");
  fprintf(stderr, "\n");

  fprintf(stderr, "Run 'yame <command> -h' for command-specific options and details.\n");
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
  else {
    fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
    return 1;
  }

  fflush(stdout);             /* not enough for remote file systems */
  fclose(stdout);

  return ret;
}
