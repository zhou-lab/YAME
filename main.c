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

#define PACKAGE_VERSION "v1.1.0"

static int usage()
{
  fprintf(stderr, "\n");
  fprintf(stderr, "Program: yame (Yet Another Methylation Encoder) - whole genome DNA methylation data management.\n");
  fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
  fprintf(stderr, "Contact: Wanding Zhou<wanding.zhou@pennmedicine.upenn.edu>\n\n");
  fprintf(stderr, "Usage:   yame <command> [options]\n\n");
  fprintf(stderr, "Available commands:\n");
  fprintf(stderr, "     pack         - Pack data into a cx file.\n");
  fprintf(stderr, "     unpack       - Unpack data from a cx file.\n");
  fprintf(stderr, "     hprint       - Print data horizontally (format 6).\n");
  fprintf(stderr, "     subset       - Subset samples or terms from a cx file.\n");
  fprintf(stderr, "     rowsub       - Subset rows a cx file using an index list file.\n");
  fprintf(stderr, "     info         - Display basic parameter of the cx file.\n");
  fprintf(stderr, "     summary      - calculate summary, with or without masks.\n");
  fprintf(stderr, "     index        - Index samples in a cx file.\n");
  fprintf(stderr, "     split        - Split multi-sample data into single-sample data.\n");
  fprintf(stderr, "     pairwise     - Pairwise differential methylation calling.\n");
  fprintf(stderr, "     chunk        - Chunk data into smaller fragments.\n");
  fprintf(stderr, "     chunkchar    - Chunk text data into smaller fragments.\n");
  fprintf(stderr, "     rowop        - Perform operations on rows, e.g., sum binary values.\n");
  fprintf(stderr, "     mask         - Mask methylation data by setting masked record to M=U=0.\n");
  fprintf(stderr, "     binarize     - Binarize methylation data by methylation level (beta).\n");
  fprintf(stderr, "     dsample      - Downsample methylation data (works with fmt 3 or 6).\n");
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
