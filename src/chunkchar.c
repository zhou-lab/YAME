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

#include <sys/stat.h>
#include <sys/types.h>
#include "cfile.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: yame chunkchar [options] <in.txt>\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -v        verbose\n");
  fprintf(stderr, "    -s        chunk size\n");
  fprintf(stderr, "    -h        This help\n");
  fprintf(stderr, "\n");

  return 1;
}

int main_chunkchar(int argc, char *argv[]) {

  int c, verbose = 0;
  uint64_t chunk_size = 1000000;
  while ((c = getopt(argc, argv, "s:vh"))>=0) {
    switch (c) {
    case 'v': verbose = 1; break;
    case 's': chunk_size = atoi(optarg); break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 1 > argc) {
    usage(); 
    wzfatal("Please supply input file.\n");
  }

  char *fname = argv[optind];
  char *outdir = malloc(strlen(fname)+1000);
  strcpy(outdir, fname);
  strcat(outdir, "_chunks");
  mkdir(outdir, 0777);

  gzFile fh = wzopen(fname, 1);
  char **lines = NULL; uint64_t n = 0;
  char *line = NULL;
  while (gzFile_read_line(fh, &line) > 0) {
    lines = realloc(lines, (n+1)*sizeof(char*));
    lines[n++] = strdup(line);
  }
  free(line);
  gzclose(fh);

  uint64_t u,i;
  for (u=0; u<=n/chunk_size; ++u) {
    char *tmp = malloc(strlen(outdir) + 1000);
    sprintf(tmp, "%s/%lu.txt", outdir, u);
    if (verbose) fprintf(stdout, "%s\n", tmp);
    FILE *fh = fopen(tmp, "w");
    for (i=u*chunk_size; i<(u+1)*chunk_size; ++i) {
      if (i>=n) break;
      fprintf(fh, "%s\n", lines[i]);
      free(lines[i]);
    }
    fclose(fh);
    free(tmp);
  }
  free(lines);
  free(outdir);
  
  return 0;
}
