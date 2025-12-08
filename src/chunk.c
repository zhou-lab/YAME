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
  fprintf(stderr, "Usage: yame chunk [options] <in.cx> <outdir>\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -v        verbose\n");
  fprintf(stderr, "    -s        chunk size\n");
  fprintf(stderr, "    -h        This help\n");
  fprintf(stderr, "\n");

  return 1;
}

int main_chunk(int argc, char *argv[]) {

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
  char *outdir;
  if (argc >= optind + 2)
    outdir = strdup(argv[optind+1]);
  else {
    outdir = malloc(strlen(fname)+1000);
    strcpy(outdir, fname);
    strcat(outdir, "_chunks");
  }
  mkdir(outdir, 0777);

  cfile_t cf = open_cfile(fname);
  uint64_t i=0, k;
  for (k=0; ; ++k) {
    cdata_t c = read_cdata1(&cf);
    if (c.n == 0) break;
    
    cdata_t c2 = {0};
    decompress(&c, &c2);
    cdata_t c3 = {0};
    for (i=0; i<=(c2.n/chunk_size); ++i) {
      c3.s = NULL;
      slice(&c2, i*chunk_size, (i+1)*chunk_size-1, &c3);
      cdata_compress(&c3);
      char *tmp = malloc(strlen(outdir) + 1000);
      sprintf(tmp, "%s/%lu.cx", outdir, i);
      if (verbose) fprintf(stdout, "%s\n", tmp);
      if (k) cdata_write(tmp, &c3, "a", verbose);
      else cdata_write(tmp, &c3, "w", verbose);
      free(c3.s);
      free(tmp);
    }
    free(c2.s); free(c.s);
  }
  free(outdir);
  bgzf_close(cf.fh);
  
  return 0;
}
