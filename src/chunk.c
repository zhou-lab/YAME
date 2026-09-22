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

#include <sys/stat.h>
#include "yame_ui.h"
#include <sys/types.h>
#include "cfile.h"

static int usage() {
  yame_usage_head("yame chunk [options] <in.cx> <outdir>");
  yame_usage_sec("Options:");
  yame_usage_opt("-v", "verbose");
  yame_usage_opt("-s", "chunk size");
  yame_usage_opt("-h", "This help");
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
    outdir = wzstrdup(argv[optind+1]);
  else {
    outdir = wzmalloc(strlen(fname)+1000);
    strcpy(outdir, fname);
    strcat(outdir, "_chunks");
  }
  mkdir(outdir, 0777);

  cfile_t cf = open_cfile(fname);
  uint64_t i=0, k;
  for (k=0; ; ++k) {
    cdata_t c = read_cdata1(&cf);
    if (c.n == 0) break;
    
    cdata_t c2 = decompress(c);
    cdata_t c3 = {0};
    /* Ceiling, not floor-plus-one: when n is an exact multiple of the chunk
     * size, the extra iteration began at row n, slice() clamped end below
     * beg, and the run aborted on "Slicing negative span". */
    uint64_t n_chunks = (c2.n + chunk_size - 1) / chunk_size;
    for (i=0; i<n_chunks; ++i) {
      c3.s = NULL;
      slice(&c2, i*chunk_size, (i+1)*chunk_size-1, &c3);
      cdata_compress(&c3);
      char *tmp = wzmalloc(strlen(outdir) + 1000);
      sprintf(tmp, "%s/%"PRIu64".cx", outdir, i);
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
