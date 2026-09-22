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
  yame_usage_head("yame chunkchar [options] <in.txt>");
  yame_usage_sec("Options:");
  yame_usage_opt("-v", "verbose");
  yame_usage_opt("-s", "chunk size");
  yame_usage_opt("-h", "This help");
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
  char *outdir = wzmalloc(strlen(fname)+1000);
  strcpy(outdir, fname);
  strcat(outdir, "_chunks");
  mkdir(outdir, 0777);

  gzFile fh = wzopen(fname, 1);
  char **lines = NULL; uint64_t n = 0;
  char *line = NULL;
  while (gzFile_read_line(fh, &line) > 0) {
    lines = wzrealloc(lines, (n+1)*sizeof(char*));
    lines[n++] = wzstrdup(line);
  }
  free(line);
  gzclose(fh);

  uint64_t u,i;
  for (u=0; u<=n/chunk_size; ++u) {
    char *tmp = wzmalloc(strlen(outdir) + 1000);
    sprintf(tmp, "%s/%"PRIu64".txt", outdir, u);
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
