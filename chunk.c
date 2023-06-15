#include <sys/stat.h>
#include <sys/types.h>
#include "kycg.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: kycg chunk [options] <in.cg> <outdir>\n");
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

  cgfile_t cgf = open_cgfile(fname);
  uint64_t i=0, k;
  for (k=0; ; ++k) {
    cgdata_t cg = read_cg(&cgf);
    if (cg.n == 0) break;
    
    cgdata_t cg2 = {0};
    decompress(&cg, &cg2);
    cgdata_t cg3 = {0};
    for (i=0; i<=(cg2.n/chunk_size); ++i) {
      cg3.s = NULL;
      slice(&cg2, i*chunk_size, (i+1)*chunk_size-1, &cg3);
      recompress(&cg3);
      char *tmp = malloc(strlen(outdir) + 1000);
      sprintf(tmp, "%s/%lu.cg", outdir, i);
      if (verbose) fprintf(stdout, "%s\n", tmp);
      if (k) cgdata_write(tmp, &cg3, "ab", verbose);
      else cgdata_write(tmp, &cg3, "wb", verbose);
      free(cg3.s);
      free(tmp);
    }
    free(cg2.s); free(cg.s);
  }
  free(outdir);
  bgzf_close(cgf.fh);
  
  return 0;
}
