#include <sys/stat.h>
#include "kycg.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: kycg chunk [options] <in.cg> out_prefix\n");
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
  while ((c = getopt(argc, argv, "n:vh"))>=0) {
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
  strcat(outdir, "_chunks");
  mkdir(outdir, 777);

  cgfile_t cgf = open_cgfile(fname);
  cgdata_t cg = read_cg(&cgf);
  cgdata_t cg2 = {0};
  decompress(&cg, &cg2);
  uint64_t i=0;
  cgdata_t cg3 = {0};
  for (i=0; i<=(cg2.n/chunk_size); ++i) {
    slice(&cg2, i*chunk_size, (i+1)*chunk_size-1, &cg3);
    recompress(&cg3);
    char *tmp = malloc(strlen(outdir) + 1000);
    sprintf(tmp, "%s/%lu", outdir, i);
    cgdata_write(tmp, &cg3, verbose);
    free(cg3.s); free(tmp);
  }
  
  free(outdir);
  
  return 0;
}
