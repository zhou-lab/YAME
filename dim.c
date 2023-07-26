#include "cgfile.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: yame dim [options] <in.cg>\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -h        This help\n");
  fprintf(stderr, "\n");

  return 1;
}

int main_dim(int argc, char *argv[]) {

  int c;
  while ((c = getopt(argc, argv, "vh"))>=0) {
    switch (c) {
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 1 > argc) { 
    usage(); 
    wzfatal("Please supply input file.\n"); 
  }

  cgfile_t cgf = open_cgfile(argv[optind]);
  int i = 0;
  for (i=0; ; ++i) {
    cgdata_t cg = read_cg(&cgf);
    if (cg.n == 0) break;
    cgdata_t expanded = {0};
    decompress(&cg, &expanded);
    fprintf(stdout, "%d\t%"PRIu64"\n", i+1, expanded.n);
    free(expanded.s); free(cg.s);
  }
  bgzf_close(cgf.fh);
  
  return 0;
}
