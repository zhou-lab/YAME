#include "cfile.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: yame dim [options] <in.cx>\n");
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

  cfile_t cf = open_cfile(argv[optind]);
  int i = 0;
  for (i=0; ; ++i) {
    cdata_t c = read_cdata1(&cf);
    if (c.n == 0) break;
    cdata_t expanded = {0};
    decompress(&c, &expanded);
    fprintf(stdout, "%d\t%"PRIu64"\n", i+1, expanded.n);
    free(expanded.s); free(c.s);
  }
  bgzf_close(cf.fh);
  
  return 0;
}
