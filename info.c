#include "cfile.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: yame info [options] <in.cx>\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -h        This help\n");
  fprintf(stderr, "\n");

  return 1;
}

int main_info(int argc, char *argv[]) {

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

  char *fname_in = argv[optind];
  cfile_t cf = open_cfile(fname_in);
  snames_t snames = loadSampleNamesFromIndex(fname_in);
  int i = 0;
  fprintf(stdout, "Sample\tN\tFormat\tUnitBytes\n");
  for (i=0; ; ++i) {
    cdata_t c = read_cdata1(&cf);
    if (c.n == 0) break;
    cdata_t expanded = {0};
    decompress(&c, &expanded);
    if (snames.n) {
      fputs(snames.s[i], stdout);
    } else {
      fprintf(stdout, "%d", i+1);
    }
    fprintf(stdout, "\t%"PRIu64"\t%c\t%u\n", expanded.n, expanded.fmt, expanded.unit);
    free(expanded.s); free(c.s);
  }
  bgzf_close(cf.fh);
  
  return 0;
}
