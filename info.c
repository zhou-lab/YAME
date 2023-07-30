#include "cfile.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: yame info [options] <in.cx>\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -1        Report one record per file.\n");
  fprintf(stderr, "    -h        This help\n");
  fprintf(stderr, "\n");

  return 1;
}

int main_info(int argc, char *argv[]) {

  int c; int report1 = 0;
  while ((c = getopt(argc, argv, "1hv"))>=0) {
    switch (c) {
    case '1': report1 = 1; break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 1 > argc) { 
    usage(); 
    wzfatal("Please supply input file.\n"); 
  }

  fprintf(stdout, "File\tSample\tN\tFormat\tUnitBytes\n");
  for (int j = optind; j < argc; ++j) {
    char *fname_in = argv[j];
    cfile_t cf = open_cfile(fname_in);
    snames_t snames = loadSampleNamesFromIndex(fname_in);
    int i = 0;
    for (i=0; ; ++i) {
      cdata_t c = read_cdata1(&cf);
      if (c.n == 0) break;
      cdata_t expanded = {0};
      decompress(&c, &expanded);
      fprintf(stdout, "%s\t", fname_in);
      if (snames.n) {
        fputs(snames.s[i], stdout);
      } else {
        fprintf(stdout, "%d", i+1);
      }
      fprintf(stdout, "\t%"PRIu64"\t%c\t%u\n", expanded.n, expanded.fmt, expanded.unit);
      free(expanded.s); free(c.s);
      if (report1) break;
    }
    cleanSampleNames2(snames);
    bgzf_close(cf.fh);
  }
  return 0;
}
