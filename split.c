#include "kycg.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: kycg split [options] <in.cg> out_prefix\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -v        verbose\n");
  fprintf(stderr, "    -h        This help\n");
  fprintf(stderr, "\n");

  return 1;
}

int main_split(int argc, char *argv[]) {

  int c, verbose = 0;
  while ((c = getopt(argc, argv, "vh"))>=0) {
    switch (c) {
    case 'v': verbose = 1; break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 2 > argc) { 
    usage(); 
    wzfatal("Please supply input file.\n");
  }

  cgfile_t cgf = open_cgfile(argv[optind++]);
  char *prefix = argv[optind];

  int i = 0;
  for (i=0; ; ++i) {
    cgdata_t cg = read_cg(&cgf);
    if (cg.n == 0) break;
    char *tmp = malloc(strlen(prefix) + 1000);
    sprintf(tmp, "%s_split_%i.cg", prefix, i+1);
    cgdata_write(tmp, &cg, verbose);
  }
  
  return 0;
}
