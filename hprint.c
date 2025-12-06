#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include "cfile.h"

#define ANSI_COLOR_YELLOW   "\x1b[33m"
#define ANSI_COLOR_BLUE     "\x1b[34m"
#define ANSI_COLOR_GREY     "\x1b[90m" // light background: \x1b[37m
#define ANSI_COLOR_RESET    "\x1b[0m"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: yame hprint [options] <binary.cg>\n");
  fprintf(stderr, "Print data transposed / horizontally.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -h         This help\n");
  fprintf(stderr, "\n");

  return 1;
}

int main_hprint(int argc, char *argv[]) {

  int c;
  while ((c = getopt(argc, argv, "h"))>=0) {
    switch (c) {
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 1 > argc) {
    usage(); 
    wzfatal("Please supply input file.\n");
  }

  char *fname = argv[optind];
  cfile_t cf = open_cfile(fname);
  while (1) {
    cdata_t c = read_cdata1(&cf);
    if (c.n == 0) break;
    decompress2(&c);
    if (c.fmt != '6') {
      fprintf(stderr, "[%s:%d] Only format 6 (given %d) files are supported.\n", __func__, __LINE__, c.fmt);
      fflush(stderr);
      exit(1);
    }

    for (uint64_t i=0; i<c.n; ++i) {
      if (FMT6_IN_UNI(c,i)) {
        if (FMT6_IN_SET(c,i)) {
          fprintf(stdout, ANSI_COLOR_YELLOW "1" ANSI_COLOR_RESET);
        } else {
          fprintf(stdout, ANSI_COLOR_BLUE "0" ANSI_COLOR_RESET);
        }
      } else {
        fprintf(stdout, ANSI_COLOR_GREY "2" ANSI_COLOR_RESET);
      }
    }
    fputc('\n', stdout);
    free_cdata(&c);
  }

  bgzf_close(cf.fh);
  return 0;
}
