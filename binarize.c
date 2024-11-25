#include <sys/stat.h>
#include <sys/types.h>
#include "cfile.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: yame binarize [options] <mu.cg>\n");
  fprintf(stderr, "If Beta>0, then 1, else 0. M+U>0 is used as universe.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -t [T]    1 if Beta>T else 0. (default: 0.5).\n");
  fprintf(stderr, "    -o        output cx file name (format 6). if missing, output to stdout.\n");
  fprintf(stderr, "    -h        This help\n");
  fprintf(stderr, "\n");

  return 1;
}

int main_binarize(int argc, char *argv[]) {

  int c; double T = 0.5;
  char *fname_out = NULL;
  while ((c = getopt(argc, argv, "o:t:h"))>=0) {
    switch (c) {
    case 'o': fname_out = strdup(optarg); break;
    case 't': T = atof(optarg); break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 1 > argc) {
    usage(); 
    wzfatal("Please supply input file.\n");
  }

  char *fname = argv[optind];

  BGZF *fp_out;
  if (fname_out) fp_out = bgzf_open2(fname_out, "w");
  else fp_out = bgzf_dopen(fileno(stdout), "w");
  if (fp_out == NULL) {
    fprintf(stderr, "Error opening file for writing: %s\n", fname_out);
    exit(1);
  }

  cfile_t cf = open_cfile(fname);
  while (1) {
    cdata_t c = read_cdata1(&cf);
    if (c.n == 0) break;
    decompress2(&c);
    if (c.fmt != '3') {
      fprintf(stderr, "[%s:%d] Only format %d files are supported.\n", __func__, __LINE__, c.fmt);
      fflush(stderr);
      exit(1);
    }

    cdata_t c6 = {.fmt = '6', .n = c.n};
    c6.s = calloc((c6.n+3)/4, sizeof(uint8_t));
    for (uint64_t i=0; i<c6.n; ++i) {
      uint64_t mu = f3_get_mu(&c, i);
      if (mu>0) {
        if (MU2beta(mu)>T) FMT6_SET1(c6, i);
        else FMT6_SET0(c6, i);
      }
    }
    cdata_compress(&c6);
    cdata_write1(fp_out, &c6);
    free_cdata(&c6);
    free_cdata(&c);
  }

  if (fname_out) free(fname_out);
  bgzf_close(fp_out);
  bgzf_close(cf.fh);
  return 0;
}
