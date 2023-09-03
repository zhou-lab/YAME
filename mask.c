#include <sys/stat.h>
#include <sys/types.h>
#include "cfile.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: yame mask [options] <in.cg> <mask.cx>\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -v        verbose\n");
  fprintf(stderr, "    -h        This help\n");
  fprintf(stderr, "\n");

  return 1;
}

static void prepare_mask(cdata_t *c) {
  if (c->fmt < '2') {
    convertToFmt0(c);
  } else {
    decompress2(c);
  }
}

int main_mask(int argc, char *argv[]) {

  int c, verbose = 0;
  while ((c = getopt(argc, argv, "vh"))>=0) {
    switch (c) {
    case 'v': verbose = 1; break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 2 > argc && verbose) {
    usage(); 
    wzfatal("Please supply input file.\n");
  }

  char *fname = argv[optind];
  char *fname_mask = argv[optind+1];

  cfile_t cf_mask = open_cfile(fname_mask);
  cdata_t c_mask = read_cdata1(&cf_mask);
  prepare_mask(&c_mask);
  BGZF *fp_out = bgzf_dopen(fileno(stdout), "w");

  cfile_t cf = open_cfile(fname);
  while (1) {
    cdata_t c = read_cdata1(&cf);
    if (c.n == 0) break;
    decompress2(&c);

    if (c.n != c_mask.n) {
      fprintf(stderr, "[%s:%d] mask (n=%"PRIu64") and query (N=%"PRIu64") are of different lengths.\n", __func__, __LINE__, c_mask.n, c.n);
      fflush(stderr);
      exit(1);
    }

    if (c.fmt != '3') {
      fprintf(stderr, "[%s:%d] Only format 3 files are supported.\n", __func__, __LINE__);
      fflush(stderr);
      exit(1);
    }
    
    for (uint64_t i=0; i<c.n; ++i) {
      if (c_mask.s[i>>3]&(1<<(i&0x7))) {
        f3_set_mu(&c, i, 0, 0);
      }
    }
    cdata_compress(&c);
    cdata_write1(fp_out, &c);
    free_cdata(&c);
  }

  bgzf_close(fp_out);
  free_cdata(&c_mask);
  bgzf_close(cf_mask.fh);
  return 0;
}
