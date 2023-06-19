#include <sys/stat.h>
#include <sys/types.h>
#include "kycg.h"
#include "snames.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: kycg rowops [options] <in.cg> <out.cg>\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -v        verbose\n");
  fprintf(stderr, "    -o        operations: {sum, binsum}\n");
  fprintf(stderr, "    -s        chunk size\n");
  fprintf(stderr, "    -h        This help\n");
  fprintf(stderr, "\n");

  return 1;
}

int main_rowops(int argc, char *argv[]) {

  int c, verbose = 0;
  char *op = NULL;
  while ((c = getopt(argc, argv, "s:o:vh"))>=0) {
    switch (c) {
    case 'v': verbose = 1; break;
    case 'o': op = strdup(optarg); break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 1 > argc) {
    usage(); 
    wzfatal("Please supply input file.\n");
  }

  char *fname = argv[optind];
  char *fname_out = NULL;
  if (argc >= optind + 2)
    fname_out = strdup(argv[optind+1]);

  cgfile_t cgf = open_cgfile(fname);
  cgdata_t cg_out = {0};
  uint64_t i=0, k;
  for (k=0; ; ++k) {
    cgdata_t cg = read_cg(&cgf);
    if (cg.n == 0) break;

    if (cg.fmt != '3') {
      fprintf(stderr, "[%s:%d] Only fmt 3 is supported. Got file format: %c\n", __func__, __LINE__, cg.fmt);
      fflush(stderr);
      exit(1);
    }

    cgdata_t cg2 = {0};
    decompress(&cg, &cg2);

    if (!op || strcmp(op, "sum")==0) {
      if (k) {
        uint64_t *src = (uint64_t*) cg2.s;
        uint64_t *dest = (uint64_t*) cg_out.s;
        for (i=0; i<cg2.n; ++i)
          dest[i] = sumMUpair(src[i], dest[i]);
      } else {
        cg_out = cg2;
        cg_out.s = malloc(sizeof(uint64_t)*cg2.n);
        memcpy(cg_out.s, cg2.s, sizeof(uint64_t)*cg2.n);
      }
    } else if (strcmp(op, "binsum")==0) {
      if (!k) {
        cg_out = cg2;
        cg_out.s = calloc(cg2.n, sizeof(uint64_t));
      }
      uint64_t *src = (uint64_t*) cg2.s;
      uint64_t *dest = (uint64_t*) cg_out.s;
      for (i=0; i<cg2.n; ++i)
        dest[i] += MUbinarize(src[i]);

    } else {
      fprintf(stderr, "[%s:%d] Unsupported operation: %s\n", __func__, __LINE__, op);
      fflush(stderr);
      exit(1);
    }

    free(cg2.s); free(cg.s);
  }
  cgdata_write(fname_out, &cg_out, "wb", verbose);
  free(cg_out.s);
  bgzf_close(cgf.fh);
  
  return 0;
}
