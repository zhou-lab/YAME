#include <sys/stat.h>
#include <sys/types.h>
#include "cfile.h"
#include "snames.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: yame rowops [options] <in.cx> <out.cx>\n");
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

  cfile_t cf = open_cfile(fname);
  cdata_t c_out = {0};
  uint64_t i=0, k;
  for (k=0; ; ++k) {
    cdata_t c = read_cdata1(&cf);
    if (c.n == 0) break;

    if (c.fmt != '3') {
      fprintf(stderr, "[%s:%d] Only fmt 3 is supported. Got file format: %c\n", __func__, __LINE__, c.fmt);
      fflush(stderr);
      exit(1);
    }

    cdata_t c2 = {0};
    decompress(&c, &c2);

    if (!op || strcmp(op, "sum")==0) {
      if (k) {
        uint64_t *src = (uint64_t*) c2.s;
        uint64_t *dest = (uint64_t*) c_out.s;
        for (i=0; i<c2.n; ++i)
          dest[i] = sumMUpair(src[i], dest[i]);
      } else {
        c_out = c2;
        c_out.s = malloc(sizeof(uint64_t)*c2.n);
        memcpy(c_out.s, c2.s, sizeof(uint64_t)*c2.n);
      }
    } else if (strcmp(op, "binsum")==0) {
      if (!k) {
        c_out = c2;
        c_out.s = calloc(c2.n, sizeof(uint64_t));
      }
      uint64_t *src = (uint64_t*) c2.s;
      uint64_t *dest = (uint64_t*) c_out.s;
      for (i=0; i<c2.n; ++i)
        dest[i] += MUbinarize(src[i]);

    } else {
      fprintf(stderr, "[%s:%d] Unsupported operation: %s\n", __func__, __LINE__, op);
      fflush(stderr);
      exit(1);
    }

    free(c2.s); free(c.s);
  }
  cdata_write(fname_out, &c_out, "wb", verbose);
  free(c_out.s);
  bgzf_close(cf.fh);
  
  return 0;
}
