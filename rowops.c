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

static void sumFmt0(cdata_t *cout, cdata_t *c) {

  
  
  if (k) {
    f3_unpack_mu(&c2, i);
    uint64_t *src = (uint64_t*) c2.s;
    uint64_t *dest = (uint64_t*) cout.s;
    for (i=0; i<c2.n; ++i)
      dest[i] = sumMUpair(src[i], dest[i]);
  } else {
    cout = c2;
    cout.s = malloc(sizeof(uint64_t)*c2.n);
    memcpy(cout.s, c2.s, sizeof(uint64_t)*c2.n);
  }
  free(c2.s);
}

static cdata_t rowops_sum(cfile_t cf) {
  cdata_t c = read_cdata1(&cf);
  if (c.n == 0) return cout; // nothing in cfile
  char fmt = c.fmt;
  uint64_t n = cdata_n(&c);
  cdata_t cout = {0};
  cout = malloc(sizeof(uint64_t)*n);
  
  uint64_t i=0, k;
  for (k=0; ; ++k) {
    if (k) c = read_cdata1(&cf); // skip 1st cdata
    if (c.n == 0) break;
    if (fmt != c.fmt) {
      fprintf(stderr, "[%s:%d] File formats are inconsistent: %c vs %c.\n", __func__, __LINE__, fmt, c.fmt);
      fflush(stderr);
      exit(1);
    }
    cdata_t c2 = {0};
    decompress(&c, &c2);
    if (c2.n != n) {
      fprintf(stderr, "[%s:%d] Data dimensions are inconsistent: %"PRIu64" vs %"PRIu64"\n", __func__, __LINE__, c, c2.n);
      fflush(stderr);
      exit(1);
    }
    
    switch (fmt) {
    case '0': {
      sumFmt0(&cout, &c2);
      break;
    }
    case '1': {
      sumFmt1(&cout, &c2);
      break;
    }
    case '3': {
      sumFmt3(&cout, &c2);
      break;
    }
    default: {
      fprintf(stderr, "[%s:%d] File format: %c unsupported.\n", __func__, __LINE__, c.fmt);
      fflush(stderr);
      exit(1);
    }}

    free(c.s); free(c2.s);
  }
  return cout;
}

/* if (!k) { */
/*         cout = c2; */
/*         cout.s = calloc(c2.n, sizeof(uint64_t)); */
/*       } */
/*       uint64_t *src = (uint64_t*) c2.s; */
/*       uint64_t *dest = (uint64_t*) cout.s; */
/*       for (i=0; i<c2.n; ++i) */
/*         dest[i] += MUbinarize(src[i]); */

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
  if (!op || strcmp(op, "sum") == 0) { // default
    cout = rowops_sum(cf);
  } else if (strcmp(op, "binsum") == 0) {
    cout = rowops_binsum(cf);
  } else {
    fprintf(stderr, "[%s:%d] Unsupported operation: %s\n", __func__, __LINE__, op);
    fflush(stderr);
    exit(1);
  }
  cdata_write(fname_out, &cout, "wb", verbose);
  free(cout.s);
  bgzf_close(cf.fh);
  
  return 0;
}


