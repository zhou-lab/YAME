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
  fprintf(stderr, "    -o        operations: {binasum}\n");
  fprintf(stderr, "              binasum: sum binary data to M and U (format 3).\n");
  fprintf(stderr, "    -s        chunk size\n");
  fprintf(stderr, "    -h        This help\n");
  fprintf(stderr, "\n");

  return 1;
}

static void binasumFmt0(cdata_t *cout, cdata_t *c) {
  for (uint64_t i=0; i<c->n; ++i) {
    uint64_t mu = f3_get_mu(cout, i);
    if (c->s[i]) {
      f3_set_mu(cout, i, ((mu>>32)+1), (mu<<32>>32));
    } else {
      f3_set_mu(cout, i, (mu>>32), ((mu<<32>>32)+1));
    }
  }
}

static void binasumFmt3(cdata_t *cout, cdata_t *c) {
  for (uint64_t i=0; i<c->n; ++i) {
    uint64_t mu0 = f3_get_mu(c, i);
    uint64_t mu = f3_get_mu(cout, i);
    if ((mu0>>32) > (mu0<<32>>32)) {
      f3_set_mu(cout, i, ((mu>>32)+1), (mu<<32>>32));
    } else {
      f3_set_mu(cout, i, (mu>>32), ((mu<<32>>32)+1));
    }
  }
}

static cdata_t rowops_binasum(cfile_t cf) {
  cdata_t c = read_cdata1(&cf);
  cdata_t cout = {0};
  if (c.n == 0) return cout;    // nothing in cfile
  char fmt = c.fmt;
  cout.n = cdata_n(&c);
  cout.compressed = 0;
  cout.fmt = '3';
  cout.unit = 8;                // max-size result
  cout.s = malloc(sizeof(uint64_t)*cout.n);
  
  for (uint64_t k=0; ; ++k) {
    if (k) c = read_cdata1(&cf); // skip 1st cdata
    if (c.n == 0) break;
    if (fmt != c.fmt) {
      fprintf(stderr, "[%s:%d] File formats are inconsistent: %c vs %c.\n", __func__, __LINE__, fmt, c.fmt);
      fflush(stderr);
      exit(1);
    }
    cdata_t c2 = {0};
    decompress(&c, &c2);
    if (c2.n != cout.n) {
      fprintf(stderr, "[%s:%d] Data dimensions are inconsistent: %"PRIu64" vs %"PRIu64"\n", __func__, __LINE__, cout.n, c2.n);
      fflush(stderr);
      exit(1);
    }
    
    switch (fmt) {
    case '0': binasumFmt0(&cout, &c2); break;
    case '1': binasumFmt0(&cout, &c2); break;
    case '3': binasumFmt3(&cout, &c2); break;
    default: {
      fprintf(stderr, "[%s:%d] File format: %c unsupported.\n", __func__, __LINE__, c.fmt);
      fflush(stderr);
      exit(1);
    }}

    free(c.s); free(c2.s);
  }
  return cout;
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
  cdata_t cout = {0};
  if (!op || strcmp(op, "binasum") == 0) { // default
    cout = rowops_binasum(cf);
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


