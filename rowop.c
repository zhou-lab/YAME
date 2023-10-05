#include <sys/stat.h>
#include <sys/types.h>
#include "cfile.h"
#include "snames.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: yame rowop [options] <in.cx> <out.cx>\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -v        verbose\n");
  fprintf(stderr, "    -o        operations: {binasum}\n");
  fprintf(stderr, "              binasum: sum binary data to M and U (format 3).\n");
  fprintf(stderr, "              musum: sum M and U separately (format 3).\n");
  fprintf(stderr, "              mean: mean beta and counts of data points.\n");
  fprintf(stderr, "    -c        minimum sequencing depth for rowops (default 1).\n");
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

static void binasumFmt3(cdata_t *cout, cdata_t *c, unsigned mincov) {
  for (uint64_t i=0; i<c->n; ++i) {
    uint64_t mu0 = f3_get_mu(c, i);
    if (!mu0) continue; // 0-0 is skipped
    if (((mu0>>32) + (mu0<<32>>32)) < mincov) continue;
    uint64_t mu = f3_get_mu(cout, i);
    if ((mu0>>32) > (mu0<<32>>32)) {
      f3_set_mu(cout, i, ((mu>>32)+1), (mu<<32>>32));
    } else {
      f3_set_mu(cout, i, (mu>>32), ((mu<<32>>32)+1));
    }
  }
}

static cdata_t rowop_binasum(cfile_t cf, unsigned mincov) {
  cdata_t c = read_cdata1(&cf);
  cdata_t cout = {0};
  if (c.n == 0) return cout;    // nothing in cfile
  char fmt = c.fmt;
  cout.n = cdata_n(&c);
  cout.compressed = 0;
  cout.fmt = '3';
  cout.unit = 8;                // max-size result
  cout.s = calloc(cout.n, sizeof(uint64_t));
  
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
    case '3': binasumFmt3(&cout, &c2, mincov); break;
    default: {
      fprintf(stderr, "[%s:%d] File format: %c unsupported.\n", __func__, __LINE__, c.fmt);
      fflush(stderr);
      exit(1);
    }}

    free(c.s); free(c2.s);
  }
  return cout;
}

static void musumFmt3(cdata_t *cout, cdata_t *c) {
  for (uint64_t i=0; i<c->n; ++i) {
    uint64_t mu0 = f3_get_mu(c, i);
    if (!mu0) continue; // 0-0 is skipped
    uint64_t mu = f3_get_mu(cout, i);
    f3_set_mu(cout, i, ((mu>>32)+(mu0>>32)), ((mu<<32>>32)+(mu0<<32>>32)));
  }
}

static cdata_t rowop_musum(cfile_t cf) {
  cdata_t c = read_cdata1(&cf);
  cdata_t cout = {0};
  if (c.n == 0) return cout;    // nothing in cfile
  char fmt = c.fmt;
  cout.n = cdata_n(&c);
  cout.compressed = 0;
  cout.fmt = '3';
  cout.unit = 8;                // max-size result
  cout.s = calloc(cout.n, sizeof(uint64_t));
  
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
    case '3': musumFmt3(&cout, &c2); break;
    default: {
      fprintf(stderr, "[%s:%d] File format: %c unsupported.\n", __func__, __LINE__, c.fmt);
      fflush(stderr);
      exit(1);
    }}

    free(c.s); free(c2.s);
  }
  return cout;
}

static void meanFmt3(uint32_t *cnts, double *fracs, cdata_t *c, unsigned mincov) {
  for (uint64_t i=0; i<c->n; ++i) {
    uint64_t mu0 = f3_get_mu(c, i);
    if (!mu0) continue; // 0-0 is skipped
    if (((mu0>>32) + (mu0<<32>>32)) < mincov) continue;
    uint64_t M = mu0>>32;
    uint64_t U = (mu0<<32>>32);
    if (M+U >= mincov) {
      fracs[i] += (double) M / (M+U);
      cnts[i]++;
    }
  }
}

static void rowop_mean(cfile_t cf, char *fname_out, unsigned mincov) {
  cdata_t c = read_cdata1(&cf);
  if (c.n == 0) return;    // nothing in cfile
  uint64_t n = cdata_n(&c);
  uint32_t *cnts = calloc(n, sizeof(uint32_t));
  double *fracs = calloc(n, sizeof(double));
  
  for (uint64_t k=0; ; ++k) {
    if (k) c = read_cdata1(&cf); // skip 1st cdata
    if (c.n == 0) break;
    cdata_t c2 = {0};
    decompress(&c, &c2);
    
    switch (c.fmt) {
    case '3': meanFmt3(cnts, fracs, &c2, mincov); break;
    default: {
      fprintf(stderr, "[%s:%d] File format: %c unsupported.\n", __func__, __LINE__, c.fmt);
      fflush(stderr);
      exit(1);
    }}

    free(c.s); free(c2.s);
  }

  FILE *out;
  if (fname_out) {
    out = fopen(fname_out, "w");
  } else {
    out = stdout;
  }
  for (uint64_t i=0; i<n; ++i) {
    if (!cnts[i]) {
      fputs("NA\t0\n", out);
    } else {
      fprintf(out, "%1.3f\t%d\n", fracs[i] / cnts[i], cnts[i]);
    }
  }
  if (fname_out) fclose(out);
}

int main_rowop(int argc, char *argv[]) {

  int c, verbose = 0;
  unsigned mincov = 1;
  char *op = NULL;
  while ((c = getopt(argc, argv, "vo:c:h"))>=0) {
    switch (c) {
    case 'v': verbose = 1; break;
    case 'o': op = strdup(optarg); break;
    case 'c': mincov = atoi(optarg); break;
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
    cout = rowop_binasum(cf, mincov);
    cdata_write(fname_out, &cout, "wb", verbose);
    free(cout.s);
  } else if (strcmp(op, "mean") == 0) {
    rowop_mean(cf, fname_out, mincov);
  } else if (strcmp(op, "musum") == 0) {
    cout = rowop_musum(cf);
    cdata_write(fname_out, &cout, "wb", verbose);
    free(cout.s);
  } else {
    fprintf(stderr, "[%s:%d] Unsupported operation: %s\n", __func__, __LINE__, op);
    fflush(stderr);
    exit(1);
  }
  bgzf_close(cf.fh);
  if (fname_out) free(fname_out);
  
  return 0;
}


