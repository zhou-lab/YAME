#include <sys/stat.h>
#include <sys/types.h>
#include "cfile.h"
#include "snames.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: yame rowop [options] <in.cx> <out>\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -o        Operations (choose one):\n");
  fprintf(stderr, "              binasum     Sum binarized data to M and U (format 3).\n");
  fprintf(stderr, "                          Can also take format 0 or 1 binary data.\n");
  fprintf(stderr, "                          Output: new cx file.\n");
  fprintf(stderr, "              musum       Sum M and U separately (format 3).\n");
  fprintf(stderr, "                          Output: new cx file.\n");
  fprintf(stderr, "              mean        Mean beta and counts of data points (format 3).\n");
  fprintf(stderr, "                          Output: plain text (two columns).\n");
  fprintf(stderr, "              std         Standard deviation. Requires format 3 cx.\n");
  fprintf(stderr, "                          Output: plain text (std, counts).\n");
  fprintf(stderr, "              binstring   Binarize data to row-wise string (format 3).\n");
  fprintf(stderr, "                          Output: plain text file with binary strings.\n");
  fprintf(stderr, "              cometh      Co-methylation of neighboring CGs.\n");
  fprintf(stderr, "                          Output: plain text in uint64_t U0U1-U0M1-M0U1-M0M1,\n");
  fprintf(stderr, "                          U0U2-U0M2-M0U2-M0M2, etc. '0' is target CG,\n");
  fprintf(stderr, "                          followed by 1, 2, etc for neighboring CGs.\n");
  fprintf(stderr, "                          Each pair occupies 16 bits. For visual, use -v.\n");
  fprintf(stderr, "                          Intermediate methylations (0.3-0.7) are excluded.\n");
  fprintf(stderr, "    -w        Number of neighboring CGs for cometh (default: 5).\n");
  fprintf(stderr, "    -c        Minimum sequencing depth for rowops (default 1).\n");
  fprintf(stderr, "    -v        Verbose mode\n");
  fprintf(stderr, "    -h        Display this help message\n\n");

  return 1;
}

static void binasumFmt0(cdata_t *cout, cdata_t *c) {
  for (uint64_t i=0; i<c->n; ++i) {
    uint64_t mu = f3_get_mu(cout, i);
    if (FMT0_IN_SET(*c, i)) {
      f3_set_mu(cout, i, ((mu>>32)+1), (mu<<32>>32));
    } else {
      f3_set_mu(cout, i, (mu>>32), ((mu<<32>>32)+1));
    }
  }
}

static void binasumFmt1(cdata_t *cout, cdata_t *c) {
  for (uint64_t i=0; i<c->n; ++i) {
    uint64_t mu = f3_get_mu(cout, i);
    if (c->s[i]-'0') {
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
    case '1': binasumFmt1(&cout, &c2); break;
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
  free(cnts); free(fracs);
  if (fname_out) fclose(out);
}

static void sumsqFmt3(uint32_t *cnts, double *sum, double *sum_sq, cdata_t *c, unsigned mincov) {
  for (uint64_t i=0; i<c->n; ++i) {
    uint64_t mu0 = f3_get_mu(c, i);
    if (!mu0) continue; // 0-0 is skipped
    if (((mu0>>32) + (mu0<<32>>32)) < mincov) continue;
    uint64_t M = mu0>>32;
    uint64_t U = (mu0<<32>>32);
    if (M+U >= mincov) {
      double x = (double) M / (M+U);
      sum[i] += x;
      sum_sq[i] += x * x;
      cnts[i]++;
    }
  }
}


// the following doesn't work for large numbers but should be ok for meth levels
// see https://www.strchr.com/standard_deviation_in_one_pass
static void rowop_std(cfile_t cf, char *fname_out, unsigned mincov) {

  cdata_t c = read_cdata1(&cf);
  if (c.n == 0) return; // nothing in cfile, output nothing
  uint64_t n = cdata_n(&c);
  uint32_t *cnts = calloc(n, sizeof(uint32_t));
  double *sum = calloc(n, sizeof(double));
  double *sum_sq = calloc(n, sizeof(double));
  
  for (uint64_t k = 0; ; ++k) {
    if (k) c = read_cdata1(&cf); // skip 1st cdata
    if (c.n == 0) break;
    cdata_t c2 = {0};
    decompress(&c, &c2);

    switch (c.fmt) {
    case '3': sumsqFmt3(cnts, sum, sum_sq, &c2, mincov); break;
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
      double mean = sum[i] / cnts[i];
      fprintf(out, "%1.3f\t%d\n", sqrt((sum_sq[i] / cnts[i]) - mean * mean), cnts[i]);
    }
  }
  free(cnts); free(sum_sq); free(sum);
  if (fname_out) fclose(out);
}

static void rowop_binstring(cfile_t cf, char *fname_out) {
  cdata_t c = read_cdata1(&cf);
  if (c.n == 0) return;    // nothing in cfile
  uint64_t n = cdata_n(&c);
  uint64_t binstring_bytes = 0;
  uint8_t *binstring = NULL;
  uint64_t k=0;
  for (k=0; ; ++k) {
    if (k) c = read_cdata1(&cf); // skip 1st cdata
    if (c.n == 0) break;
    cdata_t c2 = {0};
    decompress(&c, &c2);

    if (binstring_bytes*8 <= k) {
      binstring_bytes++;
      binstring = realloc(binstring, (binstring_bytes*n));
      memset(binstring + (binstring_bytes-1)*n, 0, sizeof(n));
    }
    
    switch (c.fmt) {
    case '3': {
      for (uint64_t i=0; i<c2.n; ++i) {
        uint64_t mu = f3_get_mu(&c2, i);
        if ((mu>>32) > (mu<<32>>32)) {
          binstring[(k>>3)*n+i] |= (1<<(k&0x7));
        }
      }
      break;
    }
    default: {
      fprintf(stderr, "[%s:%d] File format: %c unsupported.\n", __func__, __LINE__, c.fmt);
      fflush(stderr);
      exit(1);
    }}

    free(c.s); free(c2.s);
  }

  FILE *out;
  if (fname_out) { out = fopen(fname_out, "w");
  } else { out = stdout; }
  for (uint64_t i=0; i<n; ++i) {
    for (uint64_t kk=0; kk<k; ++kk) {
      fputc('0'+((binstring[(kk>>3)*n+i] >> (kk&0x7))&0x1), out);
    }
    fputc('\n', out);
  }
  free(binstring);
  if (fname_out) fclose(out);
}

void rowop_cometh(cfile_t cf, char *fname_out, unsigned mincov, int cometh_window, int verbose) {

  uint64_t *cnts = NULL; uint64_t ncnts = 0;
  for (uint64_t k=0; ;++k) {
    cdata_t c0 = read_cdata1(&cf);
    if (c0.n == 0) break;
    cdata_t c = {0};
    decompress(&c0, &c);
    if (!k) {                   /* first data, initialize */
      cnts = calloc(c.n*cometh_window, sizeof(uint64_t));
      ncnts = c.n;
    }
    assert(c.fmt == '3');
    for (uint64_t i=0; i<ncnts-cometh_window; ++i) {
      for (uint64_t j=i+1; j<=min(ncnts-1, i+cometh_window); ++j) {
        uint64_t mu = f3_get_mu(&c, i);
        uint64_t M = mu>>32; uint64_t U = (mu<<32>>32);
        uint64_t mu1 = f3_get_mu(&c, j);
        uint64_t M1 = mu1>>32; uint64_t U1 = (mu1<<32>>32);
        if (M+U >= mincov && M1+U1 >= mincov) {
          // also skip intermediate values too close to 0.5
          double b = M/(M+U); double b1 = M1/(M1+U1);
          if (fabs(b - 0.5) >= 0.2 || fabs(b1 - 0.5) >= 0.2) {
            int shift = (M<U?2:0) + (M1<U1?1:0);
            cnts[i*cometh_window+(j-i-1)] += (1ul<<(shift*16));
          }
        }
      }
    }
    free(c0.s); free(c.s);
  }

  FILE *out;
  if (fname_out) out = fopen(fname_out, "w");
  else out = stdout;
  for (uint64_t i=0; i<ncnts; ++i) {
    fprintf(out, "%"PRIu64, i+1);
    for (uint64_t j=0; j<(unsigned) cometh_window; ++j) {
      uint64_t data = cnts[i*cometh_window+j];
      fputc('\t', out);
      if (verbose) {
        fprintf(out, "%"PRIu64"-%"PRIu64"-%"PRIu64"-%"PRIu64,
                data>>(16*3), data<<(16)>>(16*3), data<<(16*2)>>(16*3), data<<(16*3)>>(16*3));
      } else {
        fprintf(out, "%"PRIu64, data);
      }
    }
    fputc('\n', out);
  }
  if (ncnts) free(cnts);
  if (fname_out) fclose(out);
}

int main_rowop(int argc, char *argv[]) {

  int c, verbose = 0;
  unsigned mincov = 1;
  char *op = NULL; int cometh_window = 5;
  while ((c = getopt(argc, argv, "vo:c:w:h"))>=0) {
    switch (c) {
    case 'v': verbose = 1; break;
    case 'o': op = strdup(optarg); break;
    case 'c': mincov = atoi(optarg); break;
    case 'w': cometh_window = atoi(optarg); break;
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
  } else if (strcmp(op, "std") == 0) {
    rowop_std(cf, fname_out, mincov);
  } else if (strcmp(op, "musum") == 0) {
    cout = rowop_musum(cf);
    cdata_write(fname_out, &cout, "wb", verbose);
    free(cout.s);
  } else if (strcmp(op, "binstring") == 0) {
    rowop_binstring(cf, fname_out);
  } else if (strcmp(op, "cometh") == 0) {
    rowop_cometh(cf, fname_out, mincov, cometh_window, verbose);
  } else {
    fprintf(stderr, "[%s:%d] Unsupported operation: %s\n", __func__, __LINE__, op);
    fflush(stderr);
    exit(1);
  }
  bgzf_close(cf.fh);
  if (fname_out) free(fname_out);
  free(op);
  
  return 0;
}


