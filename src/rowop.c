// SPDX-License-Identifier: AGPL-3.0-or-later
/**
 * This file is part of YAME.
 *
 * Copyright (C) 2021-present Wanding Zhou
 *
 * YAME is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * YAME is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with YAME.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include "cfile.h"
#include "snames.h"

/**
 * yame rowop
 * ==========
 *
 * Overview
 * --------
 * rowop performs row-wise aggregation across multiple records (samples) in a CX file.
 * Records are read sequentially, and most operations assume:
 *   - identical format across records
 *   - identical row dimension across records
 *
 * Operations (-o; default = "binasum")
 * -----------------------------------
 *
 * 1) binasum  (CX output; fmt3)
 *    Purpose:
 *      Convert per-sample values into per-row sample counts:
 *        M = #samples called methylated
 *        U = #samples called unmethylated
 *      (counts reflect samples, not sequencing depth).
 *
 *    Supported inputs:
 *      - fmt0: bitset (1->M, 0->U)
 *      - fmt1: ASCII '0'/'1' (nonzero->M, zero->U)
 *      - fmt3: MU counts; beta thresholds (-p/-q) define calls:
 *          * skip if mu==0 or cov < mincov
 *          * beta > beta1 => M++
 *          * beta < beta0 => U++
 *          * otherwise ignored
 *
 * 2) musum  (CX output; fmt3)
 *    Purpose:
 *      Sum MU sequencing counts across samples.
 *    Input:
 *      fmt3 only.
 *
 * 3) stat  (text output; fmt3)
 *    Purpose:
 *      Compute per-row summary statistics across samples.
 *    Filters:
 *      skip mu==0 and cov < mincov.
 *    Output columns:
 *      count        number of contributing samples
 *      mean_beta   mean beta value across samples
 *      sd_beta     standard deviation of beta
 *      delta_beta  min(beta>0.5) - max(beta<0.5)
 *      min_n       min(#beta<0.5, #beta>0.5)
 *
 *    Notes:
 *      - delta_beta is only defined when both sides exist; otherwise printed as NA.
 *      - sd is computed as sqrt(E[x^2] - E[x]^2).
 *
 * 4) binstring  (text output; fmt3)
 *    Purpose:
 *      Emit a row-wise binary string across samples.
 *    Behavior:
 *      For each sample/row, output '1' if beta > beta_threshold (-b), else '0'.
 *    Notes:
 *      Current implementation checks mu!=0 but does not enforce mincov.
 *
 * 5) cometh  (text output; fmt3)
 *    Purpose:
 *      Summarize co-methylation between each row and its neighbors (i+1..i+W).
 *    Behavior:
 *      - requires cov >= mincov at both sites
 *      - skips intermediate methylation near 0.5 (|beta-0.5| < 0.2)
 *    Output:
 *      One line per row, with packed 4-way counts (UU, UM, MU, MM) per neighbor.
 *      With -v, counts are printed as "UU-UM-MU-MM".
 *
 * I/O
 * ---
 * - <in.cx> is required.
 * - [out] is optional:
 *     * text operations write to stdout if omitted
 *     * CX-output operations write to stdout via cdata_write() when out is NULL
 */

typedef struct config_rowop_t {
  double beta0; // lower threshold
  double beta1; // higher threshold
  unsigned mincov;
  double beta_threshold;   // default to 0.5
  int cometh_window;
  int verbose;
  unsigned seed;
} config_rowop_t;

static int usage(void) {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "  yame rowop [options] <in.cx> [out]\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Purpose:\n");
  fprintf(stderr, "  Perform row-wise operations across multiple records (samples) in a CX file.\n");
  fprintf(stderr, "  Depending on the operation, output is either a new CX file or plain text.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Operation:\n");
  fprintf(stderr, "  -o <op>      Operation name (default: binasum)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "CX-output operations:\n");
  fprintf(stderr, "  binasum      Convert per-sample values into per-row sample counts (M/U) as format 3.\n");
  fprintf(stderr, "              Input: fmt0, fmt1, or fmt3.\n");
  fprintf(stderr, "              For fmt3, beta thresholds (-p/-q) define methylated vs unmethylated calls.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  musum        Sum MU sequencing counts across samples.\n");
  fprintf(stderr, "              Input: fmt3 only. Output: one fmt3 record.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Text-output operations:\n");
  fprintf(stderr, "  stat         Per-row summary statistics across samples.\n");
  fprintf(stderr, "              Input: fmt3 only.\n");
  fprintf(stderr, "              Output columns:\n");
  fprintf(stderr, "                count  mean_beta  sd_beta  delta_beta  min_n\n");
  fprintf(stderr, "              delta_beta = min(beta>0.5) - max(beta<0.5).\n");
  fprintf(stderr, "              min_n      = min(#beta<0.5, #beta>0.5).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  binstring    Convert per-sample beta values into row-wise binary strings.\n");
  fprintf(stderr, "              Input: fmt3 only. Uses -b as the beta threshold.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  cometh       Neighbor co-methylation summary within a window.\n");
  fprintf(stderr, "              Input: fmt3 only.\n");
  fprintf(stderr, "              Output: packed 4-way counts (UU, UM, MU, MM) per neighbor offset.\n");
  fprintf(stderr, "              Use -v to print unpacked lanes.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Common filters:\n");
  fprintf(stderr, "  -c <mincov>  Minimum coverage (M+U) for a sample/row to contribute (default: 1).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "binasum (fmt3 input) thresholds:\n");
  fprintf(stderr, "  -p <beta0>   Call unmethylated if beta < beta0 (default: 0.4).\n");
  fprintf(stderr, "  -q <beta1>   Call methylated   if beta > beta1 (default: 0.6).\n");
  fprintf(stderr, "              Betas in [beta0, beta1] are ignored.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "binstring threshold:\n");
  fprintf(stderr, "  -b <beta>    Call methylated if beta > threshold (default: 0.5).\n");
  fprintf(stderr, "  -s [int]     Seed for tie breaking (default: current time).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "cometh options:\n");
  fprintf(stderr, "  -w <W>       Neighbor window size (default: 5).\n");
  fprintf(stderr, "  -v           Verbose output (print UU-UM-MU-MM instead of packed uint64).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Other:\n");
  fprintf(stderr, "  -h           Show this help message.\n");
  fprintf(stderr, "\n");
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

static void binasumFmt3(cdata_t *cout, cdata_t *c, config_rowop_t *cfg) {
  for (uint64_t i=0; i<c->n; ++i) {
    uint64_t mu = f3_get_mu(c, i);
    if (!mu) continue; // 0-0 is skipped
    if (MU2cov(mu) < cfg->mincov) continue;

    double beta = MU2beta(mu);
    uint64_t mu_out = f3_get_mu(cout, i);
    if (beta > cfg->beta1) {
      f3_set_mu(cout, i, ((mu_out>>32)+1), (mu_out<<32>>32));
    } else if (beta < cfg->beta0) {
      f3_set_mu(cout, i, (mu_out>>32), ((mu_out<<32>>32)+1));
    }
  }
}

static cdata_t rowop_binasum(cfile_t cf, config_rowop_t *cfg) {
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
    cdata_t c2 = decompress(c);
    if (c2.n != cout.n) {
      fprintf(stderr, "[%s:%d] Data dimensions are inconsistent: %"PRIu64" vs %"PRIu64"\n", __func__, __LINE__, cout.n, c2.n);
      fflush(stderr);
      exit(1);
    }
    
    switch (fmt) {
    case '0': binasumFmt0(&cout, &c2); break;
    case '1': binasumFmt1(&cout, &c2); break;
    case '3': binasumFmt3(&cout, &c2, cfg); break;
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
    cdata_t c2 = decompress(c);
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

static void collect_stat_fmt3(uint32_t *cnts, double *sum, double *sum_sq, double *b0max, double *b1min, int *b0n, int *b1n, cdata_t *c, config_rowop_t *cfg) {
  for (uint64_t i=0; i<c->n; ++i) {
    uint64_t mu0 = f3_get_mu(c, i);
    if (!mu0) continue; // 0-0 is skipped
    if (((mu0>>32) + (mu0<<32>>32)) < cfg->mincov) continue;
    uint64_t M = mu0>>32;
    uint64_t U = (mu0<<32>>32);
    if (M+U >= cfg->mincov) {
      double x = (double) M / (M+U);
      sum[i] += x;
      sum_sq[i] += x * x;
      cnts[i]++;
      if (x < 0.5) {
        b0n[i]++;
        if (x > b0max[i])
          b0max[i] = x;
      }
      if (x > 0.5) {
        b1n[i]++;
        if (x < b1min[i])
          b1min[i] = x;
      }
    }
  }
}


// the following standard deviation doesn't work for large numbers but should be ok for meth levels
// see https://www.strchr.com/standard_deviation_in_one_pass
static void rowop_stat(cfile_t cf, char *fname_out, config_rowop_t *cfg) {

  cdata_t c = read_cdata1(&cf);
  if (c.n == 0) return; // nothing in cfile, output nothing
  uint64_t n = cdata_n(&c);
  uint32_t *cnts = calloc(n, sizeof(uint32_t));
  double *sum = calloc(n, sizeof(double));
  double *sum_sq = calloc(n, sizeof(double));
  double *b0max = calloc(n, sizeof(double));
  double *b1min = calloc(n, sizeof(double));
  int *b0n = calloc(n, sizeof(int));
  int *b1n = calloc(n, sizeof(int));

  srand(cfg->seed);
  for (uint64_t i = 0; i < n; ++i) b1min[i] = 1.0;
  
  for (uint64_t k = 0; ; ++k) {
    if (k) c = read_cdata1(&cf); // skip 1st cdata
    if (c.n == 0) break;
    cdata_t c2 = decompress(c);

    switch (c.fmt) {
    case '3': collect_stat_fmt3(cnts, sum, sum_sq, b0max, b1min, b0n, b1n, &c2, cfg); break;
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

  fputs("count\tmean_beta\tsd_beta\tdelta_beta\tmin_n\n", out);
  for (uint64_t i = 0; i < n; ++i) {
    if (cnts[i] == 0) {
      fputs("0\tNA\tNA\tNA\t0\n", out);
      continue;
    }

    double mean = sum[i] / cnts[i];
    double sd   = sqrt((sum_sq[i] / cnts[i]) - mean * mean);

    /* delta_beta = b1min - b0max, but only meaningful if both sides exist */
    double delta_beta = (b0n[i] > 0 && b1n[i] > 0) ? (b1min[i] - b0max[i]) : -1.0;

    /* min_n = min(#beta<0.5, #beta>0.5) */
    uint32_t min_n = (b1n[i] < b0n[i]) ? b1n[i] : b0n[i];

    if (delta_beta < 0) {
      fprintf(out, "%u\t%1.3f\t%1.3f\tNA\t%u\n", cnts[i], mean, sd, min_n);
    } else {
      fprintf(out, "%u\t%1.3f\t%1.3f\t%1.3f\t%u\n", cnts[i], mean, sd, delta_beta, min_n);
    }
  }
  free(cnts);
  free(sum_sq);
  free(sum);
  if (fname_out) fclose(out);
}

static double random_zero_to_one() {
  // rand() returns an integer in the range [0, RAND_MAX]
  // Casting to double ensures floating-point division
  return (double)rand() / ((double) RAND_MAX + 1.0);
}

static void rowop_binstring(cfile_t cf, char *fname_out, config_rowop_t *cfg) {
  cdata_t c = read_cdata1(&cf);
  if (c.n == 0) return;    // nothing in cfile
  uint64_t n = cdata_n(&c);
  uint64_t binstring_bytes = 0;
  uint8_t *binstring = NULL;
  uint64_t k=0;
  for (k=0; ; ++k) {
    if (k) c = read_cdata1(&cf); // skip 1st cdata
    if (c.n == 0) break;
    cdata_t c2 = decompress(c);

    if (binstring_bytes*8 <= k) {
      binstring_bytes++;
      binstring = realloc(binstring, (binstring_bytes*n));
      memset(binstring + (binstring_bytes-1)*n, 0, sizeof(n));
    }
    
    switch (c.fmt) {
    case '3': {
      for (uint64_t i=0; i<c2.n; ++i) {
        uint64_t mu = f3_get_mu(&c2, i);
        /* if ((mu>>32) > (mu<<32>>32)) { */
        if (mu) {
          
          if (MU2beta(mu) > cfg->beta_threshold) {
            binstring[(k>>3)*n+i] |= (1<<(k&0x7));
          } else if (MU2beta(mu) == cfg->beta_threshold) {
            if (random_zero_to_one()>0.5)
              binstring[(k>>3)*n+i] |= (1<<(k&0x7));
          }
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

void rowop_cometh(cfile_t cf, char *fname_out, config_rowop_t *cfg) {

  uint64_t *cnts = NULL; uint64_t ncnts = 0;
  int cometh_window = cfg->cometh_window;
  for (uint64_t k=0; ;++k) {
    cdata_t c0 = read_cdata1(&cf);
    if (c0.n == 0) break;
    cdata_t c = decompress(c0);
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
        if (M+U >= cfg->mincov && M1+U1 >= cfg->mincov) {
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
      if (cfg->verbose) {
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

  int c;
  config_rowop_t config = {
    .beta0 = 0.4,
    .beta1 = 0.6,
    .mincov = 1,
    .beta_threshold = 0.5,
    .cometh_window = 5,
    .seed = (unsigned) time(NULL),
    .verbose = 0};
    
  char *op = NULL;
  while ((c = getopt(argc, argv, "vo:p:q:c:b:w:s:h"))>=0) {
    switch (c) {
    case 'o': op = strdup(optarg); break;
    case 'p': config.beta0 = atof(optarg); break;
    case 'q': config.beta1 = atof(optarg); break;
    case 'c': config.mincov = atoi(optarg); break;
    case 'b': config.beta_threshold = atof(optarg); break;
    case 'w': config.cometh_window = atoi(optarg); break;
    case 's': config.seed = atoi(optarg); break;
    case 'v': config.verbose = 1; break;
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
    cout = rowop_binasum(cf, &config);
    cdata_write(fname_out, &cout, "wb", config.verbose);
    free(cout.s);
  } else if (strcmp(op, "stat") == 0) {
    rowop_stat(cf, fname_out, &config);
  } else if (strcmp(op, "musum") == 0) {
    cout = rowop_musum(cf);
    cdata_write(fname_out, &cout, "wb", config.verbose);
    free(cout.s);
  } else if (strcmp(op, "binstring") == 0) {
    rowop_binstring(cf, fname_out, &config);
  } else if (strcmp(op, "cometh") == 0) {
    rowop_cometh(cf, fname_out, &config);
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


