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
#include "cfile.h"
#include "snames.h"

/**
 * yame rowop
 * ==========
 *
 * Overview
 * --------
 * rowop applies row-wise aggregation across multiple records (samples) in a CX file.
 * The input file is read sequentially record-by-record. Most operations require:
 *   - consistent format across records
 *   - consistent row dimension N across records
 *
 * Operations (-o; default = "binasum")
 * -----------------------------------
 *
 * 1) binasum  (CX output; fmt3)
 *    Purpose:
 *      Convert per-sample values into per-row sample counts:
 *        M = #samples called methylated
 *        U = #samples called unmethylated
 *      (these are counts of samples, not read-depth counts).
 *
 *    Input formats:
 *      - fmt0: bitset. If bit=1 -> M++, else U++.
 *      - fmt1: ASCII '0'/'1'. If nonzero -> M++, else U++.
 *      - fmt3: MU counts. For each sample and row:
 *          * skip if mu==0 (M=U=0) or cov < mincov
 *          * compute beta = M/(M+U)
 *          * if beta > beta1 => M++
 *            else if beta < beta0 => U++
 *            else ignored (ambiguous / intermediate methylation)
 *
 *    Options used:
 *      -c mincov, -p beta0, -q beta1
 *
 * 2) musum  (CX output; fmt3)
 *    Purpose:
 *      Sum MU sequencing counts across samples.
 *    Input:
 *      fmt3 only.
 *    Behavior:
 *      For each row, add sample M and U into output (skip mu==0).
 *
 * 3) stat  (text output; fmt3)
 *    Purpose:
 *      Compute per-row statistics across samples using beta values.
 *    Input:
 *      fmt3 only.
 *    Filters:
 *      skip mu==0 and cov < mincov.
 *    Output columns (tab-delimited):
 *      count    mean_beta    sd_beta    b0max    b1min
 *    Definitions:
 *      - count: number of samples contributing to this row
 *      - mean_beta: average beta across contributing samples
 *      - sd_beta: standard deviation of beta across contributing samples
 *      - b0max: max(beta) among samples with beta < 0.5 (0 if none)
 *      - b1min: min(beta) among samples with beta > 0.5 (1 if none)
 *    Implementation:
 *      Uses sum(beta) and sum(beta^2) to compute sd via sqrt(E[x^2]-E[x]^2).
 *
 * 4) binstring  (text output; fmt3)
 *    Purpose:
 *      Emit a row-wise binary string across samples. Each character corresponds
 *      to a sample in file order: '1' if beta > beta_threshold else '0'.
 *    Input:
 *      fmt3 only.
 *    Notes:
 *      Current implementation does not apply mincov (only checks mu!=0).
 *
 * 5) cometh  (text output; fmt3)
 *    Purpose:
 *      Summarize co-methylation between each row i and its neighbors (i+1..i+W),
 *      aggregated across samples.
 *    Input:
 *      fmt3 only.
 *    Filters:
 *      - requires cov >= mincov for both sites
 *      - skips intermediate methylation near 0.5 (|beta-0.5| < 0.2) by default logic
 *    Output:
 *      First column: 1-based row index.
 *      Then W columns (neighbor offsets). Each cell packs four 16-bit counters
 *      (UU, UM, MU, MM) into a uint64_t. With -v, prints lanes as "UU-UM-MU-MM".
 *
 * I/O
 * ---
 * - <in.cx> is required.
 * - [out] is optional; if omitted, text ops write to stdout and CX ops write to stdout
 *   (via cdata_write with fname_out == NULL).
 */

typedef struct config_rowop_t {
  double beta0; // lower threshold
  double beta1; // higher threshold
  unsigned mincov;
  double beta_threshold;   // default to 0.5
  int cometh_window;
  int verbose;
} config_rowop_t;

static int usage(void) {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "  yame rowop [options] <in.cx> [out]\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Purpose:\n");
  fprintf(stderr, "  Row-wise operations across multiple records (samples) in a CX file.\n");
  fprintf(stderr, "  Some operations write a new .cx; others write tab-delimited text.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Operation:\n");
  fprintf(stderr, "  -o <op>      Operation name (default: binasum)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "CX-output operations:\n");
  fprintf(stderr, "  binasum      Summarize per-sample calls into per-row sample counts (M/U) as format 3.\n");
  fprintf(stderr, "              Input:\n");
  fprintf(stderr, "                fmt0 : bitset, 1->M, 0->U\n");
  fprintf(stderr, "                fmt1 : ASCII '0'/'1', nonzero->M, zero->U\n");
  fprintf(stderr, "                fmt3 : MU; uses beta thresholds -p/-q, skips ambiguous betas\n");
  fprintf(stderr, "              Output: one fmt3 record (M=#methylated samples, U=#unmethylated samples).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  musum        Sum MU sequencing counts across samples.\n");
  fprintf(stderr, "              Input: fmt3 only. Output: one fmt3 record.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Text-output operations:\n");
  fprintf(stderr, "  stat         Per-row statistics across samples.\n");
  fprintf(stderr, "              Input: fmt3 only.\n");
  fprintf(stderr, "              Output columns:\n");
  fprintf(stderr, "                count  mean_beta  sd_beta  b0max  b1min\n");
  fprintf(stderr, "              where b0max is max(beta<0.5) and b1min is min(beta>0.5).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  binstring    Row-wise 0/1 string across samples using a single beta threshold.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  cometh       Neighbor co-methylation summary within a window (packed 4-way counts).\n");
  fprintf(stderr, "              With -v, prints four 16-bit lanes as UU-UM-MU-MM per neighbor offset.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Common filters:\n");
  fprintf(stderr, "  -c <mincov>  Minimum coverage (M+U) for a sample/row to contribute (default: 1).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "binasum (fmt3 input) thresholds:\n");
  fprintf(stderr, "  -p <beta0>   Call unmethylated if beta < beta0 (default: 0.4)\n");
  fprintf(stderr, "  -q <beta1>   Call methylated   if beta > beta1 (default: 0.6)\n");
  fprintf(stderr, "              Betas in [beta0, beta1] are ignored for that sample/row.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "binstring threshold:\n");
  fprintf(stderr, "  -b <beta>    Call methylated if beta > threshold (default: 0.5)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "cometh options:\n");
  fprintf(stderr, "  -w <W>       Neighbor window size (default: 5)\n");
  fprintf(stderr, "  -v           Verbose output (cometh prints UU-UM-MU-MM instead of packed uint64)\n");
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

static void collect_stat_fmt3(uint32_t *cnts, double *sum, double *sum_sq, double *b0max, double *b1min, cdata_t *c, config_rowop_t *cfg) {
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
      if (x < 0.5 && x > b0max[i]) b0max[i] = x;
      if (x > 0.5 && x < b1min[i]) b1min[i] = x;
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
  for (uint64_t i = 0; i < n; ++i) b1min[i] = 1.0;
  
  for (uint64_t k = 0; ; ++k) {
    if (k) c = read_cdata1(&cf); // skip 1st cdata
    if (c.n == 0) break;
    cdata_t c2 = decompress(c);

    switch (c.fmt) {
    case '3': collect_stat_fmt3(cnts, sum, sum_sq, b0max, b1min, &c2, cfg); break;
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
  fputs("count\tmean_beta\tsd_beta\tb0max\tb1min\n", out);
  for (uint64_t i=0; i<n; ++i) {
    if (!cnts[i]) {
      fputs("0\tNA\tNA\tNA\tNA\n", out);
    } else {
      double mean = sum[i] / cnts[i];
      double sd = sqrt((sum_sq[i] / cnts[i]) - mean * mean);
      fprintf(out, "%u\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n", cnts[i], mean, sd, b0max[i], b1min[i]);
    }
  }
  free(cnts);
  free(sum_sq);
  free(sum);
  if (fname_out) fclose(out);
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
        if (mu && ((double) (mu>>32) / ((mu>>32)+(mu<<32>>32)) > cfg->beta_threshold)) {
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
    .verbose = 0};
    
  char *op = NULL;
  while ((c = getopt(argc, argv, "vo:p:q:c:b:w:h"))>=0) {
    switch (c) {
    case 'o': op = strdup(optarg); break;
    case 'p': config.beta0 = atof(optarg); break;
    case 'q': config.beta1 = atof(optarg); break;
    case 'c': config.mincov = atoi(optarg); break;
    case 'b': config.beta_threshold = atof(optarg); break;
    case 'w': config.cometh_window = atoi(optarg); break;
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


