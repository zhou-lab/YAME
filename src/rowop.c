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
#include "yame_ui.h"
#include <sys/types.h>
#include <time.h>
#include <pthread.h>
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
 *      delta_beta  min(beta>0.5) - max(beta<0.5)   (worst-case group margin)
 *      min_n       min(#beta<0.5, #beta>0.5)
 *      delta_mean  mean(beta>0.5) - mean(beta<0.5) (group-center separation)
 *
 *    Notes:
 *      - delta_beta/delta_mean are only defined when both sides exist; else NA.
 *      - delta_beta measures the gap between the closest members of the two
 *        groups (sensitive to group extremes); delta_mean measures the
 *        separation of group means (robust to a single near-0.5 member).
 *      - sd is computed as sqrt(E[x^2] - E[x]^2).
 *
 * 4) binstring  (text output; fmt3)
 *    Purpose:
 *      Emit a row-wise binary string across samples (one line per CpG row,
 *      one character per sample).
 *    Behavior:
 *      For each sample at a CpG:
 *        - confident '1' if beta >  beta_threshold (-b)
 *        - confident '0' if beta <  beta_threshold
 *        - "ambiguous"   if mu==0, cov < mincov (-c), or beta == beta_threshold
 *      Ambiguous cells are filled deterministically with the CpG's majority
 *      confident state (more 0s -> 0, more 1s -> 1; exact tie -> 0). This
 *      replaces the previous random (-s seeded) tie-break, so output is now
 *      reproducible and -s no longer affects binstring.
 *      The fill is only trusted when the majority is "sweeping": the larger
 *      confident side must be >= -M <fold> times the smaller (unanimous counts
 *      as sweeping; default fold = 10). If a CpG has ambiguous cells but the
 *      majority is below this fold, the whole line is emitted as the '2'
 *      sentinel (see Filter). CpGs with no ambiguous cells are unaffected.
 *    Filter:
 *      A CpG's line is emitted as an all-'2' sentinel (same length as other
 *      lines, preserving positional alignment) when its ambiguous fraction
 *      exceeds -m <frac> (default 1.0 = off), or when it has ambiguous cells
 *      but no sweeping majority to fill them. '2' never appears otherwise.
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
  double max_ambig_frac;   // binstring: max fraction of ambiguous cells per CpG (default 1.0 = no filtering)
  double min_major_fold;   // binstring: min fold (hi/lo) of confident calls to trust the majority fill (default 10)
  int cometh_window;
  int verbose;
  unsigned seed;
  int threads;             // -t; 1 = the serial scan
  int decimals;            // -d; digits printed for stat's fractions
} config_rowop_t;

static int usage(void) {
  yame_usage_head("yame rowop [options] <in.cx> [out]");
  yame_usage_sec("Purpose:");
  yame_usage_text("Perform row-wise operations across multiple records (samples) in a CX file.");
  yame_usage_text("Depending on the operation, output is either a new CX file or plain text.");
  yame_usage_sec("Operation:");
  yame_usage_opt("-o <op>", "Operation name (default: binasum)");
  yame_usage_sec("CX-output operations:");
  yame_usage_text("binasum      Convert per-sample values into per-row sample counts (M/U) as format 3.");
  yame_usage_cont("Input: fmt0, fmt1, or fmt3.");
  yame_usage_cont("For fmt3, beta thresholds (-p/-q) define methylated vs unmethylated calls.");
  yame_usage_text("musum        Sum MU sequencing counts across samples.");
  yame_usage_cont("Input: fmt3 only. Output: one fmt3 record.");
  yame_usage_sec("Text-output operations:");
  yame_usage_text("stat         Per-row summary statistics across samples.");
  yame_usage_cont("Input: fmt3 only.");
  yame_usage_cont("Output columns:");
  yame_usage_cont(" count  mean_beta  sd_beta  delta_beta  min_n  delta_mean  q95_0  q05_1  delta_q");
  yame_usage_cont("q95_0 = 95th pct of beta<0.5; q05_1 = 5th pct of beta>0.5.");
  yame_usage_cont("delta_q = q05_1 - q95_0: delta_beta's worst-case idea, but tolerant of");
  yame_usage_cont("an outlier per side. Both quantiles are reported, not just the gap,");
  yame_usage_cont("because a filter usually constrains ONE side -- \"no expected-0 sample");
  yame_usage_cont("creeps toward 0.5\" is q95_0, which the difference cannot express.");
  yame_usage_cont("delta_beta = min(beta>0.5) - max(beta<0.5)   (worst-case margin).");
  yame_usage_cont("min_n      = min(#beta<0.5, #beta>0.5).");
  yame_usage_cont("delta_mean = mean(beta>0.5) - mean(beta<0.5) (group-center separation).");
  yame_usage_text("binstring    Convert per-sample beta values into row-wise binary strings.");
  yame_usage_cont("Input: fmt3 only. Uses -b as the beta threshold, -c as min coverage.");
  yame_usage_cont("Ambiguous cells (mu==0, cov<mincov, or beta==threshold) are filled");
  yame_usage_cont("with the CpG's majority state (deterministic; -s does not apply).");
  yame_usage_text("cometh       Neighbor co-methylation summary within a window.");
  yame_usage_cont("Input: fmt3 only.");
  yame_usage_cont("Output: packed 4-way counts (UU, UM, MU, MM) per neighbor offset.");
  yame_usage_cont("Use -v to print unpacked lanes.");
  yame_usage_sec("Common filters:");
  yame_usage_opt("-c <mincov>", "Minimum coverage (M+U) for a sample/row to contribute (default: 1).");
  yame_usage_opt("-d <0-9>", "Decimals printed for stat's fractions (default 6). The old");
  yame_usage_cont("default of 3 put many values on a rounding boundary: a beta is a");
  yame_usage_cont("ratio of small integers, so a mean like 9/400 = 0.0225 sits exactly");
  yame_usage_cont("between 0.022 and 0.023 and any change to the arithmetic flips it.");
  yame_usage_cont("It also matters downstream, since a threshold applied to a printed");
  yame_usage_cont("value inherits half the last digit as slop.");
  yame_usage_sec("Threads:");
  yame_usage_opt("-t <N>", "Split the records over N threads (default 1). Applies to");
  yame_usage_cont("binasum, musum and stat, whose accumulators combine; the other");
  yame_usage_cont("ops stay serial. An indexed file lets each thread seek to its");
  yame_usage_cont("own run, so the decompression parallelises too; a stream or an");
  yame_usage_cont("unindexed file is read once and dispatched to the threads, which");
  yame_usage_cont("still parallelises everything after the inflate. Each thread");
  yame_usage_cont("holds its own accumulator, so memory is N x (8 bytes/row) for");
  yame_usage_cont("binasum/musum and N x (84 bytes/row) for stat -- at hg38 scale,");
  yame_usage_cont("235 MB and 2.2 GB per thread (stat also carries a 16-bin beta");
  yame_usage_cont("histogram per row for delta_q). Output is byte-identical at every");
  yame_usage_cont("N for all three ops.");
  yame_usage_text("binasum (fmt3 input) thresholds:");
  yame_usage_opt("-p <beta0>", "Call unmethylated if beta < beta0 (default: 0.4).");
  yame_usage_opt("-q <beta1>", "Call methylated   if beta > beta1 (default: 0.6).");
  yame_usage_cont("Betas in [beta0, beta1] are ignored.");
  yame_usage_text("binstring options:");
  yame_usage_opt("-b <beta>", "Call methylated if beta > threshold (default: 0.5).");
  yame_usage_opt("-m <frac>", "Max ambiguous fraction per CpG; above this the line is emitted");
  yame_usage_cont("as an all-'2' sentinel (default: 1.0 = no filtering).");
  yame_usage_opt("-M <fold>", "Min majority fold (hi/lo of confident calls) to trust the fill");
  yame_usage_cont("(default: 10). If a CpG has ambiguous cells but the majority is");
  yame_usage_cont("below this fold, its line is emitted as the all-'2' sentinel.");
  yame_usage_text("cometh options:");
  yame_usage_opt("-w <W>", "Neighbor window size (default: 5).");
  yame_usage_opt("-v", "Verbose output (print UU-UM-MU-MM instead of packed uint64).");
  yame_usage_sec("Other:");
  yame_usage_opt("-h", "Show this help message.");
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

/* The sufficient statistics stat() carries per row.
 *
 * The sums are fixed point, not floating point: a beta in [0,1] is carried as
 * round(beta * 2^31), so accumulation is integer addition. That matters
 * because integer addition is associative and floating-point addition is not
 * -- with doubles, splitting the records across threads regroups the sums and
 * can move the last decimal, and delta_mean feeds threshold filters (the MRMP
 * >=0.6 cut) where that changes which probes are selected. In fixed point
 * every thread count gives bit-identical output by construction.
 *
 * Widths: beta_fx <= 2^31 fits uint32; a sum of N of them fits int64 for N up
 * to 2^32 samples. sum_sq keeps beta^2 * 2^31 by shifting the product back
 * down, so it fits the same width. Per-term rounding is 2^-32 ~ 2e-10, seven
 * orders below the three decimals printed.
 *
 * Every field combines -- sums add, the extremes take min/max -- which is what
 * lets the record list be split at all. */
#define STAT_FX_BITS 31
#define STAT_FX_ONE  (1u << STAT_FX_BITS)          /* 1.0 in fixed point */
#define STAT_FX_SCALE ((double) STAT_FX_ONE)

typedef struct {
  uint32_t *cnts;
  int64_t *sum, *sum_sq, *b0sum, *b1sum;   /* fixed point, scale 2^31 */
  uint32_t *b0max, *b1min;                 /* fixed point, scale 2^31 */
  int *b0n, *b1n;
  uint16_t *hist;                          /* n * STAT_NBINS beta histogram */
} statacc_t;

/* Order statistics need more than running sums, so carry a small per-row
 * histogram of beta. 16 bins put the boundary exactly at 0.5 (bins 0-7 below,
 * 8-15 above), matching the group split, and 0.0625 resolution is finer than
 * the data: a pooled single-cell reference has beta = M/(M+U) over ~11 cells,
 * so its true granularity is ~1/11. Counts add, so this merges across threads
 * exactly like the sums and preserves the byte-identical guarantee. */
#define STAT_NBINS 16

static void statacc_alloc(statacc_t *a, uint64_t n) {
  a->hist   = calloc(n * STAT_NBINS, sizeof(uint16_t));
  a->cnts   = calloc(n, sizeof(uint32_t));
  a->sum    = calloc(n, sizeof(int64_t));
  a->sum_sq = calloc(n, sizeof(int64_t));
  a->b0max  = calloc(n, sizeof(uint32_t));
  a->b1min  = calloc(n, sizeof(uint32_t));
  a->b0sum  = calloc(n, sizeof(int64_t));
  a->b1sum  = calloc(n, sizeof(int64_t));
  a->b0n    = calloc(n, sizeof(int));
  a->b1n    = calloc(n, sizeof(int));
  for (uint64_t i = 0; i < n; ++i) a->b1min[i] = STAT_FX_ONE;
}

static void statacc_free(statacc_t *a) {
  free(a->cnts); free(a->sum); free(a->sum_sq); free(a->b0max);
  free(a->b1min); free(a->b0sum); free(a->b1sum); free(a->b0n); free(a->b1n);
  free(a->hist);
}


/* The <0.5 / >0.5 split still tests the double, so which side a beta lands on
 * is decided exactly as before; only the accumulation moved to fixed point. */
static void collect_stat_fmt3(statacc_t *a, cdata_t *c, config_rowop_t *cfg) {
  for (uint64_t i=0; i<c->n; ++i) {
    uint64_t mu0 = f3_get_mu(c, i);
    if (!mu0) continue; // 0-0 is skipped
    if (((mu0>>32) + (mu0<<32>>32)) < cfg->mincov) continue;
    uint64_t M = mu0>>32;
    uint64_t U = (mu0<<32>>32);
    if (M+U >= cfg->mincov) {
      double x = (double) M / (M+U);
      uint32_t xf = (uint32_t) llround(x * STAT_FX_SCALE);
      /* round the square back down rather than truncate: truncation biases
       * sum_sq low, which drives the variance below zero on rows where every
       * sample agrees and turns sd into NaN */
      uint64_t sq = (uint64_t) xf * xf;
      int64_t xf2 = (int64_t)((sq + (1ull << (STAT_FX_BITS-1))) >> STAT_FX_BITS);
      a->sum[i] += xf;
      a->sum_sq[i] += xf2;
      a->cnts[i]++;
      if (x != 0.5) {                      /* exactly 0.5 joins neither group */
        int bin = (int)(x * STAT_NBINS);
        if (bin >= STAT_NBINS) bin = STAT_NBINS - 1;
        uint16_t *h = a->hist + i * STAT_NBINS;
        if (h[bin] < UINT16_MAX) h[bin]++;  /* saturate rather than wrap */
      }
      if (x < 0.5) {
        a->b0n[i]++;
        a->b0sum[i] += xf;
        if (xf > a->b0max[i])
          a->b0max[i] = xf;
      }
      if (x > 0.5) {
        a->b1n[i]++;
        a->b1sum[i] += xf;
        if (xf < a->b1min[i])
          a->b1min[i] = xf;
      }
    }
  }
}


static void stat_emit(statacc_t *a, uint64_t n, char *fname_out, int dec);

// the following standard deviation doesn't work for large numbers but should be ok for meth levels
// see https://www.strchr.com/standard_deviation_in_one_pass
static void rowop_stat(cfile_t cf, char *fname_out, config_rowop_t *cfg) {

  cdata_t c = read_cdata1(&cf);
  if (c.n == 0) return; // nothing in cfile, output nothing
  uint64_t n = cdata_n(&c);
  statacc_t a; statacc_alloc(&a, n);

  for (uint64_t k = 0; ; ++k) {
    if (k) c = read_cdata1(&cf); // skip 1st cdata
    if (c.n == 0) break;
    cdata_t c2 = decompress(c);

    switch (c.fmt) {
    case '3': collect_stat_fmt3(&a, &c2, cfg); break;
    default: {
      fprintf(stderr, "[%s:%d] File format: %c unsupported.\n", __func__, __LINE__, c.fmt);
      fflush(stderr);
      exit(1);
    }}

    free(c.s); free(c2.s);
  }

  stat_emit(&a, n, fname_out, cfg->decimals);
  statacc_free(&a);
}

static void stat_emit(statacc_t *a, uint64_t n, char *fname_out, int dec) {

  uint32_t *cnts = a->cnts;
  int64_t *sum = a->sum, *sum_sq = a->sum_sq, *b0sum = a->b0sum, *b1sum = a->b1sum;
  uint32_t *b0max = a->b0max, *b1min = a->b1min;
  int *b0n = a->b0n, *b1n = a->b1n;

  FILE *out;
  if (fname_out) {
    out = fopen(fname_out, "w");
  } else {
    out = stdout;
  }

  fputs("count\tmean_beta\tsd_beta\tdelta_beta\tmin_n\tdelta_mean"
        "\tq95_0\tq05_1\tdelta_q\n", out);
  for (uint64_t i = 0; i < n; ++i) {
    if (cnts[i] == 0) {
      fputs("0\tNA\tNA\tNA\t0\tNA\tNA\tNA\tNA\n", out);
      continue;
    }

    /* fixed point back to a fraction only here, once per row, after all the
     * combining is done -- so the arithmetic that had to be exact was */
    double mean = (double) sum[i] / STAT_FX_SCALE / cnts[i];
    /* a true variance of zero (one sample, or every sample equal) can come out
     * a hair negative from the rounding; clamp rather than emit NaN */
    double var  = ((double) sum_sq[i] / STAT_FX_SCALE / cnts[i]) - mean * mean;
    double sd   = sqrt(var > 0 ? var : 0.0);

    /* delta_beta = b1min - b0max, but only meaningful if both sides exist */
    double delta_beta = (b0n[i] > 0 && b1n[i] > 0)
      ? ((double) b1min[i] - (double) b0max[i]) / STAT_FX_SCALE : -1.0;

    /* delta_mean = mean(beta>0.5) - mean(beta<0.5); robust group separation,
       only defined when both sides exist */
    int both = (b0n[i] > 0 && b1n[i] > 0);
    double delta_mean = both
      ? ((double) b1sum[i] / b1n[i] - (double) b0sum[i] / b0n[i]) / STAT_FX_SCALE
      : -1.0;

    /* min_n = min(#beta<0.5, #beta>0.5) */
    uint32_t min_n = (b1n[i] < b0n[i]) ? b1n[i] : b0n[i];

    if (delta_beta < 0) {
      fprintf(out, "%u\t%1.*f\t%1.*f\tNA\t%u\t",
              cnts[i], dec, mean, dec, sd, min_n);
    } else {
      fprintf(out, "%u\t%1.*f\t%1.*f\t%1.*f\t%u\t",
              cnts[i], dec, mean, dec, sd, dec, delta_beta, min_n);
    }
    if (both) fprintf(out, "%1.*f\t", dec, delta_mean);
    else      fputs("NA\t", out);

    /* delta_q = q05(beta > 0.5) - q95(beta < 0.5): the same worst-case idea as
     * delta_beta but tolerant of a few outliers on each side. delta_beta uses
     * min and max, so ONE class sitting near the middle collapses it; delta_mean
     * uses group centres and cannot see a straggler at all. The quantile gap
     * sits between them, which is where a useful separation filter belongs.
     * Read off the histogram, so the reported value is a bin edge (0.0625
     * granularity) -- finer than a pooled single-cell reference resolves. */
    if (both) {
      const uint16_t *h = a->hist + i * STAT_NBINS;
      uint32_t lo_n = 0, hi_n = 0;
      for (int b = 0; b < STAT_NBINS/2; ++b) lo_n += h[b];
      for (int b = STAT_NBINS/2; b < STAT_NBINS; ++b) hi_n += h[b];
      if (!lo_n || !hi_n) { fputs("NA\tNA\tNA\n", out); continue; }
      /* q95 of the low group: first bin whose cumulative count reaches 95% */
      double need = 0.95 * lo_n; uint32_t c = 0; int qb = STAT_NBINS/2 - 1;
      for (int b = 0; b < STAT_NBINS/2; ++b) {
        c += h[b];
        if ((double)c >= need) { qb = b; break; }
      }
      /* q05 of the high group: same from the bottom of the upper half */
      double need1 = 0.05 * hi_n; uint32_t c1 = 0; int qb1 = STAT_NBINS/2;
      for (int b = STAT_NBINS/2; b < STAT_NBINS; ++b) {
        c1 += h[b];
        if ((double)c1 >= need1) { qb1 = b; break; }
      }
      /* upper edge of the low bin against the lower edge of the high bin, so
       * the gap is never optimistic about how close the two groups came */
      double q95_lo = (double)(qb + 1) / STAT_NBINS;
      double q05_hi = (double)qb1 / STAT_NBINS;
      /* Both quantiles are reported, not only their gap: a filter usually wants
       * to constrain ONE side -- "no expected-0 class creeps toward 0.5" is
       * q95_0, and says something the difference cannot, since a large gap can
       * hide a high q95_0 sitting under an even higher q05_1. */
      fprintf(out, "%1.*f\t%1.*f\t%1.*f\n", dec, q95_lo, dec, q05_hi,
              dec, q05_hi - q95_lo);
    } else fputs("NA\tNA\tNA\n", out);
  }
  if (fname_out) fclose(out);   /* the accumulator belongs to the caller */
}

/**
 * Write one sample's column of the bit planes from its decompressed record.
 *
 * The planes are laid out per group of 8 samples -- sample k lives in plane
 * k>>3 at bit k&7, and each plane is n bytes, one per row.
 */
static void binstring_fill(uint8_t *bs, uint8_t *amb, uint64_t k,
                           uint64_t n, cdata_t *c2, config_rowop_t *cfg) {
  for (uint64_t i = 0; i < c2->n; ++i) {
    uint64_t mu = f3_get_mu(c2, i);
    if (!mu || MU2cov(mu) < cfg->mincov) {            // no / low depth
      amb[(k>>3)*n+i] |= (1<<(k&0x7));
    } else {
      double beta = MU2beta(mu);
      if (beta > cfg->beta_threshold) {               // confident methylated
        bs[(k>>3)*n+i] |= (1<<(k&0x7));
      } else if (beta == cfg->beta_threshold) {       // M==U tie
        amb[(k>>3)*n+i] |= (1<<(k&0x7));
      }                                               // else confident 0
    }
  }
}

/**
 * Emit one line per row: the per-row majority, then a character per sample.
 *
 * Split out from the fill because the two phases cost very differently -- the
 * fill is about 30% of a binstring run and this is the rest. That is why
 * threading the fill alone is pointless (measured 0.99 s serial against
 * 1.01 s on 8 threads), and why anyone returning to speed binstring up should
 * start here. Both costs scale as rows x samples, so the ratio does not
 * improve with a bigger input.
 */
static void binstring_emit(uint8_t *bs, uint8_t *amb, uint64_t ncells,
                           uint64_t n, config_rowop_t *cfg, char *fname_out) {
  FILE *out;
  if (fname_out) { out = fopen(fname_out, "w"); }
  else { out = stdout; }

  for (uint64_t i=0; i<n; ++i) {
    // Per-CpG majority over confidently-called cells; ambiguous cells filled with it.
    uint64_t n1 = 0, namb = 0;
    for (uint64_t kk=0; kk<ncells; ++kk) {
      int is_amb = (amb[(kk>>3)*n+i] >> (kk&0x7)) & 0x1;
      if (is_amb) { namb++; }
      else if ((bs[(kk>>3)*n+i] >> (kk&0x7)) & 0x1) { n1++; }
    }
    uint64_t n0 = ncells - n1 - namb;
    char fill = (n1 > n0) ? '1' : '0';  // exact tie -> 0

    // Only trust the majority fill when it is "sweeping" (larger side >= fold x smaller;
    // unanimous counts as sweeping). Filling from a near-even split fabricates calls.
    uint64_t hi = (n1 > n0) ? n1 : n0;
    uint64_t lo = (n1 > n0) ? n0 : n1;
    int sweeping = (hi > 0) && (lo == 0 || (double)hi >= cfg->min_major_fold * (double)lo);

    if ((ncells && (double)namb > cfg->max_ambig_frac * (double)ncells) ||
        (namb > 0 && !sweeping)) {
      for (uint64_t kk=0; kk<ncells; ++kk) fputc('2', out);  // filtered sentinel
    } else {
      for (uint64_t kk=0; kk<ncells; ++kk) {
        if ((amb[(kk>>3)*n+i] >> (kk&0x7)) & 0x1) {
          fputc(fill, out);
        } else {
          fputc('0'+((bs[(kk>>3)*n+i] >> (kk&0x7))&0x1), out);
        }
      }
    }
    fputc('\n', out);
  }
  if (fname_out) fclose(out);
}

static void rowop_binstring(cfile_t cf, char *fname_out, config_rowop_t *cfg) {
  cdata_t c = read_cdata1(&cf);
  if (c.n == 0) return;    // nothing in cfile
  uint64_t n = cdata_n(&c);
  uint64_t binstring_bytes = 0;
  uint8_t *binstring = NULL;   // confident methylated (1) bits
  uint8_t *ambig = NULL;       // ambiguous cells (tie / low-or-no depth)
  uint64_t k=0;
  for (k=0; ; ++k) {
    if (k) c = read_cdata1(&cf); // skip 1st cdata
    if (c.n == 0) break;
    cdata_t c2 = decompress(c);

    if (binstring_bytes*8 <= k) {
      binstring_bytes++;
      binstring = realloc(binstring, (binstring_bytes*n));
      ambig = realloc(ambig, (binstring_bytes*n));
      memset(binstring + (binstring_bytes-1)*n, 0, n);
      memset(ambig + (binstring_bytes-1)*n, 0, n);
    }

    switch (c.fmt) {
    case '3':
      binstring_fill(binstring, ambig, k, n, &c2, cfg);
      break;
    default: {
      fprintf(stderr, "[%s:%d] File format: %c unsupported.\n", __func__, __LINE__, c.fmt);
      fflush(stderr);
      exit(1);
    }}

    free(c.s); free(c2.s);
  }

  uint64_t ncells = k;
  binstring_emit(binstring, ambig, ncells, n, cfg, fname_out);
  free(binstring);
  free(ambig);
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

/* ---- parallel reduction over records ------------------------------------
 *
 * Every reducing op walks the records one at a time into an accumulator that
 * is indexed by genomic row, so the accumulator's size is O(rows) and never
 * O(samples x rows), and combining two of them is just the op's merge. That
 * makes the record list splittable: give each thread a contiguous run of
 * records and its own accumulator, then merge in chunk order.
 *
 * Threads need to seek to their first record, so this wants the .cg index.
 * Without one -- or reading a stream -- there is nothing to partition and the
 * serial scan runs instead, with a word on stderr rather than silently.
 *
 * Exactness: binasum and musum accumulate integer counts, so any thread count
 * gives byte-identical output. stat sums doubles, and floating-point addition
 * is not associative, so a different -t can move the last bits (well below the
 * 3 decimals printed). Merging in chunk order at least makes a given -t
 * reproducible.
 */

typedef enum { ROWOP_BINASUM, ROWOP_MUSUM, ROWOP_STAT } rowop_kind_t;

/* Wall-clock seconds, for the -v phase breakdown. Threading only helps the
 * reduce phase, so knowing how much time is spent outside it is what tells
 * you whether more threads will buy anything. */
static double rowop_now(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double) ts.tv_sec + (double) ts.tv_nsec * 1e-9;
}

typedef struct {
  char *fname;
  int64_t *off;                 /* record virtual offsets, file order */
  rowop_kind_t kind;
  config_rowop_t *cfg;
  uint64_t n;                   /* rows */
  char fmt;                     /* format of record 0 */
  cdata_t acc;                  /* binasum/musum */
  statacc_t st;                 /* stat */
  struct recq_t *q;             /* non-NULL: take records from the queue
                                 * instead of seeking (stream input) */
  struct rowdisp_t *d;          /* shared cursor over the record list */
  double secs;                  /* -v: this worker's own reduce time */
  uint64_t nrec, bytes;         /* -v: records and compressed bytes handled */
  int err;                      /* non-zero on failure */
  char errmsg[256];
} rowop_worker_t;

/**
 * Which record a thread does next, decided when it asks rather than up front.
 *
 * Splitting the list into equal *counts* is not an equal split of *work*:
 * records differ in coverage, so they differ in compressed size, and time
 * tracks size closely. Measured over 200 single cells on 8 threads, an equal
 * count of 25 records each gave 279.6 to 389.8 MB per thread and 11.48 to
 * 14.29 s, leaving 13% of the reduce phase idle behind the straggler.
 *
 * Weighting the split by bytes would still be predicting. Handing out one
 * record at a time measures instead: a thread that draws a heavy record only
 * delays itself, and the worst imbalance falls to a single record rather than
 * a whole share. It is also what the stream path does already, which is why
 * that path balanced itself.
 */
typedef struct rowdisp_t {
  int next, n;
  pthread_mutex_t mu;
} rowdisp_t;

/* -1 when the list is exhausted */
static int rowdisp_take(rowdisp_t *d) {
  pthread_mutex_lock(&d->mu);
  int k = (d->next < d->n) ? d->next++ : -1;
  pthread_mutex_unlock(&d->mu);
  return k;
}

/**
 * A bounded queue of *compressed* records, for input that cannot be seeked.
 *
 * With an index each thread seeks to its own run, so the BGZF inflate happens
 * in parallel. A pipe has no index and no way to find record boundaries
 * without inflating, so one reader has to walk the stream and hand records
 * out. Only the inflate is serialised that way; decompressing the CX payload
 * and accumulating still spread across the workers.
 *
 * The queue holds records still compressed. A decompressed fmt3 record is 8
 * bytes a row -- 235 MB at hg38 scale -- so queueing those would cost more
 * memory than the whole reduction.
 */
typedef struct recq_t {
  cdata_t *slot;
  int cap, head, tail, count, done;
  pthread_mutex_t mu;
  pthread_cond_t not_empty, not_full;
} recq_t;

static void recq_init(recq_t *q, int cap) {
  q->slot = calloc(cap, sizeof(cdata_t));
  q->cap = cap; q->head = q->tail = q->count = q->done = 0;
  pthread_mutex_init(&q->mu, NULL);
  pthread_cond_init(&q->not_empty, NULL);
  pthread_cond_init(&q->not_full, NULL);
}

static void recq_free(recq_t *q) {
  free(q->slot);
  pthread_mutex_destroy(&q->mu);
  pthread_cond_destroy(&q->not_empty);
  pthread_cond_destroy(&q->not_full);
}

static void recq_push(recq_t *q, cdata_t c) {
  pthread_mutex_lock(&q->mu);
  while (q->count == q->cap) pthread_cond_wait(&q->not_full, &q->mu);
  q->slot[q->tail] = c;
  q->tail = (q->tail + 1) % q->cap;
  q->count++;
  pthread_cond_signal(&q->not_empty);
  pthread_mutex_unlock(&q->mu);
}

/* 1 with a record in *out, 0 when the stream is finished and drained */
static int recq_pop(recq_t *q, cdata_t *out) {
  pthread_mutex_lock(&q->mu);
  while (q->count == 0 && !q->done) pthread_cond_wait(&q->not_empty, &q->mu);
  if (q->count == 0) { pthread_mutex_unlock(&q->mu); return 0; }
  *out = q->slot[q->head];
  q->head = (q->head + 1) % q->cap;
  q->count--;
  pthread_cond_signal(&q->not_full);
  pthread_mutex_unlock(&q->mu);
  return 1;
}

static void recq_close(recq_t *q) {
  pthread_mutex_lock(&q->mu);
  q->done = 1;
  pthread_cond_broadcast(&q->not_empty);
  pthread_mutex_unlock(&q->mu);
}


/* Merging serially costs (nt-1) passes over every row, so the tail grew with
 * the thread count -- measured 0.57 s at -t 2 and 3.96 s at -t 8, which is
 * the wrong direction. Rows are independent, so instead each thread folds all
 * the accumulators down for one slice of the rows. Same arithmetic, same
 * order (accumulator 0 upward), so the result is unchanged. */
typedef struct {
  rowop_worker_t *w;
  int nt;
  rowop_kind_t kind;
  uint64_t beg, end;            /* rows */
} rowop_merger_t;

static void *rowop_merge_rows(void *arg) {
  rowop_merger_t *m = (rowop_merger_t*) arg;
  for (int t = 1; t < m->nt; ++t) {
    if (m->kind == ROWOP_STAT) {
      statacc_t *d = &m->w[0].st, *s = &m->w[t].st;
      for (uint64_t i = m->beg; i < m->end; ++i) {
        d->cnts[i]   += s->cnts[i];
        d->sum[i]    += s->sum[i];
        d->sum_sq[i] += s->sum_sq[i];
        d->b0sum[i]  += s->b0sum[i];
        d->b1sum[i]  += s->b1sum[i];
        d->b0n[i]    += s->b0n[i];
        d->b1n[i]    += s->b1n[i];
        if (s->b0max[i] > d->b0max[i]) d->b0max[i] = s->b0max[i];
        if (s->b1min[i] < d->b1min[i]) d->b1min[i] = s->b1min[i];
        { uint16_t *dh = d->hist + i*STAT_NBINS, *sh = s->hist + i*STAT_NBINS;
          for (int b = 0; b < STAT_NBINS; ++b) {
            unsigned v = (unsigned)dh[b] + sh[b];
            dh[b] = v > UINT16_MAX ? UINT16_MAX : (uint16_t)v; } }
      }
    } else {
      cdata_t *d = &m->w[0].acc, *s = &m->w[t].acc;
      for (uint64_t i = m->beg; i < m->end; ++i) {
        uint64_t a = f3_get_mu(d, i), b = f3_get_mu(s, i);
        f3_set_mu(d, i, (a>>32) + (b>>32), (a<<32>>32) + (b<<32>>32));
      }
    }
  }
  return NULL;
}

/* Fold one already-read record into this worker's accumulator. */
static void rowop_absorb(rowop_worker_t *w, cdata_t c) {
    if (c.fmt != w->fmt) {
      w->err = 1;
      snprintf(w->errmsg, sizeof(w->errmsg),
               "formats are inconsistent: %c vs %c", w->fmt, c.fmt);
      return;
    }
    cdata_t c2 = decompress(c);
    if (c2.n != w->n) {
      w->err = 1;
      snprintf(w->errmsg, sizeof(w->errmsg),
               "dimensions are inconsistent: %"PRIu64" vs %"PRIu64, w->n, c2.n);
      free(c2.s);
      return;
    }
    switch (w->kind) {
    case ROWOP_BINASUM:
      switch (w->fmt) {
      case '0': binasumFmt0(&w->acc, &c2); break;
      case '1': binasumFmt1(&w->acc, &c2); break;
      case '3': binasumFmt3(&w->acc, &c2, w->cfg); break;
      default:
        w->err = 1;
        snprintf(w->errmsg, sizeof(w->errmsg), "format %c unsupported", w->fmt);
      }
      break;
    case ROWOP_MUSUM:
      if (w->fmt == '3') musumFmt3(&w->acc, &c2);
      else {
        w->err = 1;
        snprintf(w->errmsg, sizeof(w->errmsg), "format %c unsupported", w->fmt);
      }
      break;
    case ROWOP_STAT:
      if (w->fmt == '3')
        collect_stat_fmt3(&w->st, &c2, w->cfg);
      else {
        w->err = 1;
        snprintf(w->errmsg, sizeof(w->errmsg), "format %c unsupported", w->fmt);
      }
      break;
    }
    free(c2.s);
}

static void *rowop_worker(void *arg) {
  rowop_worker_t *w = (rowop_worker_t*) arg;
  double t0 = rowop_now();

  if (w->q) {                   /* stream: take whatever the reader hands out */
    cdata_t c;
    while (!w->err && recq_pop(w->q, &c)) {
      w->nrec++; w->bytes += cdata_nbytes(&c);
      rowop_absorb(w, c);
      free(c.s);
    }
    w->secs = rowop_now() - t0;
    return NULL;
  }

  cfile_t cf = open_cfile(w->fname);
  cdata_t c = {0};
  for (;;) {                    /* take the next record whenever free */
    int k = rowdisp_take(w->d);
    if (k < 0) break;
    if (bgzf_seek(cf.fh, w->off[k], SEEK_SET) < 0) {
      w->err = 1;
      snprintf(w->errmsg, sizeof(w->errmsg), "cannot seek record %d", k);
      break;
    }
    if (!read_cdata2(&cf, &c)) {
      w->err = 1;
      snprintf(w->errmsg, sizeof(w->errmsg), "short read at record %d", k);
      break;
    }
    w->nrec++; w->bytes += cdata_nbytes(&c);
    rowop_absorb(w, c);
    if (w->err) break;
  }
  free(c.s);
  bgzf_close(cf.fh);
  w->secs = rowop_now() - t0;
  return NULL;
}

/* Record offsets in file order, or NULL when the input cannot be partitioned
 * (no index). Caller frees. */
static int64_t *rowop_record_offsets(char *fname, int *n_rec) {
  *n_rec = 0;
  if (strcmp(fname, "-") == 0) return NULL;
  char *fname_index = get_fname_index(fname);
  index_t *idx = loadIndex(fname_index);
  free(fname_index);
  if (!idx) return NULL;
  int npairs = 0;
  index_pair_t *pairs = index_pairs(idx, &npairs);
  if (npairs <= 0) { clean_index_pairs(pairs, npairs); cleanIndex(idx); return NULL; }
  int64_t *off = malloc(npairs * sizeof(int64_t));
  for (int i = 0; i < npairs; ++i) off[i] = pairs[i].value;
  clean_index_pairs(pairs, npairs);
  cleanIndex(idx);
  /* index_pairs() gives hash order; the file order is by offset */
  for (int i = 1; i < npairs; ++i) {   /* insertion sort: n is sample count */
    int64_t v = off[i]; int j = i - 1;
    while (j >= 0 && off[j] > v) { off[j+1] = off[j]; --j; }
    off[j+1] = v;
  }
  *n_rec = npairs;
  return off;
}

/**
 * Run a reducing op across `threads` workers. Returns 1 when it ran, 0 when
 * the input could not be partitioned and the caller should fall back.
 * On success the reduced result is left in *out (binasum/musum) or *stout
 * (stat), which the caller emits exactly as the serial path does.
 */
static int rowop_parallel(char *fname, rowop_kind_t kind, config_rowop_t *cfg,
                          cdata_t *out, statacc_t *stout, uint64_t *n_out) {

  int n_rec = 0;
  int64_t *off = rowop_record_offsets(fname, &n_rec);
  /* No index: fall back to reading the stream once and handing records out. */
  int streaming = (off == NULL);
  if (!streaming && n_rec < 2) { free(off); return 0; }  /* nothing to split */

  /* rows and format come from the first record, as in the serial path. For a
   * stream that record cannot be re-read, so it is kept and queued first. */
  cfile_t cf0 = open_cfile(fname);
  cdata_t c0 = read_cdata1(&cf0);
  if (c0.n == 0) { free(c0.s); bgzf_close(cf0.fh); free(off); return 0; }
  char fmt = c0.fmt;
  uint64_t n = cdata_n(&c0);
  if (!streaming) { free(c0.s); bgzf_close(cf0.fh); }

  int nt = cfg->threads;
  if (!streaming && nt > n_rec) nt = n_rec;

  recq_t q;
  if (streaming) recq_init(&q, nt + 2);
  rowdisp_t disp = { .next = 0, .n = n_rec };
  pthread_mutex_init(&disp.mu, NULL);

  rowop_worker_t *w = calloc(nt, sizeof(rowop_worker_t));
  pthread_t *tid = calloc(nt, sizeof(pthread_t));
  for (int t = 0; t < nt; ++t) {
    w[t].fname = fname; w[t].off = off; w[t].kind = kind; w[t].cfg = cfg;
    w[t].n = n; w[t].fmt = fmt;
    w[t].q = streaming ? &q : NULL;
    w[t].d = streaming ? NULL : &disp;
    if (kind == ROWOP_STAT) statacc_alloc(&w[t].st, n);
    else {
      w[t].acc.n = n; w[t].acc.compressed = 0;
      w[t].acc.fmt = '3'; w[t].acc.unit = 8;
      w[t].acc.s = calloc(n, sizeof(uint64_t));
    }
  }
  if (cfg->verbose) {
    if (streaming)
      fprintf(stderr, "[rowop] no index: one reader feeding %d threads.\n", nt);
    else
      fprintf(stderr, "[rowop] %d records over %d threads, taken one at a "
              "time.\n", n_rec, nt);
  }

  double t_reduce = rowop_now();
  for (int t = 0; t < nt; ++t) pthread_create(&tid[t], NULL, rowop_worker, &w[t]);
  if (streaming) {
    /* this thread is the reader: inflate records and hand them out */
    recq_push(&q, c0);          /* record 0, already in hand */
    for (;;) {
      cdata_t c = {0};
      if (!read_cdata2(&cf0, &c)) { free(c.s); break; }
      recq_push(&q, c);
    }
    recq_close(&q);
  }
  for (int t = 0; t < nt; ++t) pthread_join(tid[t], NULL);
  if (streaming) { recq_free(&q); bgzf_close(cf0.fh); }
  t_reduce = rowop_now() - t_reduce;

  for (int t = 0; t < nt; ++t) {
    if (w[t].err) {
      fprintf(stderr, "[%s:%d] %s\n", __func__, __LINE__, w[t].errmsg);
      fflush(stderr);
      exit(1);
    }
  }

  /* merge into worker 0, split by rows so the tail does not grow with -t */
  double t_merge = rowop_now();
  if (nt > 1) {
    rowop_merger_t *m = calloc(nt, sizeof(rowop_merger_t));
    uint64_t rows_per = (n + nt - 1) / nt;
    for (int t = 0; t < nt; ++t) {
      m[t].w = w; m[t].nt = nt; m[t].kind = kind;
      m[t].beg = (uint64_t) t * rows_per;
      m[t].end = m[t].beg + rows_per < n ? m[t].beg + rows_per : n;
      if (m[t].beg > n) m[t].beg = n;
      pthread_create(&tid[t], NULL, rowop_merge_rows, &m[t]);
    }
    for (int t = 0; t < nt; ++t) pthread_join(tid[t], NULL);
    free(m);
    for (int t = 1; t < nt; ++t) {
      if (kind == ROWOP_STAT) statacc_free(&w[t].st);
      else free(w[t].acc.s);
    }
  }
  t_merge = rowop_now() - t_merge;

  if (kind == ROWOP_STAT) *stout = w[0].st;
  else *out = w[0].acc;
  *n_out = n;

  if (cfg->verbose) {
    fprintf(stderr, "[rowop] reduce %.2f s (%d threads), merge %.2f s.\n",
            t_reduce, nt, t_merge);
    /* Per-worker time against per-worker input. If the slowest thread is much
     * slower than the fastest, the split is uneven and the reduce phase ends
     * when the straggler does; if they all agree, the shortfall is elsewhere. */
    double lo = 1e30, hi = 0, tot = 0;
    for (int t = 0; t < nt; ++t) {
      if (w[t].secs < lo) lo = w[t].secs;
      if (w[t].secs > hi) hi = w[t].secs;
      tot += w[t].secs;
      fprintf(stderr, "[rowop]   thread %d: %6.2f s  %4"PRIu64" records  "
              "%6.1f MB\n", t, w[t].secs, w[t].nrec, w[t].bytes/1e6);
    }
    fprintf(stderr, "[rowop]   spread: fastest %.2f s, slowest %.2f s "
            "(%.0f%% idle waiting on the straggler)\n",
            lo, hi, nt > 0 ? 100.0 * (1.0 - tot / (hi * nt)) : 0.0);
  }

  pthread_mutex_destroy(&disp.mu);
  free(w); free(tid); free(off);
  return 1;
}

int main_rowop(int argc, char *argv[]) {

  int c;
  config_rowop_t config = {
    .beta0 = 0.4,
    .beta1 = 0.6,
    .mincov = 1,
    .beta_threshold = 0.5,
    .max_ambig_frac = 1.0,
    .min_major_fold = 10.0,
    .cometh_window = 5,
    .seed = (unsigned) time(NULL),
    .verbose = 0,
    .threads = 1,
    .decimals = 6};

  char *op = NULL;
  while ((c = getopt(argc, argv, "vo:p:q:c:b:w:s:m:M:t:d:h"))>=0) {
    switch (c) {
    case 'o': op = strdup(optarg); break;
    case 't': config.threads = atoi(optarg); break;
    case 'd': {
      char *endp = NULL;
      long v = strtol(optarg, &endp, 10);
      if (endp == optarg || *endp || v < 0 || v > 9)
        wzfatal("-d takes 0-9 decimals, given \"%s\".\n", optarg);
      config.decimals = (int) v;
      break;
    }
    case 'p': config.beta0 = atof(optarg); break;
    case 'q': config.beta1 = atof(optarg); break;
    case 'c': config.mincov = atoi(optarg); break;
    case 'b': config.beta_threshold = atof(optarg); break;
    case 'w': config.cometh_window = atoi(optarg); break;
    case 's': config.seed = atoi(optarg); break;
    case 'm': config.max_ambig_frac = atof(optarg); break;
    case 'M': config.min_major_fold = atof(optarg); break;
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

  if (config.threads < 1) wzfatal("-t takes a thread count >= 1, given %d.\n",
                                  config.threads);

  /* The reducing ops can split the record list; the rest are serial. */
  if (config.threads > 1) {
    rowop_kind_t kind; int reducing = 1;
    if (!op || strcmp(op, "binasum") == 0) kind = ROWOP_BINASUM;
    else if (strcmp(op, "musum") == 0)     kind = ROWOP_MUSUM;
    else if (strcmp(op, "stat") == 0)      kind = ROWOP_STAT;
    else { reducing = 0; kind = ROWOP_BINASUM; }

    if (!reducing) {
      fprintf(stderr, "[rowop] -t applies to binasum, musum and stat; %s runs "
              "single-threaded.\n", op);
    } else {
      cdata_t cout = {0}; statacc_t st = {0}; uint64_t n = 0;
      if (rowop_parallel(fname, kind, &config, &cout, &st, &n)) {
        double t_out = rowop_now();
        if (kind == ROWOP_STAT) { stat_emit(&st, n, fname_out, config.decimals); statacc_free(&st); }
        else { cdata_write(fname_out, &cout, "wb", config.verbose); free(cout.s); }
        /* The write is single-threaded whatever -t says: for a CX result it is
         * a BGZF deflate of the whole accumulator. When this is a large share
         * of the run, more threads cannot help. */
        if (config.verbose)
          fprintf(stderr, "[rowop] write %.2f s (single-threaded).\n",
                  rowop_now() - t_out);
        if (fname_out) free(fname_out);
        free(op);
        return 0;
      }
    }
  }

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


