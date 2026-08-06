#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rowstat.h"

/* Fixed point, scale 2^31. Betas are accumulated as integers so a row's totals
 * do not depend on the order samples arrived -- which is what makes the
 * multi-threaded merge exact rather than merely close. */
#define STAT_FX_BITS 31
#define STAT_FX_ONE  (1u << STAT_FX_BITS)
#define STAT_FX_SCALE ((double) STAT_FX_ONE)

#define NB YAME_ROWSTAT_NBINS

void yame_rowstat_init(yame_rowstat_t *a, uint64_t n) {
  a->n      = n;
  a->hist   = calloc(n * NB, sizeof(uint16_t));
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

void yame_rowstat_free(yame_rowstat_t *a) {
  free(a->cnts); free(a->sum); free(a->sum_sq); free(a->b0max);
  free(a->b1min); free(a->b0sum); free(a->b1sum); free(a->b0n); free(a->b1n);
  free(a->hist);
  memset(a, 0, sizeof(*a));
}

/* The <0.5 / >0.5 split still tests the double, so which side a beta lands on
 * is decided exactly as before; only the accumulation is fixed point. */
void yame_rowstat_add_fmt3(yame_rowstat_t *a, cdata_t *c, uint32_t mincov) {
  for (uint64_t i = 0; i < c->n; ++i) {
    uint64_t mu0 = f3_get_mu(c, i);
    if (!mu0) continue;                       /* 0-0 is skipped */
    if (((mu0>>32) + (mu0<<32>>32)) < mincov) continue;
    uint64_t M = mu0>>32;
    uint64_t U = (mu0<<32>>32);
    if (M+U >= mincov) {
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
      if (x != 0.5) {                         /* exactly 0.5 joins neither group */
        int bin = (int)(x * NB);
        if (bin >= NB) bin = NB - 1;
        uint16_t *h = a->hist + i * NB;
        if (h[bin] < UINT16_MAX) h[bin]++;    /* saturate rather than wrap */
      }
      if (x < 0.5) {
        a->b0n[i]++;
        a->b0sum[i] += xf;
        if (xf > a->b0max[i]) a->b0max[i] = xf;
      }
      if (x > 0.5) {
        a->b1n[i]++;
        a->b1sum[i] += xf;
        if (xf < a->b1min[i]) a->b1min[i] = xf;
      }
    }
  }
}

void yame_rowstat_merge(yame_rowstat_t *d, const yame_rowstat_t *s,
                        uint64_t beg, uint64_t end) {
  for (uint64_t i = beg; i < end; ++i) {
    d->cnts[i]   += s->cnts[i];
    d->sum[i]    += s->sum[i];
    d->sum_sq[i] += s->sum_sq[i];
    d->b0sum[i]  += s->b0sum[i];
    d->b1sum[i]  += s->b1sum[i];
    d->b0n[i]    += s->b0n[i];
    d->b1n[i]    += s->b1n[i];
    if (s->b0max[i] > d->b0max[i]) d->b0max[i] = s->b0max[i];
    if (s->b1min[i] < d->b1min[i]) d->b1min[i] = s->b1min[i];
    uint16_t *dh = d->hist + i*NB;
    const uint16_t *sh = s->hist + i*NB;
    for (int b = 0; b < NB; ++b) {
      unsigned v = (unsigned)dh[b] + sh[b];
      dh[b] = v > UINT16_MAX ? UINT16_MAX : (uint16_t)v;
    }
  }
}

int yame_rowstat_get(const yame_rowstat_t *a, uint64_t i,
                     yame_rowstat_row_t *o) {
  memset(o, 0, sizeof(*o));
  if (!a->cnts[i]) return 0;
  o->count = a->cnts[i];

  /* fixed point back to a fraction only here, once per row, after all the
   * combining is done -- so the arithmetic that had to be exact was */
  o->mean = (double) a->sum[i] / STAT_FX_SCALE / a->cnts[i];
  /* a true variance of zero (one sample, or every sample equal) can come out
   * a hair negative from the rounding; clamp rather than emit NaN */
  double var = ((double) a->sum_sq[i] / STAT_FX_SCALE / a->cnts[i])
             - o->mean * o->mean;
  o->sd = sqrt(var > 0 ? var : 0.0);
  o->min_n = (uint32_t)(a->b1n[i] < a->b0n[i] ? a->b1n[i] : a->b0n[i]);

  o->has_groups = (a->b0n[i] > 0 && a->b1n[i] > 0);
  if (!o->has_groups) return 1;

  o->delta_beta = ((double) a->b1min[i] - (double) a->b0max[i]) / STAT_FX_SCALE;
  o->delta_mean = ((double) a->b1sum[i] / a->b1n[i]
                 - (double) a->b0sum[i] / a->b0n[i]) / STAT_FX_SCALE;

  /* delta_q = q05(beta > 0.5) - q95(beta < 0.5): the same worst-case idea as
   * delta_beta but tolerant of a few outliers on each side. delta_beta uses
   * min and max, so ONE sample sitting near the middle collapses it; delta_mean
   * uses group centres and cannot see a straggler at all. The quantile gap sits
   * between them, which is where a useful separation filter belongs. Read off
   * the histogram, so the value is a bin edge (0.0625 granularity) -- and a
   * threshold falling between edges therefore behaves as the nearer edge. */
  const uint16_t *h = a->hist + i * NB;
  uint32_t lo_n = 0, hi_n = 0;
  for (int b = 0; b < NB/2; ++b) lo_n += h[b];
  for (int b = NB/2; b < NB; ++b) hi_n += h[b];
  if (!lo_n || !hi_n) return 1;               /* has_quantiles stays 0 */

  /* q95 of the low group: first bin whose cumulative count reaches 95% */
  double need = 0.95 * lo_n; uint32_t c = 0; int qb = NB/2 - 1;
  for (int b = 0; b < NB/2; ++b) {
    c += h[b];
    if ((double)c >= need) { qb = b; break; }
  }
  /* q05 of the high group: same from the bottom of the upper half */
  double need1 = 0.05 * hi_n; uint32_t c1 = 0; int qb1 = NB/2;
  for (int b = NB/2; b < NB; ++b) {
    c1 += h[b];
    if ((double)c1 >= need1) { qb1 = b; break; }
  }
  /* upper edge of the low bin against the lower edge of the high bin, so the
   * gap is never optimistic about how close the two groups came */
  o->has_quantiles = 1;
  o->q95_0 = (double)(qb + 1) / NB;
  o->q05_1 = (double)qb1 / NB;
  o->delta_q = o->q05_1 - o->q95_0;
  return 1;
}
