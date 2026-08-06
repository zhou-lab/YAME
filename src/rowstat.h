/* Per-row summary statistics across the samples of a format-3 (M/U) store.
 *
 * This is the engine behind `yame rowop -o stat`, exported because it is not
 * only a reporting statistic: the CpG filter that picks reference positions for
 * an MRMP is a predicate over these same numbers, and methscope links libyame
 * to build MRMPs. Keeping one implementation here means the filter and the
 * printed table can never disagree, and a caller that wants the filter does not
 * have to shell out to `rowop` and parse a multi-hundred-megabyte TSV.
 *
 * Accumulation is fixed point (scale 2^31) so that a row's value does not
 * depend on the order samples were added, which is what lets several threads
 * accumulate disjointly and merge exactly. Every field of the accumulator is
 * either a sum, a count, or an extremum, so yame_rowstat_merge() is exact and
 * `rowop -t N` is byte-identical at every N.
 */
#ifndef _YAME_ROWSTAT_H
#define _YAME_ROWSTAT_H

#include <stdint.h>
#include "cdata.h"

/* Beta histogram bins per row. 16 puts a boundary exactly at 0.5 -- matching
 * the <0.5 / >0.5 group split the other fields use -- and 1/16 = 0.0625 is
 * finer than the data: a pooled single-cell reference beta is M/(M+U) over
 * ~11 cells, so its true granularity is ~1/11. */
#define YAME_ROWSTAT_NBINS 16

/* Open, because rowop embeds one per worker and merges them pairwise. Callers
 * should treat the fields as opaque and read values through yame_rowstat_get(),
 * which is where fixed point turns back into a fraction. */
typedef struct {
  uint64_t  n;                             /* rows */
  uint32_t *cnts;
  int64_t  *sum, *sum_sq, *b0sum, *b1sum;  /* fixed point, scale 2^31 */
  uint32_t *b0max, *b1min;                 /* fixed point, scale 2^31 */
  int      *b0n, *b1n;
  uint16_t *hist;                          /* n * YAME_ROWSTAT_NBINS */
} yame_rowstat_t;

/* One row, decoded. `has_groups` is false when every observed sample fell on
 * the same side of 0.5, which leaves delta_beta / delta_mean / the quantiles
 * undefined -- `rowop` prints NA for exactly those. `has_quantiles` can be
 * false even when has_groups is true, in the degenerate case where the
 * histogram holds no counts on one side (only reachable via saturation at
 * 65535 samples). */
typedef struct {
  uint32_t count;        /* samples observed at this row */
  double   mean, sd;
  uint32_t min_n;        /* min(#beta<0.5, #beta>0.5) */
  int      has_groups;
  double   delta_beta;   /* min(beta>0.5) - max(beta<0.5): worst-case margin */
  double   delta_mean;   /* mean(beta>0.5) - mean(beta<0.5): group centres */
  int      has_quantiles;
  double   q95_0;        /* 95th pct of the beta<0.5 group (a bin edge) */
  double   q05_1;        /* 5th pct of the beta>0.5 group (a bin edge) */
  double   delta_q;      /* q05_1 - q95_0 */
} yame_rowstat_row_t;

/* Allocate for `n` rows, zeroed and with b1min primed to 1.0. */
void yame_rowstat_init(yame_rowstat_t *a, uint64_t n);
void yame_rowstat_free(yame_rowstat_t *a);

/* Fold one decompressed format-3 record in. Samples with total coverage below
 * `mincov`, and the 0/0 no-data entry, are skipped. */
void yame_rowstat_add_fmt3(yame_rowstat_t *a, cdata_t *c, uint32_t mincov);

/* d[beg,end) += s[beg,end). Exact: every field is a sum, count or extremum. */
void yame_rowstat_merge(yame_rowstat_t *d, const yame_rowstat_t *s,
                        uint64_t beg, uint64_t end);

/* Decode row i. Returns 0 and zeroes *o when nothing was observed there. */
int yame_rowstat_get(const yame_rowstat_t *a, uint64_t i, yame_rowstat_row_t *o);

#endif /* _YAME_ROWSTAT_H */
