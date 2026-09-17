/* summary_multi.h -- accumulate N masks in ONE pass over a query record.
 *
 * SPDX-License-Identifier: AGPL-3.0-or-later
 *
 * summarize1() takes one mask per call, so scoring a query against N masks
 * walks the query N times. The walk decompresses; the accumulation is a few
 * adds. Measured on a 29.4M-row methylome of 9 samples: 1 mask 1.6 s, 4 masks
 * 4.8 s, 16 masks 16.8 s -- linear in the mask count, which for a 178-set bank
 * is about half a minute per cell.
 *
 * THE MASK SIDE IS AN ENUMERATOR, NOT A RECORD. A caller that already holds
 * membership as runs -- methscope walks them inside a model file and never
 * materialises a .cm -- feeds them straight in. Requiring cdata_t would have
 * cost that caller a conversion it measured at 0.86 s per set, minutes on a
 * bank, and it would have kept its own copy instead. yame_summarize_multi_cx()
 * is the wrapper that feeds the same kernel from records, which is what
 * `yame summary -m` uses.
 *
 * The result is bit-identical to summarize1(), not approximately equal: within
 * each mask the rows are still accumulated in increasing order, so the same
 * doubles are added in the same sequence. Only the interleaving across masks
 * is new.
 */
#ifndef _YAME_SUMMARY_MULTI_H
#define _YAME_SUMMARY_MULTI_H

#include <stdint.h>
#include "cdata.h"

#ifdef __cplusplus
extern "C" {
#endif

/* One mask's totals over one query record. n_u and n_q do not appear: the
 * universe is the query length and the covered count is the same for every
 * mask, so both are returned once by the call rather than N times. */
typedef struct {
  uint64_t n_m;        /* rows this mask claims, whatever the query holds */
  uint64_t n_o;        /* rows this mask claims AND the query covers       */
  uint64_t sum_depth;  /* summed coverage over those rows                  */
  double   sum_beta;   /* summed beta over those rows                      */
} yame_acc_t;

/* A membership run: rows [start, start+len) of mask `mask`, in state `state`.
 * A binary mask emits state 0 for every run; a state mask emits the term. */
typedef void (*yame_emit_fn)(void *emit_ctx, uint32_t mask,
                             uint64_t start, uint64_t len, uint32_t state);

/* Called once. It must emit every run of every mask; order does not matter,
 * the kernel sorts. */
typedef int (*yame_runs_fn)(void *ctx, yame_emit_fn emit, void *emit_ctx);

/*
 * Accumulate `n_acc` accumulators over one query record in a single pass.
 *
 * `acc` is indexed by the (mask, state) slot the caller chose when emitting:
 * a caller with one state per mask passes n_acc == n_masks and emits state 0;
 * a caller with states passes n_acc == sum of states and emits the flat slot.
 * The kernel does not interpret `state` beyond adding it to `mask`'s base,
 * which the caller gives in `state_base` (NULL means one slot per mask).
 *
 * Returns 0, or -1 when the query format is one the kernel does not cover.
 * Format 3 today: the M/U methylome, which is where the cost is and what both
 * callers measured. A format 6 QUERY counts a different universe and derives
 * beta differently, so it needs its own kernel rather than a branch of this
 * one; a caller that gets -1 falls back to summarize1() per mask.
 *
 * `n_q_out` receives the covered-row count, computed once.
 */
int yame_summarize_multi(const cdata_t *query,
                         yame_runs_fn runs, void *runs_ctx,
                         uint32_t n_masks, const uint32_t *state_base,
                         yame_acc_t *acc, uint32_t n_acc,
                         uint64_t *n_q_out);

/* The same kernel, fed from mask RECORDS rather than from runs. Binary and
 * set+universe masks; returns -1 for anything else, or for a mask whose length
 * disagrees with the query, so the caller falls back to the per-mask path and
 * that path reports the mismatch in its own words. */
int yame_summarize_multi_cx(const cdata_t *query, cdata_t *masks,
                            uint32_t n_masks, yame_acc_t *acc,
                            uint64_t *n_q_out);

#ifdef __cplusplus
}
#endif
#endif /* _YAME_SUMMARY_MULTI_H */
