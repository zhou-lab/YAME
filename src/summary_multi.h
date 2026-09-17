/* summary_multi.h -- summary accumulation that reads the query ONCE.
 *
 * SPDX-License-Identifier: AGPL-3.0-or-later
 *
 * WHAT THE OLD PATH COSTS. summarize1() takes one mask per call, so N masks
 * walk the query N times. Measured on the TFBS knowledgebase against a
 * 29.4M-row methylome, machine idle: 147 ms per set, of which about half is the
 * query walk. One query pass alone is 300 ms.
 *
 * WHAT IS SAVED, AND HOW. Two separate things:
 *
 *   1. The query is read once. A bitmap of its covered rows is built in one
 *      pass and reused for every mask. That is the saving the request was
 *      about, and it is about half the per-mask cost.
 *
 *   2. The mask is read 64 rows at a time. The first version of this kernel
 *      tested one bit per row, which left the other half of the cost in place
 *      and capped the gain at about 2x. Word-wise, the claimed count is a
 *      popcount of the mask words, the overlap count is a popcount of the mask
 *      AND the query bitmap, and only the surviving bits are turned back into
 *      row numbers. 460 thousand words per mask instead of 29 million rows.
 *
 * ONE MASK AT A TIME. Because the query bitmap is what carries state between
 * masks, nothing has to hold N masks in memory: `summary` streams them as it
 * always did and still reads the query once. The first version needed every
 * mask resident, which put a 5 GB floor under the full 1359-set TFBS set and
 * confined the speedup to -M.
 *
 * THE MASK SIDE CAN BE AN ENUMERATOR. A caller that holds membership as runs --
 * methscope walks them inside a model file and never materialises a .cm --
 * feeds runs rather than records. Requiring cdata_t would have cost it a
 * conversion it measured at 0.86 s per set, minutes on a 178-set bank.
 *
 * BIT-IDENTICAL, not approximately equal. Rows are still visited in increasing
 * order within a mask, so the same doubles are added in the same sequence.
 */
#ifndef _YAME_SUMMARY_MULTI_H
#define _YAME_SUMMARY_MULTI_H

#include <stdint.h>
#include "cdata.h"

#ifdef __cplusplus
extern "C" {
#endif

/* One mask's totals over one query record.
 *
 * n_u and n_q are here rather than returned once per query, because a format 6
 * query makes both depend on the MASK: its universe is the query's universe
 * intersected with the mask's, and every count is restricted to that. For a
 * format 3 query they are the same in every row of the output.
 *
 * `beta` is filled by the kernel rather than derived by the caller, because the
 * two query formats define it differently: a mean of per-row betas for format
 * 3, and the overlap over the mask count for format 6. */
typedef struct {
  uint64_t n_u;        /* the universe these counts are taken over        */
  uint64_t n_q;        /* rows the query claims, within that universe     */
  uint64_t n_m;        /* rows the mask claims, within that universe      */
  uint64_t n_o;        /* rows both claim                                 */
  uint64_t sum_depth;  /* summed coverage; 0 for a format 6 query         */
  double   sum_beta;   /* summed beta; 0 for a format 6 query             */
  double   beta;       /* as that query format defines it                 */
} yame_acc_t;

/* The query, read once, as bitmaps. Built per query record and reused for every
 * mask.
 *
 * For format 3 only `cov` is filled: bit i says row i has coverage. For format
 * 6 `cov` holds the universe and `qset` the set, because both are needed and
 * neither can be recovered from the other. */
typedef struct {
  const cdata_t *q;
  uint64_t *cov;       /* fmt3: covered.  fmt6: in the universe */
  uint64_t *qset;      /* fmt6 only: in the set */
  uint64_t  n_words;
  uint64_t  n_q;       /* fmt3: covered rows. fmt6: unused, it is per-mask */
} yame_qbits_t;

/* 0 on success, -1 for a query format this does not cover.
 *
 * Formats 3 and 6. They are different measurements, not one measurement over
 * two layouts: format 6 takes its universe from the query (and from the mask,
 * when the mask has one), restricts every count to it, and derives beta as
 * n_o/n_m rather than as a mean of betas. Both are implemented; neither
 * borrows the other's definitions, which is what made it worth keeping them
 * apart rather than adding a branch.
 *
 * A format 6 query is only handled in the DEFAULT view. The 2-bit and meth
 * views count different things again, and fall back. */
int  yame_qbits_build(const cdata_t *query, yame_qbits_t *out);
void yame_qbits_free(yame_qbits_t *qb);

/* One mask against those bitmaps, word-wise. Binary (fmt 0/1) and set+universe
 * (fmt 6) masks; -1 for anything else, or for a length that disagrees with the
 * query, so the caller falls back and the old path reports it in its own
 * words. `acc` is OVERWRITTEN, not added to. */
int yame_summarize_one_cx(const yame_qbits_t *qb, cdata_t *mask, yame_acc_t *acc);

/* --------------------------------------------- the enumerator entry point -- */

/* Rows [start, start+len) of mask `mask` are members, in state `state`. A
 * binary mask emits state 0. */
typedef void (*yame_emit_fn)(void *emit_ctx, uint32_t mask,
                             uint64_t start, uint64_t len, uint32_t state);

/* Called once; must emit every run of every mask. Order does not matter. */
typedef int (*yame_runs_fn)(void *ctx, yame_emit_fn emit, void *emit_ctx);

/*
 * N masks from runs, in one pass over the query.
 *
 * `acc` is indexed by the slot the caller chose when emitting: with one state
 * per mask pass n_acc == n_masks and emit state 0; with states pass the flat
 * total and give `state_base`, the first slot of each mask.
 *
 * `n_q_out` receives the covered-row count, computed once.
 */
int yame_summarize_multi(const cdata_t *query,
                         yame_runs_fn runs, void *runs_ctx,
                         uint32_t n_masks, const uint32_t *state_base,
                         yame_acc_t *acc, uint32_t n_acc,
                         uint64_t *n_q_out);

/* The same, fed from mask RECORDS. Kept for callers that already hold an array
 * of masks; `summary` uses yame_summarize_one_cx() in a loop instead, so it
 * never has to hold them all. */
int yame_summarize_multi_cx(const cdata_t *query, cdata_t *masks,
                            uint32_t n_masks, yame_acc_t *acc,
                            uint64_t *n_q_out);

#ifdef __cplusplus
}
#endif
#endif /* _YAME_SUMMARY_MULTI_H */
