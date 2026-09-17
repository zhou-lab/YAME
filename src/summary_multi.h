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

/* One mask's totals over one query record. The universe is the query length
 * and the covered count is the same for every mask, so neither appears here. */
typedef struct {
  uint64_t n_m;        /* rows this mask claims, whatever the query holds */
  uint64_t n_o;        /* rows this mask claims AND the query covers       */
  uint64_t sum_depth;  /* summed coverage over those rows                  */
  double   sum_beta;   /* summed beta over those rows                      */
} yame_acc_t;

/* The query, read once: which rows it covers, as a bitmap, plus the record
 * itself for the values. Built per query record and reused for every mask. */
typedef struct {
  const cdata_t *q;
  uint64_t *cov;       /* bit i set when row i is covered */
  uint64_t  n_words;
  uint64_t  n_q;       /* covered rows, counted during the build */
} yame_qbits_t;

/* 0 on success, -1 for a query format this does not cover.
 *
 * Format 3 only. A format 6 QUERY counts its own universe, restricts the mask
 * count to that universe, and derives beta as n_o/n_m rather than as a mean of
 * betas, so the same machinery would quietly give different numbers. It is a
 * separate kernel, not a branch of this one. */
int  yame_qbits_build(const cdata_t *query, yame_qbits_t *out);
void yame_qbits_free(yame_qbits_t *qb);

/* One mask against that bitmap, word-wise. Binary (fmt 0/1) and set+universe
 * (fmt 6) masks; -1 for anything else, or for a length that disagrees with the
 * query, so the caller falls back and the old path reports it in its own
 * words. `acc` is added to, not reset. */
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
