/* summary_multi.c -- N masks, one pass. See summary_multi.h.
 *
 * SPDX-License-Identifier: AGPL-3.0-or-later
 *
 * The shape is an interval sweep, not a per-row membership index. An index
 * from row to memberships is what a sparse bank wants -- 96.5% of CpGs belong
 * to nothing -- but it grows with the rows a DENSE mask claims, and
 * `yame summary -m` is routinely given masks that cover most of the genome.
 * The sweep costs one array of runs and an active list whose size is the
 * number of masks covering the current row, which is small in the sparse case
 * and correct in the dense one.
 */
#include <stdlib.h>
#include <string.h>
#include "summary_multi.h"
#include "wzmisc.h"

typedef struct {
  uint64_t start, end;      /* [start, end) */
  uint32_t slot;            /* the accumulator this run feeds */
} run_t;

typedef struct {
  run_t   *v;
  size_t   n, cap;
  uint32_t n_masks;
  const uint32_t *state_base;
  uint32_t n_acc;
  int      bad;             /* a run named a slot outside the accumulators */
} collect_t;

static void collect_emit(void *ctx, uint32_t mask,
                         uint64_t start, uint64_t len, uint32_t state) {
  collect_t *k = ctx;
  if (!len) return;
  if (mask >= k->n_masks) { k->bad = 1; return; }
  uint32_t slot = k->state_base ? k->state_base[mask] + state : mask;
  if (slot >= k->n_acc) { k->bad = 1; return; }
  if (k->n == k->cap) {
    k->cap = k->cap ? k->cap * 2 : 256;
    k->v = wzrealloc(k->v, k->cap * sizeof(run_t));
  }
  k->v[k->n].start = start;
  k->v[k->n].end   = start + len;
  k->v[k->n].slot  = slot;
  k->n++;
}

static int by_start(const void *a, const void *b) {
  const run_t *x = a, *y = b;
  if (x->start < y->start) return -1;
  if (x->start > y->start) return  1;
  return 0;
}

int yame_summarize_multi(const cdata_t *query,
                         yame_runs_fn runs, void *runs_ctx,
                         uint32_t n_masks, const uint32_t *state_base,
                         yame_acc_t *acc, uint32_t n_acc,
                         uint64_t *n_q_out) {
  if (!query || !runs || !acc || !n_acc) return -1;
  /* Format 3 only. A format 6 QUERY accumulates differently -- n_u counts the
   * query's own universe, n_m counts mask members RESTRICTED to that universe,
   * and beta is n_o/n_m rather than a mean of betas -- so the same sweep would
   * quietly produce different numbers. It is a separate kernel, not a branch
   * of this one, and until it is written and checked against
   * summarize1_queryfmt6() the caller falls back. */
  if (query->fmt != '3') return -1;

  collect_t k = {0};
  k.n_masks = n_masks; k.state_base = state_base; k.n_acc = n_acc;
  if (runs(runs_ctx, collect_emit, &k) != 0 || k.bad) { free(k.v); return -1; }

  /* n_m does not depend on the query at all: it is the rows a mask claims.
   * Counted here so the walk below only has to ask about coverage. A run that
   * runs past the query is clipped, the same as a per-mask call would see. */
  for (size_t i = 0; i < k.n; ++i) {
    if (k.v[i].start >= query->n) { k.v[i].end = k.v[i].start; continue; }
    if (k.v[i].end > query->n) k.v[i].end = query->n;
    acc[k.v[i].slot].n_m += k.v[i].end - k.v[i].start;
  }

  qsort(k.v, k.n, sizeof(run_t), by_start);

  /* The sweep. `act` holds the runs covering the current row; a run joins when
   * the row reaches its start and leaves when the row reaches its end. */
  size_t   *act = NULL, n_act = 0, cap_act = 0;
  size_t    next = 0;
  uint64_t  n_q = 0;

  /* `min_end` is what keeps this linear. Scanning the active list for expired
   * runs on EVERY row costs O(active) per row, which on 29.4M rows and 128
   * overlapping sets made the kernel SLOWER than the per-mask loop it replaces
   * (21.8 s against 20.0 s, measured on TFBS). Runs only expire at their own
   * end, so the list is swept when the row reaches the earliest of them and
   * not before. */
  uint64_t min_end = UINT64_MAX;

  for (uint64_t i = 0; i < query->n; ++i) {
    while (next < k.n && k.v[next].start <= i) {
      if (k.v[next].end > i) {
        if (n_act == cap_act) {
          cap_act = cap_act ? cap_act * 2 : 16;
          act = wzrealloc(act, cap_act * sizeof(size_t));
        }
        act[n_act++] = next;
        if (k.v[next].end < min_end) min_end = k.v[next].end;
      }
      next++;
    }
    if (i >= min_end) {
      min_end = UINT64_MAX;
      for (size_t a = 0; a < n_act; ) {
        if (k.v[act[a]].end <= i) act[a] = act[--n_act];
        else { if (k.v[act[a]].end < min_end) min_end = k.v[act[a]].end; ++a; }
      }
    }

    /* Nothing is active and nothing starts before the next run: the rows in
     * between only have to advance the covered count, so skip the per-row work
     * for the 96.5% of a sparse bank that belongs to no set at all. */

    /* Nothing claims this row: the common case for a sparse bank, and the
     * reason a sweep is cheap there. The covered count still has to advance. */
    uint64_t mu = f3_get_mu((cdata_t *)query, i);
    if (mu) n_q++;
    if (!n_act || !mu) continue;
    uint64_t cov = MU2cov(mu);
    double   b   = MU2beta(mu);
    for (size_t a = 0; a < n_act; ++a) {
      yame_acc_t *t = &acc[k.v[act[a]].slot];
      t->n_o++; t->sum_depth += cov; t->sum_beta += b;
    }
  }

  free(act); free(k.v);
  if (n_q_out) *n_q_out = n_q;
  return 0;
}

/* ---------------------------------------------------- fed from cdata_t ---- */

/*
 * The convenience wrapper: membership runs read out of mask RECORDS, so
 * `yame summary -m multi.cm` uses the same kernel as a caller that walks its
 * own runs. Binary masks (fmt 0/1) and set+universe masks (fmt 6) only; a
 * state mask (fmt 2) has a per-term shape the caller must flatten itself, and
 * gets -1 here so it falls back.
 */
typedef struct { cdata_t *m; uint32_t n; } cx_ctx_t;

static int cx_runs(void *ctx, yame_emit_fn emit, void *emit_ctx) {
  cx_ctx_t *c = ctx;
  for (uint32_t j = 0; j < c->n; ++j) {
    cdata_t *m = &c->m[j];
    if (m->fmt != '0' && m->fmt != '1' && m->fmt != '6') return -1;
    /* Runs rather than single rows: a mask that claims a contiguous block
     * becomes one entry instead of thousands, which is what keeps the sweep's
     * array small for the dense masks `summary -m` is usually given. */
    uint64_t start = 0; int in = 0;
    for (uint64_t i = 0; i < m->n; ++i) {
      int here = (m->fmt == '6')
        ? (FMT6_IN_UNI(*m, i) && FMT6_IN_SET(*m, i))
        : (FMT0_IN_SET(*m, i) ? 1 : 0);
      if (here && !in) { start = i; in = 1; }
      else if (!here && in) { emit(emit_ctx, j, start, i - start, 0); in = 0; }
    }
    if (in) emit(emit_ctx, j, start, m->n - start, 0);
  }
  return 0;
}

int yame_summarize_multi_cx(const cdata_t *query, cdata_t *masks, uint32_t n_masks,
                            yame_acc_t *acc, uint64_t *n_q_out) {
  for (uint32_t j = 0; j < n_masks; ++j)
    if (masks[j].n != query->n) return -1;   /* length mismatch: let the
                                              * per-mask path report it */
  cx_ctx_t ctx = { masks, n_masks };
  return yame_summarize_multi(query, cx_runs, &ctx, n_masks, NULL,
                              acc, n_masks, n_q_out);
}
