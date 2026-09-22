/* summary_index.h -- an INVERTED INDEX over the mask side of `yame summary`.
 *
 * SPDX-License-Identifier: LicenseRef-CHOP-Academic-BSD-2-Clause
 * Copyright (C) 2021-present The Children's Hospital of Philadelphia; see LICENSE.
 *
 * THE TWO SHAPES. The walk (summary_multi.h) is keyed by MASK: for each mask it
 * reads 64 rows at a time and popcounts against the query's coverage. Its cost
 * is the total mask content, the number of (row, mask) memberships, and it
 * pays that cost again for every query record.
 *
 * The index is keyed by ROW: for each row it lists the slots that claim it.
 * Its cost per record is one visit per row plus one add per membership at a
 * covered row. The mask count does not enter, and the build is paid once for
 * the whole run however many records follow.
 *
 * WHAT IT COSTS, measured on TFBS (1,359 masks, 29.4M rows, 390M memberships)
 * against whole-genome format 3 methylomes, machine idle:
 *
 *     records      walk      index    index RSS
 *        1        6.1 s     23.7 s      1.0 GB
 *        3       17.6 s     23.0 s      1.0 GB
 *        9       48.2 s     23.6 s      1.0 GB
 *
 * The walk costs 5.5 s per record and 45 MB. The index costs a 22 s build plus
 * 0.2 s per record and 4 bytes per row plus 2 per membership. Break-even was 4
 * records here and 2.4 on the machine methscope measured, because the build is
 * bound by decompressing the mask file and the walk by memory traffic. So YAME
 * does NOT pick between them: the walk is the default, `summary -I` asks for
 * the index, and a library caller calls the build it wants. The caller knows
 * its record count and its machine; a number baked in here would be wrong
 * somewhere, and the wrong index costs 1 GB where the wrong walk costs time.
 *
 * SCOPE. Format 3 queries. Binary (format 0/1) masks from a file; a state mask
 * claims every row, so its index would hold rows x masks entries and the build
 * declines. From RUNS, any slot layout the caller chooses, states included,
 * because an MRMP set claims a few rows however many states it has. */
#ifndef _YAME_SUMMARY_INDEX_H
#define _YAME_SUMMARY_INDEX_H

#include <stdint.h>
#include <stddef.h>
#include "cdata.h"
#include "summary_multi.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  uint64_t  n_rows;
  uint32_t  n_slots;   /* accumulators: masks, or masks x states             */
  uint32_t *off;       /* n_rows + 1, into ent                               */
  uint16_t *ent;       /* off[n_rows] slot ids, grouped by row               */
  uint64_t *n_m;       /* n_slots: rows each slot claims                     */
  uint64_t  bytes;     /* what off, ent and n_m hold, for a report           */
  double    t_build;   /* seconds                                            */
} yame_index_t;

/* The bytes an index over n_rows rows and n_memberships memberships holds,
 * before any build. A caller with the machine in view decides here; it is the
 * number it would otherwise recompute from the struct layout. The build holds
 * 4 more bytes per row while it runs, for the counts. */
uint64_t yame_index_bytes(uint64_t n_rows, uint64_t n_memberships);

/* The budget `summary -I` uses when YAME_SUMMARY_INDEX_MB is unset: 2 GB. */
#define YAME_INDEX_BUDGET_DEFAULT (2048ull << 20)

/* Read every record of `mask_path` twice (count, then fill) and invert it.
 *
 * NULL, with a sentence in `why` when given, if: the file has no record; a
 * record is not a binary mask; a record's rows disagree with n_rows; more than
 * 65,535 records (an entry is 16 bits); or the index would exceed
 * `budget_bytes` (0 means no budget). The offsets alone, 4 bytes per row, are
 * checked before anything is read. The memberships are checked as the count
 * pass runs; once the total crosses the line the pass stops scattering and
 * only popcounts the remaining records, so `why` still says the exact size
 * that would have fit. */
yame_index_t *yame_index_build(const char *mask_path, uint64_t n_rows,
                               uint64_t budget_bytes, char *why, size_t why_len);

/* The same from RUNS, the form an MRMP bank holds: `runs` is called twice,
 * once to count and once to fill. Slots as yame_summarize_multi() takes them:
 * with one state per mask pass n_slots == n_masks and emit state 0; with
 * states pass the flat total and `state_base`, the first slot of each mask.
 * n_slots must be at most 65,535.
 *
 * ONE OBLIGATION on the caller, which the kernel in summary_multi.h shares:
 * runs of ONE slot must not overlap. A row emitted twice for one slot is
 * indexed twice and counted twice. The kernel's other obligation, that runs of
 * one slot arrive in ascending order, is NOT needed here: the index visits
 * rows in ascending order whatever order the runs arrived in, so the sums are
 * bit-identical to the walk either way. A caller that satisfies the kernel
 * satisfies this. A run reaching past n_rows is clipped, and a run with a
 * mask or slot outside the layout fails the build.
 *
 * Both calls must emit the same runs; the second is trusted to match the
 * first, and a second call that emits MORE fails the build rather than
 * writing past the entries. */
yame_index_t *yame_index_build_runs(yame_runs_fn runs, void *runs_ctx,
                                    uint64_t n_rows, uint32_t n_masks,
                                    const uint32_t *state_base, uint32_t n_slots,
                                    uint64_t budget_bytes,
                                    char *why, size_t why_len);

void yame_index_free(yame_index_t *ix);

/* One query record against every slot. `acc` holds ix->n_slots entries and is
 * OVERWRITTEN. Rows are visited in ascending order, so each slot's sum is
 * added in the order the walk uses, and the result is bit-identical to it.
 * Returns 0, or -1 for a query this does not cover: not format 3, compressed,
 * or a row count other than the index's. */
int yame_index_apply(const yame_index_t *ix, const cdata_t *query,
                     yame_acc_t *acc);

/* The same, driven by the query's coverage bitmap (yame_qbits_build, or one
 * the caller filled from its own covered-row list). Only covered rows are
 * visited, in ascending order, so the numbers are the ones yame_index_apply()
 * gives and a single cell that covers a few percent of the rows costs a few
 * percent of the loop. yame_index_apply() itself reads every row once to find
 * the covered ones, which is what a bitmap already knows. Format 3 only; -1
 * otherwise, or for a row count other than the index's. */
int yame_index_apply_qb(const yame_index_t *ix, const yame_qbits_t *qb,
                        yame_acc_t *acc);

/* The same, driven by the caller's own covered-row LIST: rows[0..n) are the
 * rows to visit, each read for its M/U value as yame_index_apply() reads it.
 * For a caller that walked its record once and kept the covered rows -- an
 * MRMP featurizer draws its coverage ladder from such a list -- this is the
 * scatter with nothing in between: no second scan over the rows, and no
 * bitmap packed from the list only to be unpacked back into rows. A row with
 * no coverage in the record is skipped, not counted, so a list that
 * overstates still sums the truth.
 *
 * The list must be ASCENDING and DISTINCT, and that IS checked: a repeat
 * would count a row twice and a wrong order would add the doubles in another
 * sequence, so the sums would silently stop being bit-identical to the walk.
 * -1 for either, for a row at or past the index, or for a query this does
 * not cover (not format 3, compressed, other rows). */
int yame_index_apply_rows(const yame_index_t *ix, const cdata_t *query,
                          const uint32_t *rows, uint64_t n, yame_acc_t *acc);

#ifdef __cplusplus
}
#endif
#endif /* _YAME_SUMMARY_INDEX_H */
