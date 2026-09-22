// SPDX-License-Identifier: LicenseRef-CHOP-Academic-BSD-2-Clause
/**
 * This file is part of YAME.
 *
 * Copyright (C) 2021-present The Children's Hospital of Philadelphia
 *
 * Use of this software is available to academic and non-profit institutions
 * for research purposes under the 2-Clause BSD License; for use or transfers
 * to commercial entities, inquire with Dr. Wanding Zhou at zhouw3@chop.edu.
 * See the LICENSE file at the root of the repository for the full terms.
 */

#ifndef _SUMMARY_H
#define _SUMMARY_H

#include <stdint.h>
#include "kstring.h"
#include "cdata.h"

typedef struct stats_t {
  uint64_t sum_depth;           // sum of depth
  double sum_beta;
  double beta;
  uint64_t n_u;                 // universe
  uint64_t n_q;                 // query
  uint64_t n_m;                 // mask
  uint64_t n_o;                 // overlap
  char* sm;                     // mask name
  char* sq;                     // query name
} stats_t;

/**
 * How the two bits of a format-6 record are to be read. The same encoding is
 * produced with different meanings: `pack -f b` / feature sets mean
 * universe=tested background and set=member, while `binarize` means
 * universe=covered and set=methylated. The counts are the same either way; the
 * view decides how the summary columns are named and derived.
 */
typedef enum {
  F6_VIEW_SET = 0,              // universe = background, set = feature member
  F6_VIEW_METH,                 // universe = covered,    set = methylated
  F6_VIEW_2BIT,                 // four quaternary states counted separately
} f6_view_t;

typedef struct config_t {
  int full_name;
  int section_name;
  int in_memory;
  int use_index;              // -I: invert the mask file, see summary_index.h
  int no_header;
  f6_view_t f6_view;            // how to read format 6 (-V)
  char *fname_mask;
  char *fname_snames;
  char *fname_qry_stdin;
} config_t;

/**
 * Normalize a query/mask record into a representation summarize1() can consume:
 * fmt 0/1 are converted to an fmt0 bitset; fmt >= 2 are decompressed in place.
 */
void prepare_mask(cdata_t *c);

/**
 * Summarize one (already prepared) query record against one mask record.
 * Returns a heap array of n_st stats_t (caller frees each .sm/.sq and the array).
 * When c_mask is empty (n == 0), a whole-record summary is produced.
 */
stats_t* summarize1(cdata_t *c, cdata_t *c_mask, uint64_t *n_st,
                    char *sm, char *sq, config_t *config);

#endif /* _SUMMARY_H */
