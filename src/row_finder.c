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

#include "cdata.h"

row_finder_t init_finder(cdata_t *cr) {
  row_finder_t fdr = {0};
  fdr.h = kh_init(str2int);

  row_reader_t rdr = {0};
  char *chrm = NULL;
  chromosome_t *chrmt = NULL;
  while (row_reader_next_loc(&rdr, cr)) {
    if (rdr.chrm != chrm) { // a new chromosome
      chrm = rdr.chrm;
      if (kh_get(str2int, fdr.h, chrm) == kh_end(fdr.h)) { // key doesn't exist
        int ret;
        khiter_t k = kh_put(str2int, fdr.h, chrm, &ret);
        if (ret >= 0) kh_value(fdr.h, k) = fdr.n;
        else {
          fprintf(stderr, "[%s:%d] Cannot insert hashmap\n", __func__, __LINE__);
          fflush(stderr);
          exit(1);
        }
        fdr.n++;
        fdr.chrms = realloc(fdr.chrms, sizeof(chromosome_t)*fdr.n);
        chrmt = fdr.chrms + fdr.n - 1;
        memset(chrmt, 0, sizeof(chromosome_t));
      } else {
        fprintf(stderr, "[%s:%d] Chromosome %s appeared twice in unsorted manner.\n", __func__, __LINE__, chrm);
        fflush(stderr);
        exit(1);
      }
    }
    while ((rdr.value>>17) >= chrmt->n) {
      chrmt->locs = realloc(chrmt->locs, (chrmt->n+1)*sizeof(uint64_t));
      chrmt->vals = realloc(chrmt->vals, (chrmt->n+1)*sizeof(uint64_t));
      chrmt->inds = realloc(chrmt->inds, (chrmt->n+1)*sizeof(uint64_t));
      chrmt->locs[chrmt->n] = rdr.loc;
      chrmt->vals[chrmt->n] = rdr.value;
      chrmt->inds[chrmt->n] = rdr.index;
      chrmt->n++;
    }
  }
  return fdr;
}

uint64_t row_finder_search(char *chrm, uint64_t beg1, row_finder_t *fdr, cdata_t *cr) {

  khiter_t k = kh_get(str2int, fdr->h, chrm);
  if (k == kh_end(fdr->h)) {
    fprintf(stderr, "[%s:%d] Chromosome %s not found.\n", __func__, __LINE__, chrm);
    fflush(stderr);
    exit(1);
  }
  chromosome_t chrmt = fdr->chrms[kh_value(fdr->h, k)];
  if ((beg1>>17) >= chrmt.n) {
    fprintf(stderr, "[%s:%d] Coordinate %"PRIu64" is too big (max: %"PRIu64")\n", __func__, __LINE__, beg1, chrmt.n);
    fflush(stderr);
    exit(1);
  }

  row_reader_t rdr = {0};
  uint64_t i = (beg1>>17);
  rdr.loc = chrmt.locs[i];
  rdr.value = chrmt.vals[i];
  rdr.index = chrmt.inds[i];
  do {
    if (rdr.value == beg1) { // found
      return (int64_t) rdr.index;
    } else if (rdr.value > beg1) {
      return 0;
    }
  } while (cr->s[rdr.loc] != 0xff && row_reader_next_loc(&rdr, cr));
  return 0;
}

