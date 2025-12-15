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

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "cdata.h"
#include "summary.h"

/** ---- format 6 (binary with universe bit) -----
 *
 * Overall idea
 * ------------
 * Each position i stores a 2-bit state:
 *
 *   - "universe" bit: whether this position is in the universe
 *   - "set" bit     : value within the universe (0/1)
 *
 * Semantics per 2-bit code (UNIVERSE:SET):
 *   00  -> NA / outside universe
 *   10  -> in universe, value = 0
 *   11  -> in universe, value = 1
 *
 * Macros (see cdata.h):
 *   FMT6_IN_UNI(c,i)   -> true if universe bit is 1
 *   FMT6_IN_SET(c,i)   -> true if set bit is 1
 *   FMT6_SET_NA(c,i)   -> set to 00 (outside universe)
 *   FMT6_SET0(c,i)     -> set to 10 (in universe, value 0)
 *   FMT6_SET1(c,i)     -> set to 11 (in universe, value 1)
 *
 *
 * Inflated (uncompressed) layout
 * ------------------------------
 * - c->compressed == 0
 * - c->fmt        == '6'
 * - c->n          = number of logical positions
 * - c->unit       = 2 (conceptually 2 bits per position)
 *
 * Bits are densely packed: 4 positions per byte.
 *
 *   Byte/bit layout:
 *     byte b = c->s[i >> 2]
 *     offset = (i & 0x3) * 2
 *
 *     set-bit      = (b >> offset)     & 1
 *     universe-bit = (b >> (offset+1)) & 1
 *
 *   Schematic (one byte, 4 positions):
 *     [ u3 s3 | u2 s2 | u1 s1 | u0 s0 ]
 *      bits7-6 bits5-4 bits3-2 bits1-0
 *
 *
 * Compressed / decompressed layout
 * --------------------------------
 * Format 6 does not introduce any additional coding. The on-disk and
 * in-memory byte layouts are identical; "compression" is just a flag.
 *
 * - fmt6_compress(c):
 *     sets c->compressed = 1, but does not change c->s
 *
 * - fmt6_decompress(c):
 *     copies the packed bytes into a new cdata_t with:
 *       compressed = 0, fmt = '6', same n and byte content
 *
 *
 * Before / after schematic
 * ------------------------
 *
 *   Inflated (compressed == 0):
 *     c->s = [ packed SUSUSUSU bytes ]
 *     c->n = number of positions
 *
 *   Compressed (compressed == 1):
 *     same byte stream:
 *     [ packed SUSUSUSU bytes ]
 *
 * In other words, format 6 is always stored as 2-bit fields packed into
 * bytes; the "compression" bit only signals whether fmt6_decompress()
 * needs to be called, not a different encoding scheme.
 */

static int is_int(char *s) {
  size_t i;
  for (i=0; i<strlen(s); ++i) {
    if (!isdigit(s[i])) return 0;
  }
  return 1;
}

/* every byte is SUSUSUSU */
cdata_t* fmt6_read_raw(char *fname, int verbose) {

  gzFile fh = wzopen(fname, 1);
  char *line = NULL;
  uint8_t *s = NULL; uint64_t n = 0;
  char **fields; int nfields;
  while (gzFile_read_line(fh, &line) > 0) {
    line_get_fields(line, "\t", &fields, &nfields);
    if (nfields < 2) wzfatal("Number of fields < 2. Abort.");
    if (!is_int(fields[1])) wzfatal("The 2nd column must be integers.");
    s = realloc(s, n/4 + 1);
    // Binarize fields[0] to the first bit and fields[1] to the second bit
    if (n%4 == 0) s[n/4] = 0; // initialize
    if (atoi(fields[1])) s[n/4] |= 1<<((n%4)*2 + 1); // universe
    if (strcmp(fields[0], "1") == 0) s[n/4] |= 1<<((n%4)*2); // set
    n++;
    free_fields(fields, nfields);
  }
  free(line);
  wzclose(fh);
  if (verbose) {
    fprintf(stderr, "[%s:%d] Data of length %"PRIu64" loaded\n", __func__, __LINE__, n);
    fflush(stderr);
  }
  cdata_t *c = calloc(sizeof(cdata_t),1);
  c->s = s;
  c->n = n;
  c->compressed = 0;
  c->fmt = '6';
  c->unit = 2;
  return c;
}

// do nothing
void fmt6_compress(cdata_t *c) {
  c->compressed = 1;
}

// just copy, nothing else, one can just flip compressed bit if possible
cdata_t fmt6_decompress(const cdata_t c) {
  cdata_t expanded = {0};
  expanded.s = malloc(cdata_nbytes(&c));
  memcpy(expanded.s, c.s, cdata_nbytes(&c));
  expanded.n = c.n;
  expanded.compressed = 0;
  expanded.fmt = '6';
  expanded.unit = 2;
  return expanded;
}

// as set/universe
static stats_t* summarize1_queryfmt6_SU(
  cdata_t *c, cdata_t *c_mask, uint64_t *n_st, char *sm, char *sq, config_t *config) {
  
  stats_t *st = NULL;
  if (c_mask->n == 0) {          // no mask
    
    *n_st = 1;
    st = calloc(1, sizeof(stats_t));
    for (uint64_t i=0; i<c->n; ++i) {
      if (FMT6_IN_UNI(*c,i)) {
        st[0].n_u++;
        st[0].n_m++;
        st[0].sum_depth++;
        if (FMT6_IN_SET(*c,i)) {
          st[0].n_q++;
          st[0].n_o++;
        }
      }
    }
    st[0].sm = strdup(sm);
    st[0].sq = strdup(sq);
    st[0].beta = (double) st[0].n_q / st[0].n_u;
    
  } else if (c_mask->fmt <= '1') { // binary mask

    if (c_mask->n != c->n) wzfatal("[%s:%d] mask (N=%"PRIu64") and query (N=%"PRIu64") are of different lengths.\n", __func__, __LINE__, c_mask->n, c->n);
    
    *n_st = 1;
    st = calloc(1, sizeof(stats_t));
    for (size_t i=0; i<c->n; ++i) {
      if (FMT6_IN_UNI(*c,i)) {
        st[0].n_u++;
        int in_q = FMT6_IN_SET(*c,i);
        int in_m = FMT0_IN_SET(*c_mask,i);
        if (in_q) st[0].n_q++;
        if (in_m) st[0].n_m++;
        if (in_q && in_m) st[0].n_o++;
      }
    }
    st[0].sm = strdup(sm);
    st[0].sq = strdup(sq);
    st[0].beta = (double) st[0].n_o / st[0].n_m;

  } else if (c_mask->fmt == '2') { // state mask

    if (c_mask->n != c->n) wzfatal("[%s:%d] mask (N=%"PRIu64") and query (N=%"PRIu64") are of different lengths.\n", __func__, __LINE__, c_mask->n, c->n);

    if (!c_mask->aux) fmt2_set_aux(c_mask);
    f2_aux_t *aux = (f2_aux_t*) c_mask->aux;
    *n_st = aux->nk;
    st = calloc((*n_st), sizeof(stats_t));
    uint64_t nq = 0, nu = 0;
    for (uint64_t i=0; i<c->n; ++i) {
      uint64_t index = f2_get_uint64(c_mask, i);
      if (index >= (*n_st)) wzfatal("[%s:%d] State data is corrupted.\n", __func__, __LINE__);
      if (FMT6_IN_UNI(*c,i)) {
        nu++;
        if (FMT6_IN_SET(*c,i)) {
          nq++;
          st[index].n_o++;
        }
        st[index].n_m++;
      }
    }
    for (uint64_t k=0; k < (*n_st); ++k) {
      st[k].n_q = nq;
      st[k].n_u = nu;
      if (config->section_name) {
        kstring_t tmp = {0};
        ksprintf(&tmp, "%s-%s", sm, aux->keys[k]);
        st[k].sm = tmp.s;
      } else {
        st[k].sm = strdup(aux->keys[k]);
      }
      st[k].sq = strdup(sq);
      st[k].beta = (double) st[k].n_o / st[k].n_m;
    }
  } else if (c_mask->fmt == '6') { // binary mask with universe

    if (c_mask->n != c->n) wzfatal("[%s:%d] mask (N=%"PRIu64") and query (N=%"PRIu64") are of different lengths.\n", __func__, __LINE__, c_mask->n, c->n);
    
    *n_st = 1;
    st = calloc(1, sizeof(stats_t));
    for (size_t i=0; i<c->n; ++i) {
      if (FMT6_IN_UNI(*c,i) && FMT6_IN_UNI(*c_mask, i)) {
        st[0].n_u++;
        int in_q = FMT6_IN_SET(*c,i);
        int in_m = FMT6_IN_SET(*c_mask,i);
        if (in_q) st[0].n_q++;
        if (in_m) st[0].n_m++;
        if (in_q && in_m) st[0].n_o++;
      }
    }
    st[0].sm = strdup(sm);
    st[0].sq = strdup(sq);
    st[0].beta = (double) st[0].n_o / st[0].n_m;
    
  } else {                      // other masks
    wzfatal("[%s:%d] Mask format %c unsupported.\n", __func__, __LINE__, c_mask->fmt);
  }
  return st;
}

// as quaternary data
static stats_t* summarize1_queryfmt6_2bit(
  cdata_t *c, cdata_t *c_mask, uint64_t *n_st, char *sm, char *sq, config_t *config) {
  
  stats_t *st = NULL;
  if (c_mask->n == 0) {          // no mask
    
    *n_st = 4;
    st = calloc(4, sizeof(stats_t));
    for (uint64_t i=0; i<c->n; ++i) {
      st[FMT6_2BIT(*c,i)].n_q++;
    }
    for (uint8_t k=0; k<4; ++k) {
      st[k].n_u = c->n;
      st[k].n_m = c->n;
      st[k].n_o = st[k].n_q;
      st[k].sm = strdup(sm);
      kstring_t tmp = {0};
      ksprintf(&tmp, "%s|%u", sq, k);
      st[k].sq = tmp.s;
      st[k].beta = 1.0;
    }
    
  } else if (c_mask->fmt <= '1') { // binary mask

    if (c_mask->n != c->n) wzfatal("[%s:%d] mask (N=%"PRIu64") and query (N=%"PRIu64") are of different lengths.\n", __func__, __LINE__, c_mask->n, c->n);

    *n_st = 4;
    uint64_t *cnts = calloc(*n_st, sizeof(uint64_t));
    uint64_t *cnts_q = calloc(*n_st, sizeof(uint64_t));
    uint64_t n_m = 0; // sum of 1s in mask
    for (uint64_t i=0; i<c->n; ++i) {
      if (FMT0_IN_SET(*c_mask, i)) {
        n_m++;
        cnts[FMT6_2BIT(*c,i)]++;
      }
      cnts_q[FMT6_2BIT(*c,i)]++;
    }
    st = calloc(*n_st, sizeof(stats_t));
    for (uint8_t k=0; k<*n_st; ++k) {
      // only report masked
      st[k].n_u = c->n;
      st[k].n_q = cnts_q[k];
      st[k].n_o = cnts[k];
      st[k].n_m = n_m;
      kstring_t tmp1 = {0};
      ksprintf(&tmp1, "%s|1", sm); // the |1 means "masked"
      st[k].sm = tmp1.s;
      kstring_t tmp2 = {0};
      ksprintf(&tmp2, "%s|%u", sq, k);
      st[k].sq = tmp2.s;
    }
    free(cnts);
    free(cnts_q);
    
  } else if (c_mask->fmt == '2') { // state mask

    if (c_mask->n != c->n) wzfatal("[%s:%d] mask (N=%"PRIu64") and query (N=%"PRIu64") are of different lengths.\n", __func__, __LINE__, c_mask->n, c->n);

    if (!c_mask->aux) fmt2_set_aux(c_mask);
    f2_aux_t *aux = (f2_aux_t*) c_mask->aux;
    *n_st = aux->nk * 4;
    st = calloc((*n_st), sizeof(stats_t));
    uint64_t nu = 0;
    for (uint64_t i=0; i<c->n; ++i) {
      uint64_t index = f2_get_uint64(c_mask, i);
      if (index >= (*n_st)) wzfatal("[%s:%d] State data is corrupted.\n", __func__, __LINE__);
      nu++;
      st[index*4+FMT6_2BIT(*c,i)].n_o++;
      for (uint64_t k1 = 0; k1 < aux->nk; ++k1) st[k1*4+FMT6_2BIT(*c,i)].n_q++;
      for (uint8_t k2 = 0; k2 < 4; ++k2) st[index*4+k2].n_m++;
    }
    for (uint64_t k1=0; k1 < aux->nk; ++k1) {
      for (uint8_t k2=0; k2 < 4; ++k2) {
        uint64_t k = k1*4 + k2;
        st[k].n_u = nu;
        if (config->section_name) {
          kstring_t tmp = {0};
          ksprintf(&tmp, "%s-%s", sm, aux->keys[k1]);
          st[k].sm = tmp.s;
        } else {
          st[k].sm = strdup(aux->keys[k1]);
        }
        
        kstring_t tmp = {0};
        ksprintf(&tmp, "%s|%u", sq, k2);
        st[k].sq = tmp.s;
        st[k].beta = (double) st[k].n_o / st[k].n_m;
      }
    }
    
  } else if (c_mask->fmt == '6') { // binary mask with universe

    if (c_mask->n != c->n) wzfatal("[%s:%d] mask (N=%"PRIu64") and query (N=%"PRIu64") are of different lengths.\n", __func__, __LINE__, c_mask->n, c->n);

    *n_st = 4;
    uint64_t *cnts = calloc(*n_st, sizeof(uint64_t));
    uint64_t *cnts_q = calloc(*n_st, sizeof(uint64_t));
    uint64_t n_m = 0; // sum of 1s in mask
    for (uint64_t i=0; i<c->n; ++i) {
      if (FMT6_IN_UNI(*c_mask, i) && FMT6_IN_SET(*c_mask, i)) {
        n_m++;
        cnts[FMT6_2BIT(*c,i)]++;
      }
      cnts_q[FMT6_2BIT(*c,i)]++;
    }
    st = calloc(*n_st, sizeof(stats_t));
    for (uint8_t k=0; k<*n_st; ++k) {
      // only report masked
      st[k].n_u = c->n;
      st[k].n_q = cnts_q[k];
      st[k].n_o = cnts[k];
      st[k].n_m = n_m;
      kstring_t tmp1 = {0};
      ksprintf(&tmp1, "%s|1", sm); // the |1 means "masked"
      st[k].sm = tmp1.s;
      kstring_t tmp2 = {0};
      ksprintf(&tmp2, "%s|%u", sq, k);
      st[k].sq = tmp2.s;
    }
    free(cnts);
    free(cnts_q);
    
  } else {                      // other masks
    wzfatal("[%s:%d] Mask format %c unsupported.\n", __func__, __LINE__, c_mask->fmt);
  }
  return st;
}

stats_t* summarize1_queryfmt6(
  cdata_t *c, cdata_t *c_mask, uint64_t *n_st, char *sm, char *sq, config_t *config) {

  if (config->f6_as_2bit)
    return summarize1_queryfmt6_2bit(c, c_mask, n_st, sm, sq, config);
  else
    return summarize1_queryfmt6_SU(c, c_mask, n_st, sm, sq, config);
}
