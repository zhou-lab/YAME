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

/** ---- format 4 (float / beta values with NA runs) -----
 *
 * Overall idea
 * ------------
 * Each position stores a single float (e.g. beta value). NA is encoded
 * as a special negative value (-1.0). Compression turns long stretches
 * of NA into run-length encoded markers, while keeping non-NA floats
 * in their original 32-bit representation.
 *
 *
 * Uncompressed (inflated) layout
 * ------------------------------
 * - c->compressed == 0
 * - c->fmt        == '4'
 * - c->n          = number of positions
 * - c->unit       == 4
 *
 * Bytes:
 *   s is a float array:
 *     index i -> ((float_t*)c->s)[i]
 *
 * Value semantics:
 *   - any valid float >= 0 is interpreted as a numeric value
 *   - NA is encoded as a negative float (by fmt4_read_raw, this is -1.0)
 *
 *
 * Compressed layout (RLE on NA runs)
 * ----------------------------------
 * - c->compressed == 1
 * - c->fmt        == '4'
 * - c->s is viewed as an array of uint32_t words:
 *
 *   Two record types, both 32-bit each:
 *
 *   1) NA-run marker:
 *        bit31 == 1  (MSB)
 *        bits30..0 = run length (number of consecutive NAs)
 *
 *      layout:
 *        [ 1 | run_len(31 bits) ]
 *
 *   2) Non-NA float:
 *        bit31 == 0
 *        bits30..0 = remaining bits of the IEEE-754 float
 *
 *      layout:
 *        [ 0 | float_payload(31 bits) ]
 *
 *   Compression walks the float array and:
 *     - accumulates runs of NA (negative values) into a single NA-run marker
 *     - copies each non-NA value as its 32-bit float bit pattern
 *
 *
 * Before / after schematic
 * ------------------------
 *
 *   Uncompressed view:
 *     index:   0     1     2     3     4     5   ...
 *     value:  [0.3] [NA]  [NA]  [NA]  [0.8] [NA] ...
 *
 *   Compressed view (uint32_t words):
 *     [0 | bits(0.3)]   [1 | len=3]   [0 | bits(0.8)]   [1 | len=1] ...
 *
 *
 * Decompression (fmt4_decompress)
 * --------------------------------
 *   - interpret c.s as uint32_t s0[]
 *   - for each word w:
 *       if (w >> 31)    -> NA run:
 *                           len = w & 0x7fffffff
 *                           append 'len' copies of -1.0 (NA)
 *       else            -> non-NA:
 *                           reinterpret w as float_t and append it
 *
 *   The result is a flat float array (unit=4, compressed=0) matching
 *   the original uncompressed representation.
 */

cdata_t fmt4_decompress(cdata_t c) {
  cdata_t expanded = {0};

  uint64_t i=0, m = 1<<20,n = 0, j=0, l=0;
  uint32_t *s0 = (uint32_t*) c.s;
  float_t *s = calloc(m*sizeof(float_t), 1);

  for(i=0; i< c.n>>2; ++i) {
    if (s0[i] >> 31) {
      l = s0[i]<<1>>1;
      if (n+l+10>m) {m=n+l+10; m<<=1; s = realloc(s, m*sizeof(float_t));}
      for (j=0; j<l; ++j) s[n++] = -1.0;
    } else {
      if (n+2>m) {m<<=1; s = realloc(s, m*sizeof(float_t));}
      memcpy(s+n, s0+i, sizeof(float_t));
      n++;
    }
  }

  expanded.s = (uint8_t*) s;
  expanded.n = n;
  expanded.compressed = 0;
  expanded.fmt = '4';
  expanded.unit = 4;
  return expanded;
}

static int is_float(char *s) {
  size_t i;
  for (i=0; i<strlen(s); ++i) {
    if (!isdigit(s[i]) && s[i] != '.' && s[i] != '-')
      return 0;
  }
  return 1;
}

cdata_t* fmt4_read_raw(char *fname, int verbose) {

  gzFile fh = wzopen(fname, 1);
  char *line = NULL;
  uint64_t n = 0, m=1<<22;
  float *s = calloc(m, sizeof(float));
  while (gzFile_read_line(fh, &line) > 0) {
    if (is_float(line)) {
      s[n++] = atof(line);
    } else {
      s[n++] = -1.0;
    }
    if (n+2>m) { m<<=1; s=realloc(s, m*sizeof(float)); }
  }
  free(line);
  wzclose(fh);
  if (verbose) {
    fprintf(stderr, "[%s:%d] Vector of length %"PRIu64" loaded\n", __func__, __LINE__, n);
    fflush(stderr);
  }
  cdata_t *c = calloc(sizeof(cdata_t),1);
  c->s = (uint8_t*) s;
  c->n = n;
  c->compressed = 0;
  c->fmt = '4';
  return c;
}

/* 32 bit
   1 (1bit) + run length of NA (31 bits)
   0 (1bit) + floating number (always positive) (31bit, the sign bit is always 0)
 */
void fmt4_compress(cdata_t *c) {

  uint64_t n=0, m=1<<20;
  uint32_t *s = calloc(sizeof(uint32_t), m);
  uint64_t i = 0; uint32_t l = 0;
  uint32_t *s0 = (uint32_t*) c->s;
  for (i=0, l=0; i<c->n; ++i) {
    /* either not the same as before or reach block size max */
    if (!(s0[i] & (1ul<<31)) || l+2 >= (1ul<<31)) {
      if (l > 0) {
        if (n+2>m) { m<<=1; s = realloc(s, m*sizeof(uint32_t));}
        s[n++] = ((1<<31) | l);
        l = 0;
      }

      if (!(s0[i] & (1ul<<31))) {
        if (n+2>m) { m<<=1; s = realloc(s, m*sizeof(uint32_t));}
        memcpy(s+n, s0+i, sizeof(float_t));
        n++;
        l = 0;
      }
    } else {
      ++l;
    }
  }
  /* the last rle */
  if (l > 0) {
    if (n+2>m) { m<<=1; s = realloc(s, m*sizeof(uint32_t));}
    s[n++] = ((1<<31) | l);
  }
  
  free(c->s);
  c->s = (uint8_t*) s;
  c->n = n*4;
  c->compressed = 1;
}

stats_t* summarize1_queryfmt4(
  cdata_t *c, cdata_t *c_mask, uint64_t *n_st, char *sm, char *sq, config_t *config) {

  stats_t *st = NULL;
  float_t *vals = (float_t*) c->s;

  if (c_mask->n == 0) {          // no mask

    *n_st = 1;
    st = calloc(1, sizeof(stats_t));
    st[0].n_u = c->n;

    for (uint64_t i = 0; i < c->n; ++i) {
      double b = vals[i];
      if (b >= 0.0) {            // non-NA beta
        st[0].n_q++;
        st[0].n_o++;
        st[0].sum_beta += b;
      }
    }

    st[0].sm = strdup(sm);
    st[0].sq = strdup(sq);
    st[0].beta = st[0].n_o ? (st[0].sum_beta / st[0].n_o) : NAN;

  } else if (c_mask->fmt <= '1') { // binary mask

    if (c_mask->n != c->n) {
      fprintf(stderr, "[%s:%d] mask (N=%"PRIu64") and query (N=%"PRIu64") are of different lengths.\n",
              __func__, __LINE__, c_mask->n, c->n);
      fflush(stderr);
      exit(1);
    }

    *n_st = 1;
    st = calloc(1, sizeof(stats_t));
    st[0].n_u = c->n;

    for (uint64_t i = 0; i < c->n; ++i) {
      double b = vals[i];
      if (b >= 0.0) st[0].n_q++;
      if (FMT0_IN_SET(*c_mask, i)) {
        st[0].n_m++;
        if (b >= 0.0) {
          st[0].n_o++;
          st[0].sum_beta += b;
        }
      }
    }

    st[0].sm = strdup(sm);
    st[0].sq = strdup(sq);
    st[0].beta = st[0].n_o ? (st[0].sum_beta / st[0].n_o) : NAN;

  } else if (c_mask->fmt == '6') { // binary mask with universe

    if (c_mask->n != c->n) {
      fprintf(stderr, "[%s:%d] mask (N=%"PRIu64") and query (N=%"PRIu64") are of different lengths.\n",
              __func__, __LINE__, c_mask->n, c->n);
      fflush(stderr);
      exit(1);
    }

    *n_st = 1;
    st = calloc(1, sizeof(stats_t));
    st[0].n_u = c->n;

    for (uint64_t i = 0; i < c->n; ++i) {
      double b = vals[i];
      if (b >= 0.0) st[0].n_q++;
      if (FMT6_IN_UNI(*c_mask, i) && FMT6_IN_SET(*c_mask, i)) {
        st[0].n_m++;
        if (b >= 0.0) {
          st[0].n_o++;
          st[0].sum_beta += b;
        }
      }
    }

    st[0].sm = strdup(sm);
    st[0].sq = strdup(sq);
    st[0].beta = st[0].n_o ? (st[0].sum_beta / st[0].n_o) : NAN;

  } else if (c_mask->fmt == '2') { // state mask

    if (c_mask->n != c->n) {
      fprintf(stderr, "[%s:%d] mask (N=%"PRIu64") and query (N=%"PRIu64") are of different lengths.\n",
              __func__, __LINE__, c_mask->n, c->n);
      fflush(stderr);
      exit(1);
    }
    if (!c_mask->aux) fmt2_set_aux(c_mask);
    f2_aux_t *aux = (f2_aux_t*) c_mask->aux;

    *n_st = aux->nk;
    st = calloc((*n_st), sizeof(stats_t));
    uint64_t nq = 0;

    for (uint64_t i = 0; i < c->n; ++i) {
      uint64_t index = f2_get_uint64(c_mask, i);
      if (index >= (*n_st)) {
        fprintf(stderr, "[%s:%d] State data is corrupted.\n", __func__, __LINE__);
        fflush(stderr);
        exit(1);
      }
      double b = vals[i];
      if (b >= 0.0) {
        st[index].n_o++;
        st[index].sum_beta += b;
        nq++;
      }
      st[index].n_m++;
    }

    for (uint64_t k = 0; k < (*n_st); ++k) {
      st[k].n_q = nq;
      st[k].n_u = c->n;
      st[k].beta = st[k].n_o ? (st[k].sum_beta / st[k].n_o) : NAN;
      if (config->section_name) {
        kstring_t tmp = {0};
        ksprintf(&tmp, "%s-%s", sm, aux->keys[k]);
        st[k].sm = tmp.s;
      } else {
        st[k].sm = strdup(aux->keys[k]);
      }
      st[k].sq = strdup(sq);
    }

  } else {                      // other masks
    fprintf(stderr, "[%s:%d] Mask format %c unsupported.\n",
            __func__, __LINE__, c_mask->fmt);
    fflush(stderr);
    exit(1);
  }

  return st;
}
