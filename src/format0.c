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

/** Format 0 (binary presence/absence) bit-packing
 *
 * Inflated (uncompressed) representation
 * --------------------------------------
 * - c->compressed == 0
 * - c->n = number of CpG sites
 * - Each site stores a single binary value (0 or 1).
 * - Bits are densely packed: 8 CpGs per byte.
 *
 *   Byte layout:
 *       s[i >> 3]  holds CpGs at indices  (i&0x7) = 0..7
 *
 *   Bit test:
 *       FMT0_IN_SET(c, i)   -> true if bit (i&0x7) is 1
 *   Bit set:
 *       FMT0_SET(c, i)      -> set bit (i&0x7) to 1
 *
 *   Schematic:
 *       Inflated array of bits:
 *         index:    0 1 2 3 4 5 6 7 | 8 9 10 11 ...
 *         stored in bytes as:
 *         byte0:   [b7 b6 b5 b4 b3 b2 b1 b0]
 *         byte1:   [b15 ...                ]
 *
 *
 * Compressed representation
 * -------------------------
 * Format 0 does *not* define its own compression scheme.
 *
 * Rules:
 *   - If data is already in fmt0 form, it is considered uncompressed.
 *   - When c->compressed==1 and fmt=='0', it simply means:
 *         "data may have originated from a conversion (e.g., fmt1->fmt0,
 *          fmt3->fmt0) but storage remains the same packed-bit layout."
 *   - No variable-length encoding or RLE is used for fmt0.
 *
 * In other words, the on-disk and in-memory representations are identical:
 *
 *   Before compression:  "packed bits"
 *   After compression:   "packed bits (unchanged)"
 *
 *
 * Before / After schematic
 * ------------------------
 *
 *   Inflated (RAM):
 *       s = [ packed bytes ]
 *       c->n = number of logical CpG bits
 *
 *   Compressed (disk):
 *       identical byte stream:
 *       [ packed bytes ]
 *
 * Decompression
 * -------------
 * fmt0_decompress() simply copies the packed bytes into the expanded buffer;
 * no decoding is required.
 *
 */

/* 8 bit for 8 cpgs, each is binary */
cdata_t* fmt0_read_raw(char *fname, int verbose) {

  gzFile fh = wzopen(fname, 1);
  char *line = NULL;
  uint64_t n = 0, m=1<<22;
  uint8_t *s = calloc(m, 1);
  while (gzFile_read_line(fh, &line) > 0) {
    if (line[0] != '0') {
      s[n>>3] |= (1<<(n&0x7));
    }
    n++;
    if (n+2>m) { m<<=1; s=realloc(s,m); }
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
  c->fmt = '0';
  c->unit = 1;
  return c;
}

/* just copy, nothing done */
cdata_t fmt0_decompress(const cdata_t c) {
  cdata_t expanded = c;
  expanded.s = calloc(cdata_nbytes(&c), 1);
  memcpy(expanded.s, c.s, cdata_nbytes(&c));
  expanded.unit = 1;
  expanded.n = c.n;
  expanded.compressed = 0;
  expanded.fmt = '0';
  return expanded;
}

/**
 * - Input format is '0': no further operation is performed.
 * - Input format is '1': if the value is 0, return 0 else 1.
 * - Input format is '3': if the M+U is 0, return 0 else 1
 * - Other input formats are not allowed.
 */
cdata_t fmt3_decompress(const cdata_t c);
void convertToFmt0(cdata_t *c) {
  cdata_t c_out = {0};
  switch (c->fmt) {
  case '0': return;
  case '1': { // assume it's compressed.
    c_out.fmt = '0';
    c_out.compressed = 1;
    c_out.n=0;
    uint64_t i;
    for (i=0; i<c->n/3; ++i) {
      c_out.n += *((uint16_t*) (c->s+i*3+1));
    }
    c_out.s = calloc((c_out.n>>3)+1, 1);
    size_t sum; uint16_t l=0;
    for (i=0, sum=0; i<c->n/3; ++i, sum+=l) {
      l = *((uint16_t*) (c->s+i*3+1));
      if (c->s[i*3] > '0') {
        for(size_t j=sum; j<sum+l; ++j) {
          FMT0_SET(c_out, j);
        }
      }
    }
    break;
  }
  case '3': {
    cdata_t expanded = fmt3_decompress(*c);
    /* fmt3_decompress(c, &expanded); */

    c_out.fmt = '0';
    c_out.compressed = 1;
    c_out.n = expanded.n;
    c_out.s = calloc((c_out.n>>3)+1,1);
    for (uint64_t i=0; i<expanded.n; ++i) {
      uint64_t mu = f3_get_mu(&expanded, i);
      if (mu > 0) { /* equivalent to: 1 if M+U > 0 else 0 */
        c_out.s[i>>3] |= (1<<(i&0x7));
      }
    }
    free(expanded.s);
    break;
  }
  default: wzfatal("Format %c unsupported.\n", c->fmt);
  }
  free(c->s);
  *c = c_out;
}

// when using format 0 as query, it is treated as a set rather than a binary states.
// TODO: maybe we should have a binary mode like the quaternary mode for format 6.
stats_t* summarize1_queryfmt0(
  cdata_t *c, cdata_t *c_mask, uint64_t *n_st, char *sm, char *sq, config_t *config) {

  stats_t *st = NULL;
  if (c_mask->n == 0) {          // no mask
    
    *n_st = 1;
    st = calloc(1, sizeof(stats_t));
    st[0].n_u = c->n;
    st[0].n_q = bit_count(c[0]);
    st[0].sm = strdup(sm);
    st[0].sq = strdup(sq);
    
  } else if (c_mask->fmt <= '1') { // binary mask

    if (c_mask->n != c->n) {
      fprintf(stderr, "[%s:%d] mask (N=%"PRIu64") and query (N=%"PRIu64") are of different lengths.\n", __func__, __LINE__, c_mask->n, c->n);
      fflush(stderr);
      exit(1);
    }
    
    *n_st = 1;
    st = calloc(1, sizeof(stats_t));
    st[0].n_u = c->n;
    st[0].n_q = bit_count(c[0]);
    st[0].n_m = bit_count(c_mask[0]);
    cdata_t tmp = {0};
    tmp.s = malloc((c->n>>3)+1); tmp.n = c->n;
    memcpy(tmp.s, c->s, (c->n>>3)+1);
    for (uint64_t i=0; i<(tmp.n>>3)+1; ++i) tmp.s[i] &= c_mask->s[i];
    st[0].n_o = bit_count(tmp);
    free(tmp.s);
    st[0].sm = strdup(sm);
    st[0].sq = strdup(sq);

  } else if (c_mask->fmt == '2') { // state mask

    if (c_mask->n != c->n) {
      fprintf(stderr, "[%s:%d] mask (N=%"PRIu64") and query (N=%"PRIu64") are of different lengths.\n", __func__, __LINE__, c_mask->n, c->n);
      fflush(stderr);
      exit(1);
    }
    if (!c_mask->aux) fmt2_set_aux(c_mask);
    f2_aux_t *aux = (f2_aux_t*) c_mask->aux;
    *n_st = aux->nk;
    st = calloc((*n_st), sizeof(stats_t));
    uint64_t nq=0;
    for (uint64_t i=0; i<c->n; ++i) {
      uint64_t index = f2_get_uint64(c_mask, i);
      if (index >= (*n_st)) {
        fprintf(stderr, "[%s:%d] State data is corrupted.\n", __func__, __LINE__);
        fflush(stderr);
        exit(1);
      }
      if (FMT0_IN_SET(*c, i)) {
        st[index].n_o++;
        nq++;
      }
      st[index].n_m++;
    }
    for (uint64_t k=0; k < (*n_st); ++k) {
      st[k].n_q = nq;
      st[k].n_u = c->n;
      if (config->section_name) {
        kstring_t tmp = {0};
        ksprintf(&tmp, "%s-%s", sm, aux->keys[k]);
        st[k].sm = tmp.s;
      } else {
        st[k].sm = strdup(aux->keys[k]);
      }
      st[k].sq = strdup(sq);
    }

  } else if (c_mask->fmt == '6') { // binary mask with universe

    if (c_mask->n != c->n) {
      fprintf(stderr, "[%s:%d] mask (N=%"PRIu64") and query (N=%"PRIu64") are of different lengths.\n", __func__, __LINE__, c_mask->n, c->n);
      fflush(stderr);
      exit(1);
    }
    
    *n_st = 1;
    stats_t st1 = {0};
    for (uint64_t i=0; i<c->n; ++i) {
      if (FMT6_IN_UNI(*c_mask, i)) {
        st1.n_u++;
        int in_q = FMT0_IN_SET(*c, i);
        int in_m = FMT6_IN_SET(*c_mask, i);
        if (in_q) st1.n_q++;
        if (in_m) st1.n_m++;
        if (in_q && in_m) st1.n_o++;
      }
    }
    st = calloc(1, sizeof(stats_t));
    st[0] = st1;
    st[0].sm = strdup(sm);
    st[0].sq = strdup(sq);

  } else {                      // other masks
    fprintf(stderr, "[%s:%d] Mask format %c unsupported.\n", __func__, __LINE__, c_mask->fmt);
    fflush(stderr);
    exit(1);
  }
  return st;
}

