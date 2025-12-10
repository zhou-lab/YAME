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

/** ----- format 5 (ternary 0/1/2 with compact runs) -----
 *
 * Overall idea
 * ------------
 * Each position stores a ternary value:
 *   0 = "off"
 *   1 = "on"
 *   2 = "NA" (missing)
 *
 * Uncompressed, we keep one byte per position (0,1,2).
 * Compressed, we encode:
 *   - runs of NA (=2) as a single byte (NA-run)
 *   - runs of 0/1 as packed 0/1 bits in a byte (up to 4 values per byte)
 *
 *
 * Uncompressed (inflated) layout
 * ------------------------------
 * - c->compressed == 0
 * - c->fmt        == '5'
 * - c->n          = number of positions
 * - c->unit       = 1
 *
 * Bytes:
 *   s[i] ∈ {0,1,2} for position i
 *
 *
 * Compressed layout
 * -----------------
 * - c->compressed == 1
 * - c->fmt        == '5'
 * - c->s is a sequence of *bytes* of two possible types:
 *
 *   1) NA-run byte (MSB = 0):
 *        bit7 = 0
 *        bits6..0 = run length L (number of consecutive NA==2 values)
 *
 *      Decompression:
 *        repeat value 2, L times.
 *
 *
 *   2) packed 0/1 byte (MSB = 1):
 *        bit7 = 1  (signals this is a packed 0/1 block)
 *        remaining 7 bits are split into up to 4 (flag,value) pairs:
 *
 *          [ f3 v3 | f2 v2 | f1 v1 | f0 v0 ]
 *           ^^^^^^  ^^^^^^  ^^^^^^  ^^^^^^
 *           bits6-5 bits4-3 bits2-1 bits0-(-1)
 *
 *        For each pair (fk, vk):
 *          fk = 1  -> this slot is used, emit value vk ∈ {0,1}
 *          fk = 0  -> sentinel, stop reading further pairs
 *
 *      Decompression:
 *        starting from the highest offset (6,4,2,0):
 *          if (byte & (0x2 << offset))   // fk == 1
 *              emit (byte >> offset) & 0x1  // vk
 *          else
 *              break   // remaining slots unused
 *
 *
 * Before / after schematic
 * ------------------------
 *
 *   Uncompressed:
 *     index:  0  1  2  3  4  5  6  ...
 *     data:  [2][2][2][0][1][0][2]...
 *
 *   Compressed (bytes):
 *     [ 0b00000011 ]           // NA-run of length 3
 *     [ packed 0/1 for 0,1,0 ] // MSB=1, three used slots
 *     [ 0b00000001 ]           // NA-run of length 1
 *
 *
 * Decompression (fmt5_decompress)
 * -------------------------------
 *   for each byte b in c.s:
 *     if (b & 0x80)  // MSB=1 -> packed 0/1
 *       walk flag/value pairs at bit offsets 6,4,2,0
 *       emit 0/1 until sentinel (flag==0)
 *     else           // MSB=0 -> NA-run
 *       L = b        // 7-bit length
 *       emit value 2, L times
 *
 *   Result is a flat array s[] of length n with values in {0,1,2},
 *   i.e. the original uncompressed representation (unit = 1).
 */

cdata_t fmt5_decompress(cdata_t c) {
  cdata_t expanded = {0};
  uint64_t i = 0, m = 1<<20, n = 0, j = 0;
  uint8_t *s = malloc(m*sizeof(uint8_t));

  for (i=0; i<c.n; ++i) {
    if (c.s[i] & (1<<7)) {
      int offset = 6;
      if (n+2>m) {m<<=1; s = realloc(s, m*sizeof(uint8_t));}
      for (offset = 6; offset >= 0; offset -= 2) {
        if ((c.s[i]>>offset) & 0x2) {
          s[n++] = ((c.s[i]>>offset) & 0x1);
        } else {
          break;
        }
      }
    } else {
      if (n+c.s[i]+10>m) {m=n+c.s[i]+10; m<<=1; s = realloc(s, m*sizeof(uint8_t));}
      for (j=0; j < c.s[i]; ++j) s[n++] = 2;
    }
  }

  expanded.s = (uint8_t*) s;
  expanded.n = n;
  expanded.compressed = 0;
  expanded.fmt = '5';
  expanded.unit = 1;
  return expanded;
}

/* the input has only 0,1,2 */
cdata_t* fmt5_read_raw(char *fname, int verbose) {

  gzFile fh = wzopen(fname, 1);
  char *line = NULL;
  uint64_t n = 0, m=1<<22;
  uint8_t *s = calloc(m, 1);
  while (gzFile_read_line(fh, &line) > 0) {
    if (line[0] == '0' || line[0] == '1') s[n++] = line[0]-'0';
    else s[n++] = 2;
    if (n+2>m) { m<<=1; s=realloc(s,m); }
  }
  free(line);
  wzclose(fh);
  if (verbose) {
    fprintf(stderr, "[%s:%d] Vector of length %llu loaded\n", __func__, __LINE__, n);
    fflush(stderr);
  }
  cdata_t *c = calloc(sizeof(cdata_t),1);
  c->s = (uint8_t*) s;
  c->n = n;
  c->compressed = 0;
  c->fmt = '5';
  return c;
}

/*
  8 bits = 0 (1bit) | run length of NA (7bits)
  8 bits = 1 (1bit)|value (1bit) + 1 (1bit)|value (1bit) + ...
*/
void fmt5_compress(cdata_t *c) {
  uint64_t n = 0;
  uint8_t *s = NULL;
  uint64_t i = 0; uint16_t l = 0; int last = 0; uint8_t u = 0; int offset = 6;
  for (i=0, l=0; i<c->n; ++i) {
    if (c->s[i] == 0 || c->s[i] == 1) { /* 0 or 1 */
      u |= (1<<(offset+1));
      u |= (c->s[i]<<offset);
      offset -= 2;
      if (last <= 1) {            /* 0/1 > 0/1 */
        if (offset < 0) {
          s = realloc(s, n+1);
          s[n++] = u;
          u = 0; offset = 6;
        }
      } else if (l > 0) {       /* 2 > 0/1 */
        s = realloc(s, n+1);
        s[n++] = l;
        l = 0;
      }
      last = 1;
    } else {                    /* neither 0 nor 1, for missing value */
      if (last == 1 && u != 0) {               /* 0/1 > 2 */
        if (offset >= 0) u |= (0<<(offset+1)); /* add sentinel */
        s = realloc(s, n+1);
        s[n++] = u;
        u = 0; offset = 6;
      }
      l++;
      if (l+2 >= 1<<7) {        /* too many NA start a new count */
        s = realloc(s, n+1);
        s[n++] = l;
        l = 0;
      }
      last = 2;
    }
  }

  if (last == 1 && u != 0) {
    s = realloc(s, n+1);
    s[n++] = u;
  } else if (last == 2 && l > 0) {
    s = realloc(s, n+1);
    s[n++] = l;
  }

  free(c->s);
  c->s = s;
  c->n = n;
  c->compressed = 1;
}
