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

/** ---- format 1 (integer / ASCII with RLE) -----
 *
 * Overall idea
 * ------------
 * Each row stores a single-byte value per position (typically ASCII digits
 * '0'..'9', but not restricted). Format 1 supports a run-length encoded
 * (RLE) representation for long runs of identical values.
 *
 *
 * Uncompressed (inflated) layout
 * ------------------------------
 * - c->compressed == 0
 * - c->fmt        == '1'
 * - c->n          = number of positions
 * - c->unit       = 1 byte per position
 *
 * Bytes:
 *   s[0], s[1], ..., s[c->n-1]
 *
 *   index i -> value s[i] (1 byte, arbitrary 0–255)
 *
 *
 * Compressed layout (RLE)
 * -----------------------
 * - c->compressed == 1
 * - c->fmt        == '1'
 * - c->s is a sequence of fixed-size RLE records:
 *
 *   Each record is 3 bytes:
 *
 *     [ value (1 byte) | length (2 bytes, uint16 little-endian) ]
 *
 *   where:
 *     value  = the repeated byte
 *     length = number of times 'value' is repeated (1..32767)
 *
 *   Schematic:
 *     [ v0 | len0_lo len0_hi ]
 *     [ v1 | len1_lo len1_hi ]
 *     ...
 *
 * fmt1_compress() walks the uncompressed array and writes out these
 * (value,length) triples, splitting runs when the 16-bit length limit is hit.
 *
 *
 * Before / after schematic
 * ------------------------
 *
 *   Uncompressed view:
 *     index:  0    1    2    3    4    5   ...
 *     data:  [A]  [A]  [A]  [B]  [B]  [C]  ...
 *
 *   Compressed view:
 *     [ A | 3 ]   [ B | 2 ]   [ C | 1 ]
 *
 * fmt1_decompress():
 *   - iterates over the RLE records (value,length)
 *   - expands each into 'length' copies of 'value'
 *   - writes them into a flat 1-byte-per-position array (unit=1),
 *     restoring the uncompressed representation.
 *
 *
 * Note on binary vectors (helper fmta_tryBinary2byteRLE_ifsmaller)
 * ----------------------------------------------------------------
 * For fmt0 bit-vectors, fmta_tryBinary2byteRLE_ifsmaller() can convert them
 * into fmt1 RLE triples (with values '0'/'1') when that representation is
 * smaller than the original packed-bit layout.
 */

cdata_t fmt1_decompress(cdata_t c) {
  cdata_t expanded = {0};
  uint64_t i=0, j=0, n=0, m=1<<20;
  uint8_t *s = realloc(expanded.s, m);
  for (i=0; i<c.n; i+=3) {
    uint16_t l = ((uint16_t*) (c.s+i+1))[0];
    if (n+l+2>m) {m=n+l+2; m<<=1; s = realloc(s, m);}
    for (j=0; j<l; ++j) s[n++] = c.s[i];
  }
  expanded.s = (uint8_t*) s;
  expanded.n = n;
  expanded.compressed = 0;
  expanded.fmt = '1';
  expanded.unit = 1;
  return expanded;
}


#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "cdata.h"

cdata_t* fmt1_read_raw(char *fname, int verbose) {

  gzFile fh = wzopen(fname, 1);
  char *line = NULL;
  uint64_t n = 0, m=1<<22;
  uint8_t *s = calloc(m, 1);
  while (gzFile_read_line(fh, &line) > 0) {
    s[n++] = line[0];
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
  c->fmt = '1';
  return c;
}

/* compressed:
   3byte --- value (1byte) + run len (2bytes)
   value is unrestricted ASCII
 */
void fmt1_compress(cdata_t *c) {
  uint64_t n = 0;
  uint8_t *s = NULL;
  uint64_t i = 0; uint16_t l = 0; uint8_t u0 = 0;
  for (i=0, l=0; i<c->n; ++i) {
    /* either not the same as before or reach block size max */
    if ((l != 0 && c->s[i] != u0) || l+2 >= 1<<15) {
      s = realloc(s, n+3);
      s[n] = u0;
      *((uint16_t*) (s+n+1)) = l;
      n += 3;
      l = 1;
    } else {
      ++l;
    }
    u0 = c->s[i];
  }
  /* the last rle */
  s = realloc(s, n+3);
  s[n] = u0;
  *((uint16_t*) (s+n+1)) = l;
  n += 3;
  
  free(c->s);
  c->s = s;
  c->n = n;
  c->compressed = 1;
}

void fmta_tryBinary2byteRLE_ifsmaller(cdata_t *c) {
  uint64_t n = 0;
  uint8_t *s = NULL;
  uint64_t i=0; uint16_t l=0; uint8_t u0=0;
  for (i=0, l=0; i<c->n; ++i) {
    uint8_t u = (c->s[i>>3]>>(i&0x7))&0x1;
    /* either not the same as before or reach block size max */
    if ((l != 0 && u != u0) || l+2 >= 1<<15) {
      s = realloc(s, n+3);
      s[n] = u0+'0';
      *((uint16_t*) (s+n+1)) = l;
      n += 3;
      l = 1;
    } else {
      ++l;
    }
    u0 = u;
  }
  /* the last rle */
  s = realloc(s, n+3);
  s[n] = u0+'0';
  *((uint16_t*) (s+n+1)) = l;
  n += 3;
  if (c->n>>3 > n) {
    free(c->s);
    c->s = s;
    c->n = n;
    c->fmt = '1';
    c->compressed = 1;
  }
}


