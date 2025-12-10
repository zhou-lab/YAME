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
 * Compressed layout
 * -----------------
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
    fprintf(stderr, "[%s:%d] Data of length %llu loaded\n", __func__, __LINE__, n);
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
cdata_t fmt6_decompress(cdata_t c) {
  cdata_t expanded = {0};
  expanded.s = malloc(cdata_nbytes(&c));
  memcpy(expanded.s, c.s, cdata_nbytes(&c));
  expanded.n = c.n;
  expanded.compressed = 0;
  expanded.fmt = '6';
  expanded.unit = 2;
  return expanded;
}
