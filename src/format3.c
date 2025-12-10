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

/**
 * Format 3 (M/U counts) storage layout
 *
 * Inflated (uncompressed) representation
 * --------------------------------------
 * - c->compressed == 0
 * - c->n = number of CpG sites
 * - Each site uses c->unit bytes (little-endian) at s + i*unit:
 *
 *       [ M | U ]  packed into 'unit' bytes
 *
 *   Default: unit == 8
 *       uint64_t mu = (M << 32) | U
 *       bytes 0..3: U
 *       bytes 4..7: M
 *
 *   f3_get_mu(c,i)  -> uint64  [ M(32b) | U(32b) ]
 *   f3_set_mu(c,i,M,U) packs M,U back into unit bytes.
 *
 *
 * Compressed representation
 * -------------------------
 * - c->compressed == 1
 * - c->s is a byte stream of variable-length records.
 * - Each record type is identified by its low 2 bits:
 *
 *  2byte | U=M=0 -------------- = run len (14 bit) + 0 (2bit)
 *  1byte | U,M in [0,7] ------- = M (3bit)  | U (3bit)  + 1 (2bit)
 *  2byte | U,M in [0,127]------ = M (7bit)  | U (7bit)  + 2 (2bit)
 *  8byte | M,U in [128,2**31]-- = M (31bit) | U (31bit) + 3 (2bit)
 *
 *   Type 0: zero-run      (low bits 00)
 *     - 2 bytes (uint16, little-endian)
 *     - layout: [ run_len(14 bits) | 00 ]
 *     - represents 'run_len' consecutive sites with M=0, U=0.
 *
 *   Type 1: small M,U     (low bits 01)
 *     - 1 byte
 *     - layout: [ M(3 bits) | U(3 bits) | 01 ]
 *     - valid ranges: 0 <= M,U <= 6
 *
 *   Type 2: medium M,U    (low bits 10)
 *     - 2 bytes (uint16, little-endian)
 *     - layout: [ M(7 bits) | U(7 bits) | 10 ]
 *     - valid ranges: 0 <= M,U <= 126
 *
 *   Type 3: large M,U     (low bits 11)
 *     - 8 bytes (uint64, little-endian)
 *     - layout: [ M(31 bits) | U(31 bits) | 11 ]
 *     - M,U may be right-shifted (fitMU) to fit 31 bits.
 *
 *
 * Schematic before/after compression
 * ----------------------------------
 *
 *   Inflated view (per-site array):
 *     index i:   [ M_i, U_i ]  (stored as unit-byte packed MU)
 *     c->n = number of sites
 *
 *   Compressed view (byte stream):
 *     [ rec0 ][ rec1 ][ rec2 ] ... [ recK ]
 *
 *       where each rec* is one of:
 *         - zero-run (covers many sites with M=U=0), or
 *         - a single non-zero (M,U) entry
 *
 *   Decompression walks the compressed stream, expanding:
 *     - Type 0 into 'run_len' copies of (M=0,U=0)
 *     - Types 1/2/3 into one (M,U) pair each,
 *   and writes them back into the inflated array using f3_pack_mu().
 */

static int is_int(char *s) {
  size_t i;
  for (i=0; i<strlen(s); ++i) {
    if (!isdigit(s[i])) return 0;
  }
  return 1;
}

static int fitMU(uint64_t *M, uint64_t *U, uint64_t nbits) {
  int modified = 0;
  uint64_t max = (1ul << nbits)-1;
  while ((*M) >= max || (*U) >= max) {
    (*M) >>= 1;
    (*U) >>= 1;
    modified = 1;
  }
  return modified;
}

static void pack_value(uint8_t *data, uint64_t value, uint8_t unit) {
  for (uint8_t i=0; i<unit; ++i) {
    data[i] = (value & 0xff);
    value >>= 8;
  }
}

static uint64_t unpack_value(uint8_t *data, uint8_t unit) {
  uint64_t v = 0;
  for (uint8_t i=0; i<unit; ++i) v |= ((uint64_t) data[i] << (8*i));
  return v;
}

// TODO: add fitMU here to be safe
static void f3_pack_mu(uint8_t *data, uint64_t M, uint64_t U, uint8_t unit) {
  if (!unit) {
    fprintf(stderr, "[%s:%d] Unknown unit size.\n", __func__, __LINE__);
    fflush(stderr);
    exit(1);
  }
  pack_value(data, (M<<(unit*4)) | U, unit);
}

void f3_set_mu(cdata_t *c, uint64_t i, uint64_t M, uint64_t U) {
  f3_pack_mu(c->s+i*c->unit, M, U, c->unit);
}

// note this function generate uin32_t not uint64_t. Please fix.
uint64_t f3_get_mu(cdata_t *c, uint64_t i) {
  uint8_t *data = c->s + (c->unit*i);
  uint64_t mu = 0;
  for (uint8_t j=0; j<c->unit; ++j) {
    mu |= (((uint64_t) data[j]) << (8*j));
  }
  return (mu>>(c->unit*4)<<32) | (mu & ((1ul<<(c->unit*4))-1));
}

/* uncompressed: [ M (uint32_t) | U (uint32_t) ] */
/* usually unit == 8 by default for minimal loss */
cdata_t* fmt3_read_raw(char *fname, uint8_t unit, int verbose) {
  if (!unit) unit = 8;
  gzFile fh = wzopen(fname, 1);
  char *line = NULL;
  uint8_t *s = NULL; uint64_t n = 0;
  char **fields; int nfields;
  while (gzFile_read_line(fh, &line) > 0) {
    line_get_fields(line, "\t", &fields, &nfields);
    if (nfields < 2) wzfatal("Number of fields <2. Abort.");
    if (!is_int(fields[0]) || !is_int(fields[1]))
      wzfatal("Field 1 or 2 is not a nonnegative integer.");
    uint64_t M = atol(fields[0]);
    uint64_t U = atol(fields[1]);
    s = realloc(s, (n+1)*unit);
    fitMU(&M, &U, unit*4);
    f3_pack_mu(s+n*unit, M, U, unit);
    n++;
    free_fields(fields, nfields);
  }
  free(line);
  wzclose(fh);
  if (verbose) {
    fprintf(stderr, "[%s:%d] Vector of length %llu loaded\n", __func__, __LINE__, n);
    fflush(stderr);
  }
  cdata_t *c = calloc(sizeof(cdata_t),1);
  c->s = s;
  c->n = n;
  c->compressed = 0;
  c->fmt = '3';
  c->unit = unit;
  return c;
}

void fmt3_compress(cdata_t *c) {
  uint8_t *s = NULL;
  uint64_t n = 0;
  uint64_t i = 0;
  uint64_t l = 0;
  for (i=0; i<c->n; i++) {
    uint64_t MU = f3_get_mu(c, i);
    uint64_t M = MU>>32;
    uint64_t U = MU<<32>>32;
    if (M>0 || U>0 || l+2 >= (1ul<<14)) {
      if (l>0) {
        s = realloc(s, n+2);
        pack_value(s+n, l<<2, 2);
        n += 2;
        if (M>0 || U>0) l = 0;
        else l = 1;
      }
      if (M>0 || U>0) {
        if (M<7 && U<7) {
          s = realloc(s, n+1);
          s[n] = (M<<5) | (U<<2) | 0x1;
          n++;
        } else if (M<127 && U<127) {
          s = realloc(s, n+2);
          pack_value(s+n, (M<<9) | (U<<2) | 0x2, 2);
          n += 2;
        } else {
          fitMU(&M, &U, 31);
          s = realloc(s, n+8);
          pack_value(s+n, (M<<33) | (U<<2) | 3ul, 8);
          n += 8;
        }
      }
    } else {
      ++l;
    }
  }
  if (l>0) {
    s = realloc(s, n+2);
    *((uint16_t*) (s+n)) = (uint16_t) l<<2;
    n += 2;
  }
  free(c->s);
  c->s = s;
  c->n = n;
  c->compressed = 1;
}

static uint64_t get_data_length(cdata_t *c, uint8_t *unit) {
  uint8_t nbits = 1; // half unit nbits, M or U.
  uint64_t n = 0;
  for (uint64_t i=0; i < c->n; ) {
    if ((c->s[i] & 0x3) == 0) {
      n += (((uint16_t*) (c->s+i))[0])>>2;
      i += 2;
    } else if ((c->s[i] & 0x3) == 1) {
      uint64_t M = (c->s[i])>>5;
      uint64_t U = ((c->s[i])>>2) & 0x7;
      while (M >= (1ul << nbits) || U >= (1ul << nbits)) nbits++;
      n++; i++;
    } else if ((c->s[i] & 0x3) == 2) {
      uint64_t M = (((uint16_t*) (c->s+i))[0])>>2;
      uint64_t U = M & ((1<<7)-1);
      M >>= 7;
      while (M >= (1ul << nbits) || U >= (1ul << nbits)) nbits++;
      n++; i += 2;
    } else {
      uint64_t M = (((uint64_t*) (c->s+i))[0])>>2;
      uint64_t U = M & ((1ul<<31)-1);
      M >>= 31;
      while (M >= (1ul << nbits) || U >= (1ul << nbits)) nbits++;
      n++; i += 8;
    }
  }
  *unit = ((nbits+3)>>2); // nbits*2/8
  return n;
}

cdata_t fmt3_decompress(cdata_t c) {
  uint8_t unit = 1;
  uint64_t n0 = get_data_length(&c, &unit);
  cdata_t inflated = {0};
  if (c.unit) inflated.unit = c.unit;
  else inflated.unit = unit; // use inferred max unit if unset
  uint8_t *s = calloc(inflated.unit*n0, sizeof(uint8_t));
  uint64_t n = 0;
  for (uint64_t i=0; i < c.n; ) {
    if ((c.s[i] & 0x3) == 0) {
      uint64_t l = unpack_value(c.s+i, 2)>>2; // the length is 14 bits, so unit = 2
      memset(s+n*inflated.unit, 0, inflated.unit*l); n += l;
      i += 2;
    } else if ((c.s[i] & 0x3) == 1) {
      uint64_t M = (c.s[i])>>5;
      uint64_t U = ((c.s[i])>>2) & 0x7;
      if (inflated.unit == 1) fitMU(&M, &U, 4);
      else fitMU(&M, &U, inflated.unit<<2);
      f3_pack_mu(s+(n++)*inflated.unit, M, U, inflated.unit);
      i++;
    } else if ((c.s[i] & 0x3) == 2) {
      uint64_t M = unpack_value(c.s+i, 2)>>2;
      uint64_t U = M & ((1ul<<7)-1);
      M >>= 7;
      if (inflated.unit == 1) fitMU(&M, &U, 4);
      else fitMU(&M, &U, inflated.unit<<2);
      f3_pack_mu(s+(n++)*inflated.unit, M, U, inflated.unit);
      i += 2;
    } else {
      uint64_t M = unpack_value(c.s+i, 8)>>2;
      uint64_t U = M & ((1ul<<31)-1);
      M >>= 31;
      if (inflated.unit == 1) fitMU(&M, &U, 4);
      else fitMU(&M, &U, inflated.unit<<2);
      f3_pack_mu(s+(n++)*inflated.unit, M, U, inflated.unit);
      i += 8;
    }
  }
  inflated.s = s;
  inflated.n = n;
  inflated.compressed = 0;
  inflated.fmt = '3';
  return inflated;
}
