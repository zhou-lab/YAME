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

#ifndef _CDATA_H
#define _CDATA_H

#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <limits.h>
#include <inttypes.h>
#include "khash.h"
#include "wzmisc.h"
#include "wzmisc.h"
#include "wzbed.h"
#include "bgzf.h"

#define CDSIG 266563789635

/**
 * cdata_t
 * --------
 * Core container for all YAME data formats (0–7).  A cdata_t may hold either
 * a compressed byte stream (on-disk encoding) or an in-memory, decompressed /
 * indexed representation depending on `compressed`.
 *
 * Fields:
 *
 *   s : Pointer to the raw byte buffer.  Interpretation depends on `fmt` and
 *       whether the content is compressed or decompressed.
 *
 *   n : Size of `s`.
 *       • If compressed  (compressed = 1):  n = number of bytes in the on-disk
 *         representation.  (Exception: fmt0 stores bit-packed vectors; see below.)
 *       • If uncompressed (compressed = 0): n = number of data units (typically
 *         number of CpGs / rows), excluding key bytes.
 *
 *   compressed :
 *       1 → `s` is a compressed stream (formats 0–7 as written on disk).
 *       0 → `s` is an uncompressed / indexed representation (formats 3,6,7
 *           after decompression, where each logical unit has fixed width).
 *
 *   fmt :
 *       One-character format code: '0'..'7'.
 *       Determines the encoding, interpretation of `s`, and the meaning of `unit`.
 *
 *   unit :
 *       Size in bytes of a *single decompressed data unit*.
 *
 *       • fmt0, fmt1, fmt6 : bit-packed encodings → unit = 0
 *       • fmt2, fmt3, fmt4, fmt5 : unit = byte-width of each value (1,2,4,8)
 *       • fmt7 (decompressed via fmt7_decompress): unit = 8 (fixed 8-byte entries)
 *
 *       When compressed = 1, `unit` describes the *intended* size of each
 *       decompressed item but does not apply to the raw bytes on disk.
 *
 *   aux :
 *       Optional auxiliary pointer, format-specific.  Set on-demand.
 *       Examples:
 *         • fmt2: f2_aux_t (keys[] and pointer to state data)
 *         • fmt3: counters for M/U decoding
 *         • fmt6: universe bitmask accessors
 *         • fmt7: row_reader_t for streamed BED iteration
 *
 *
 * Special notes:
 *   • Format 0: n is the number of bits, stored bit-packed in s[].
 *   • Format 1: run-length-encoded integer stream; unit=0 and inflation is
 *               performed by format-specific helpers.
 *   • Format 7: only the *uncompressed* version has fixed-width entries;
 *               the compressed on-disk version uses delta-encoding.
 *
 * cdata_t is intentionally minimal; decoding, indexing, and slicing are
 * implemented in the per-format helpers.
 */
typedef struct cdata_t {
  uint8_t *s;       /* byte buffer */
  uint64_t n;       /* length in bytes (compressed) or #units (uncompressed) */
  int compressed;   /* 1=compressed stream, 0=indexed/decompressed */
  char fmt;         /* format code '0'..'7' */
  uint8_t unit;     /* size of each decompressed unit (0 for bit-packed fmt0/1/6) */
  void *aux;        /* optional per-format auxiliary structure */
} cdata_t;

static inline uint64_t cdata_nbytes(const cdata_t *c) {
  uint64_t n = 0;
  switch(c->fmt) {
  case '0': n = ((c->n+7)>>3); break;
  case '6': n = ((c->n+3)>>2); break;
  default: n = c->n;
  }

  if (!c->compressed) {
    if(c->fmt == '3') {
      n *= c->unit;
    } // TODO: add other formats
  }
  return n;
}

typedef struct f2_aux_t {
  uint64_t nk;                  // num keys
  char **keys;                  // pointer to keys, doesn't own memory
  uint8_t *data;                // pointer to data, doesn't own memory
} f2_aux_t;

static inline void free_cdata(cdata_t *c) {
  if (c->s) free(c->s);
  if (c->fmt == '2' && c->aux) {
    free(((f2_aux_t*) c->aux)->keys);
    free(c->aux);
  }
  if (c->fmt == '7' && c->aux) free(c->aux);
  c->s = NULL;
}

static inline size_t bit_count(cdata_t c) {

  /* create a look-up table */
  int byte2cnt[256]; int p;
  for (p=0; p<256; ++p) {
    unsigned char q = p; int ii, cnt = 0;
    for (ii=0; ii<8; ++ii) { if (q&1) cnt++; q>>=1; }
    byte2cnt[p] = cnt;
  }
  
  size_t i,k,m = 0;
  for (i=0; i<(c.n>>3); ++i) m += byte2cnt[c.s[i]];   // full bytes
  for (k=0; k<(c.n&0x7); ++k) m += (c.s[i]>>k) & 0x1; // last byte
  return m;
}

void cdata_compress(cdata_t *c);
cdata_t decompress(cdata_t c);
void decompress_in_situ(cdata_t *c);

static inline uint64_t cdata_n(cdata_t *c) {
  if (!c->compressed) return c->n;
  cdata_t c2 = decompress(*c);
  uint64_t n = c2.n;
  free(c2.s);
  return n;
}

void convertToFmt0(cdata_t *c);
#define FMT0_IN_SET(c, i) ((c).s[(i)>>3] & (1u<<((i)&0x7)))
#define FMT0_SET(c, i) (c.s[(i)>>3] |= (1u<<((i)&0x7)))

#define _FMT0_IN_SET(s, i) ((s)[(i)>>3] & (1u<<((i)&0x7)))
#define _FMT0_SET(s, i) ((s)[(i)>>3] |= (1u<<((i)&0x7)))

void fmt2_set_aux(cdata_t *c);
uint8_t* fmt2_get_data(const cdata_t *c);
uint64_t fmt2_get_keys_n(const cdata_t *c);
uint64_t fmt2_get_keys_nbytes(const cdata_t *c);
uint64_t f2_get_uint64(cdata_t *c, uint64_t i);
char* f2_get_string(cdata_t *c, uint64_t i);

void f3_set_mu(cdata_t *c, uint64_t i, uint64_t M, uint64_t U);
uint64_t f3_get_mu(cdata_t *c, uint64_t i);
#define MU2beta(mu) (double) ((mu)>>32) / (((mu)>>32) + ((mu)&0xffffffff))
#define MU2cov(mu) (((mu)>>32) + ((mu)&0xffffffff))

// fmt6 as a quaternary
#define FMT6_2BIT(c, i) (((c).s[i>>2]>>((i&0x3)*2)) & 0x3)
#define FMT6_IN_SET(c, i) ((c).s[i>>2] & (1<<((i&0x3)*2)))
#define FMT6_IN_UNI(c, i) ((c).s[i>>2] & (1<<((i&0x3)*2+1)))
#define FMT6_SET0(c, i) ((c).s[i>>2] = ((c).s[i>>2] & ~(3<<((i&0x3)*2))) | (2<<((i&0x3)*2))) // 10
#define FMT6_SET1(c, i) ((c).s[i>>2] |= (3<<((i&0x3)*2))) // 11
#define FMT6_SET_NA(c, i) ((c).s[i>>2] &= (~(3<<((i&0x3)*2)))) // 00

/* this doesn't work for format 2, no copy of aux */
static inline cdata_t cdata_duplicate(cdata_t c) {
  cdata_t cout = c;
  cout.s = (uint8_t*) malloc(cdata_nbytes(&c));
  if (cout.s==NULL) wzfatal("[cdata_duplicate] Cannot allocate memory.\n");
  memcpy(cout.s, c.s, cdata_nbytes(&c));
  return cout;
}

uint64_t fmt7_data_length(const cdata_t *c);
cdata_t fmt7_sliceToBlock(cdata_t *cr, uint64_t beg, uint64_t end);
cdata_t fmt7_sliceToIndices(cdata_t *cr, int64_t *row_indices, int64_t n_indices);
cdata_t fmt7_sliceToMask(cdata_t *cr, cdata_t *c_mask);

static inline void slice(cdata_t *c, uint64_t beg, uint64_t end, cdata_t *c_sliced) {

  if (c->compressed) {
    fprintf(stderr, "[%s:%d] Cannot slice compressed data.\n", __func__, __LINE__);
    fflush(stderr);
    exit(1);
  }
  if (end > c->n-1) end = c->n-1;
  if (end < beg) wzfatal("Slicing negative span.");

  c_sliced->s = realloc(c_sliced->s, (end-beg+1)*c->unit);
  memcpy(c_sliced->s, c->s+beg*c->unit, (end-beg+1)*c->unit);
  c_sliced->n = end-beg+1;
  c_sliced->compressed = 0;
  c_sliced->fmt = c->fmt;
  c_sliced->unit = c->unit;
}

/**
 * @brief A streaming cursor for iterating over row-coordinate records.
 *
 * This struct represents the current position in a row-coordinate cdata_t
 * (format 7, usually a `.cr` file). It does not own any memory; all pointers
 * refer into the underlying cdata_t::s buffer.
 *
 * Fields:
 *   - index: 1-based row index (i.e., which row in the matrix this record
 *            corresponds to). Incremented by row_reader_next_loc().
 *   - chrm : pointer to the chromosome name string inside cdata_t::s.
 *   - loc  : byte offset into cdata_t::s where this record is stored.
 *   - value: genomic coordinate (typically 1-based CpG position, "beg1").
 *
 * Typical usage:
 *   row_reader_t rdr = {0};
 *   while (row_reader_next_loc(&rdr, cr)) {
 *       // use rdr.chrm, rdr.value, rdr.index, ...
 *   }
 *
 * The helper function row_reader_next_loc() advances this cursor to the next
 * record and fills these fields accordingly.
 */
typedef struct row_reader_t {
  uint64_t index;
  char *chrm;                   // on cdata_t.s
  uint64_t loc;                 // on cdata_t.s
  uint64_t value;
} row_reader_t;

KHASH_MAP_INIT_STR(str2int, uint64_t) // Initialize a hashmap with keys as strings and values as uint64_t

/**
 * @brief Per-chromosome coarse index into a row-coordinate track.
 *
 * For each chromosome we construct a coarse-grained index in blocks of
 * fixed genomic size (2^17 bp). Each entry k in the arrays corresponds to
 * the first CpG at genomic coordinate >= (k << 17). locs map these block
 * starts to the positions in the cdata_t::s.
 *
 * Fields:
 *   - locs[k]: byte offset into the original cdata_t::s buffer at which the
 *              first record for block k (or the next non-empty block) starts.
 *   - vals[k]: genomic coordinate (1-based) of that record.
 *   - inds[k]: 1-based row index corresponding to that record.
 *   - n     : number of blocks currently stored (size of the arrays).
 *
 * The arrays grow as we encounter larger coordinates during init_finder().
 * These arrays do not own the underlying chromosome strings; they only store
 * offsets and coordinates.
 */
typedef struct chromosome_t {
  uint64_t *locs;
  uint64_t *vals;
  uint64_t *inds;
  uint64_t n;
} chromosome_t;

/**
 * @brief Global index for fast row lookup by chromosome and coordinate.
 *
 * row_finder_t is an in-memory index built from a row-coordinate cdata_t
 * (typically a format 7 `.cr` file). It supports O(1) chromosome lookup and
 * O(1) jump to a coarse genomic bin, followed by a short linear scan.
 *
 * Fields:
 *   - chrms: array of chromosome_t structures, one per chromosome seen in
 *            the coordinate file.
 *   - n    : number of chromosomes.
 *   - h    : khash mapping chromosome name (char*) to integer index
 *            into the chrms[] array.
 *
 * Construction:
 *   row_finder_t fdr = init_finder(cr);
 *
 * Query:
 *   uint64_t row = row_finder_search("chr1", 1234567, &fdr, cr);
 *   if (row == 0) {} // not found
 *
 * Memory:
 *   Use free_row_finder(&fdr) to release all internal allocations.
 */
typedef struct row_finder_t {
  chromosome_t *chrms;
  int n; // number of chromosomes
  khash_t(str2int) *h; // chromosome string > chromosome_t
} row_finder_t;

static inline void free_row_finder(row_finder_t *fdr) {
  for (int i=0; i<fdr->n; ++i) {
    free(fdr->chrms[i].locs);
    free(fdr->chrms[i].vals);
    free(fdr->chrms[i].inds);
  }
  free(fdr->chrms);
  kh_destroy(str2int, fdr->h);
}

int row_reader_next_loc(row_reader_t *rdr, const cdata_t *c);
row_finder_t init_finder(cdata_t *cr);
uint64_t row_finder_search(char *chrm, uint64_t beg1, row_finder_t *fdr, cdata_t *cr);

#endif /* _CDATA_H */
