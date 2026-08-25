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
 *       Optional auxiliary pointer, format-specific.  Set on-demand, and by
 *       exactly two formats -- the other two entries this list used to carry
 *       (fmt3 M/U counters, fmt6 universe accessors) were never allocated by
 *       anything, which made free_cdata() look incomplete when it is not:
 *         • fmt2: f2_aux_t (keys[] and pointer to state data), fmt2_set_aux()
 *         • fmt7: row_reader_t for streamed BED iteration, fmt7_next_bed()
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

  /* Decompressed, `n` counts ROWS, so a fixed-width format has to multiply by
   * its unit. Naming only format 3 here left formats 1, 4 and 5 reporting a
   * row count as a byte count -- a quarter of the truth for format 4, whose
   * unit is 4 -- which is a buffer size everywhere this is used to allocate
   * or copy. Formats 0 and 6 are bit-packed and already sized by the switch
   * above; 2 and 7 are not flat row vectors and have no answer here. */
  if (!c->compressed) {
    switch (c->fmt) {
    case '0': case '6': break;            /* bit-packed, sized above */
    case '2': case '7': break;            /* not a flat n-by-unit array */
    default: n *= c->unit ? c->unit : 1;  /* 1, 3, 4, 5 */
    }
  }
  return n;
}

typedef struct f2_aux_t {
  uint64_t nk;                  // num keys
  char **keys;                  // pointer to keys, doesn't own memory
  uint8_t *data;                // pointer to data, doesn't own memory
} f2_aux_t;

/**
 * Release a record, and leave it safe to release again.
 *
 * The `if (c->s)` guard made this look idempotent, and for `s` it was. `aux`
 * was freed and left pointed at, so a second call -- an error path added
 * above an existing one, say -- double-freed a format 2 key table. Clearing
 * the counts too means a freed record cannot go on claiming rows it no
 * longer holds; nothing reads a field after freeing today, and now nothing
 * can start to by accident.
 */
static inline void free_cdata(cdata_t *c) {
  if (c->s) free(c->s);
  if (c->fmt == '2' && c->aux) {
    free(((f2_aux_t*) c->aux)->keys);
    free(c->aux);
  }
  if (c->fmt == '7' && c->aux) free(c->aux);
  c->s = NULL;
  c->aux = NULL;
  c->n = 0;
  c->compressed = 0;
}

/* Set bits per byte value. Built once at load rather than rebuilt on every
 * bit_count() call, which used to cost 2048 iterations before counting a
 * single bit. */
static const uint8_t byte2cnt[256] = {
  0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4, 1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
  1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5, 2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5, 2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6, 3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
  1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5, 2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6, 3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6, 3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
  3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7, 4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8
};

static inline size_t bit_count(cdata_t c) {
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
/* Bit accessors. Every argument is parenthesized and `i` is evaluated more
 * than once -- pass a plain variable, never `idx++` or a call. */
#define FMT0_IN_SET(c, i) ((c).s[(i)>>3] & (1u<<((i)&0x7)))
#define FMT0_SET(c, i) ((c).s[(i)>>3] |= (1u<<((i)&0x7)))

#define _FMT0_IN_SET(s, i) ((s)[(i)>>3] & (1u<<((i)&0x7)))
#define _FMT0_SET(s, i) ((s)[(i)>>3] |= (1u<<((i)&0x7)))

void fmt2_set_aux(cdata_t *c);
uint8_t* fmt2_get_data(const cdata_t *c);
uint64_t fmt2_get_keys_n(const cdata_t *c);
uint64_t fmt2_get_keys_nbytes(const cdata_t *c);
uint64_t f2_get_uint64(cdata_t *c, uint64_t i);
char* f2_get_string(cdata_t *c, uint64_t i);

void     f3_set_mu(cdata_t *c, uint64_t i, uint64_t M, uint64_t U);
uint64_t f3_get_mu(cdata_t *c, uint64_t i);
#define MU2beta(mu) (double) ((mu)>>32) / (((mu)>>32) + ((mu)&0xffffffff))
#define MU2cov(mu) (((mu)>>32) + ((mu)&0xffffffff))

/* fmt6 as a quaternary. Same rules as the fmt0 accessors above: `i` is
 * parenthesized, and it is evaluated up to four times (FMT6_SET0), so it
 * must be a plain variable. The bare `i` these used to carry meant
 * FMT6_IN_SET(c, i+1) read c.s[i + (1>>2)] -- that is, c.s[i], the wrong
 * byte, silently, because >> binds looser than +. */
#define FMT6_2BIT(c, i) (((c).s[(i)>>2]>>(((i)&0x3)*2)) & 0x3)
#define FMT6_IN_SET(c, i) ((c).s[(i)>>2] & (1u<<(((i)&0x3)*2)))
#define FMT6_IN_UNI(c, i) ((c).s[(i)>>2] & (1u<<(((i)&0x3)*2+1)))
#define FMT6_SET0(c, i) ((c).s[(i)>>2] = ((c).s[(i)>>2] & ~(3<<(((i)&0x3)*2))) | (2<<(((i)&0x3)*2))) // 10
#define FMT6_SET1(c, i) ((c).s[(i)>>2] |= (3<<(((i)&0x3)*2))) // 11
#define FMT6_SET_NA(c, i) ((c).s[(i)>>2] &= (~(3<<(((i)&0x3)*2)))) // 00

/**
 * Copy a record's bytes.
 *
 * Formats 2 and 7 are refused rather than half-copied. A format 2 record's
 * aux holds pointers INTO the buffer being duplicated, so the copy came back
 * sharing the original's key table -- fine until either one is freed. Format
 * 7 has no flat row array to copy at all. This used to be a comment saying
 * the first one does not work; a comment does not stop the call.
 */
static inline cdata_t cdata_duplicate(cdata_t c) {
  if (c.fmt == '2' || c.fmt == '7')
    wzfatal("[cdata_duplicate] format %c is not a flat row vector; copy it "
            "with its own helper.\n", c.fmt);
  cdata_t cout = c;
  cout.aux = NULL;              /* the copy owns no auxiliary structure */
  uint64_t nb = cdata_nbytes(&c);
  cout.s = (uint8_t*) malloc(nb ? nb : 1);
  if (cout.s==NULL) wzfatal("[cdata_duplicate] Cannot allocate memory.\n");
  memcpy(cout.s, c.s, nb);
  return cout;
}

uint64_t fmt7_data_length(const cdata_t *c);
cdata_t fmt7_sliceToBlock(cdata_t *cr, uint64_t beg, uint64_t end);
cdata_t fmt7_sliceToIndices(cdata_t *cr, int64_t *row_indices, int64_t n_indices);
cdata_t fmt7_sliceToMask(cdata_t *cr, cdata_t *c_mask);

/**
 * Copy rows [beg,end] out of a decompressed record.
 *
 * `unit` is bytes per row for formats 1/3/4, and a flat memcpy is right for
 * them. It is not right for everything: format 0 packs eight rows into a byte
 * and format 6 packs four, so there `beg` is a bit offset that byte
 * arithmetic reads as a byte offset -- silently returning the wrong rows, or
 * running off the end of the buffer. Those two are unpacked and repacked
 * here instead.
 *
 * Formats 2 and 7 are refused rather than guessed at: a format 2 record
 * carries a key table ahead of its rows and a format 7 record is a delta
 * stream with no fixed row width, so neither survives a positional copy.
 * `yame rowsub` implements the structure-aware versions.
 */
static inline void slice(cdata_t *c, uint64_t beg, uint64_t end, cdata_t *c_sliced) {

  if (c->compressed) {
    fprintf(stderr, "[%s:%d] Cannot slice compressed data.\n", __func__, __LINE__);
    fflush(stderr);
    exit(1);
  }
  if (end > c->n-1) end = c->n-1;
  if (end < beg) wzfatal("Slicing negative span.");

  uint64_t n_out = end - beg + 1;
  c_sliced->n = n_out;
  c_sliced->compressed = 0;
  c_sliced->fmt = c->fmt;
  c_sliced->unit = c->unit;

  if (c->fmt == '2' || c->fmt == '7') {
    wzfatal("[%s:%d] Cannot slice format %c by row block: it is not a flat "
            "row vector. Use `yame rowsub` instead.\n",
            __func__, __LINE__, c->fmt);
  } else if (c->fmt == '0') {            /* 1 bit per row */
    uint64_t nb = (n_out + 7) >> 3;
    c_sliced->s = realloc(c_sliced->s, nb);
    memset(c_sliced->s, 0, nb);
    for (uint64_t i = 0; i < n_out; ++i) {
      uint64_t si = beg + i;
      if (c->s[si >> 3] & (1u << (si & 0x7)))
        c_sliced->s[i >> 3] |= (uint8_t)(1u << (i & 0x7));
    }
  } else if (c->fmt == '6') {            /* 2 bits per row */
    uint64_t nb = (n_out + 3) >> 2;
    c_sliced->s = realloc(c_sliced->s, nb);
    memset(c_sliced->s, 0, nb);
    for (uint64_t i = 0; i < n_out; ++i) {
      uint64_t si = beg + i;
      uint8_t v = (c->s[si >> 2] >> ((si & 0x3) * 2)) & 0x3;
      c_sliced->s[i >> 2] |= (uint8_t)(v << ((i & 0x3) * 2));
    }
  } else {
    c_sliced->s = realloc(c_sliced->s, n_out * c->unit);
    memcpy(c_sliced->s, c->s + beg * c->unit, n_out * c->unit);
  }
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
