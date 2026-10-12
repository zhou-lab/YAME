// SPDX-License-Identifier: LicenseRef-CHOP-Academic-BSD-2-Clause
/**
 * This file is part of YAME.
 *
 * Copyright (C) 2021-present The Children's Hospital of Philadelphia
 *
 * Use of this software is available to academic and non-profit institutions
 * for research purposes under the 2-Clause BSD License; for use or transfers
 * to commercial entities, inquire with Dr. Wanding Zhou at zhouw3@chop.edu.
 * See the LICENSE file at the root of the repository for the full terms.
 */

/**
 * rdx.h -- the row index, kept in the .idx as a third column.
 *
 * A format-3 record is a stream of variable-length tokens, so the BGZF block
 * holding row N is unknown without reading every token before it. The row
 * index records, for each format-3 record, one checkpoint per BGZF block: the
 * row and the virtual offset of the first token that starts in that block.
 * Each token is self-contained, so decoding can begin at any checkpoint, and
 * a row read inflates only the blocks that hold the rows asked for.
 *
 * The .cx bytes are unchanged; a record with no column, or one the column no
 * longer matches, is read in full exactly as before. A record is matched on
 * its offset, its format and its stream length, and every checkpoint used is
 * checked against the byte it names.
 *
 * Private to the YAME binary: no public header includes this one.
 */

#ifndef _YAME_RDX_H
#define _YAME_RDX_H

#include <stdint.h>
#include "bgzf.h"
#include <stdio.h>
#include "cdata.h"
#include "index.h"

typedef struct {
  uint64_t row0;   /* first row of the token */
  uint64_t spos;   /* byte position of the token in the record's stream */
  int64_t  voff;   /* BGZF virtual offset of that byte */
  uint8_t  byte0;  /* the token's first byte, checked on use */
} rdx_ckpt_t;

typedef struct {
  int64_t  rec_voff;  /* the record's offset, as the .idx gives it */
  int64_t  end_voff;  /* where the next record's read begins; -1 unknown */
  uint64_t nbytes;    /* the stream length the record header carries */
  uint64_t nrows;
  uint64_t nckpt;
  rdx_ckpt_t *ckpt;
} rdx_rec_t;

typedef struct {
  uint64_t nrec, mrec;
  rdx_rec_t *rec;     /* sorted by rec_voff once loaded */
  /* builder only: the file's BGZF blocks */
  uint64_t nblk;
  int64_t  *blk_coff;
  uint64_t *blk_upos; /* uncompressed position of the block's first byte */
  uint32_t *blk_usize;
} rdx_t;

/* Builder. rdx_open_build() walks the BGZF block headers of fname (no
 * inflate) and returns NULL if the file is not plain BGZF. rdx_add() takes a
 * record as read_cdata2() returned it, still compressed, with the offset the
 * read began at; records that are not format 3 are ignored. */
rdx_t *rdx_open_build(const char *fname);
void   rdx_add(rdx_t *x, int64_t rec_voff, const cdata_t *c);

/* A record copied byte-for-byte to another file: shift every compressed
 * address by coff_delta and file it under new_voff, with the reader landing
 * at new_end_voff after it. */
void   rdx_add_moved(rdx_t *x, const rdx_rec_t *r, int64_t new_voff,
                     int64_t coff_delta, int64_t new_end_voff);

/* Write an .idx: name, offset and, when x has any record, the row index as
 * a third column on every line ("." for a record it does not cover). With
 * no row index this is writeIndex(), byte for byte. */
void   rdx_write_index(FILE *fp, index_t *idx, rdx_t *x);

/* The row index in <fname_cx>.idx, parsed now; NULL when there is none. */
rdx_t *rdx_load_index(const char *fname_cx);
const rdx_rec_t *rdx_find(const rdx_t *x, int64_t rec_voff);
void   rdx_free(rdx_t *x);

/* Read the 1-based rows[0..n) of the record r describes into out (format 3,
 * inflated, unit 8, in the order given). Every row must be in range; the
 * caller checks against r->nrows. Returns 1 on success, 0 when the file no
 * longer matches r -- the caller then reads the record in full. */
int    rdx_read_rows(BGZF *fh, const rdx_rec_t *r, const int64_t *rows,
                     int64_t n, cdata_t *out);

#endif /* _YAME_RDX_H */
