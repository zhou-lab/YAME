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
 * The row index. See rdx.h for what it is and why.
 *
 * It lives in the .idx, as a third tab-separated column after the record's
 * name and offset (rdx_encode() below spells it out). Every reader of the
 * two-column form takes the first two fields and stops, so old code reads a
 * new .idx unchanged; a line without the column, or with ".", is a record
 * the reader reads in full.
 *
 * Building needs the virtual offset of a byte INSIDE a record's stream, which
 * the reader that inflated it does not report. It is arithmetic instead: a
 * BGZF block's header gives its compressed size and its footer its inflated
 * size, so one pass over the headers -- no inflate -- places every
 * uncompressed byte of the file. The record's own offset places its stream.
 */

#include <errno.h>
#include <zlib.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include "rdx.h"
#include "cfile.h"
#include "index.h"
#include "kstring.h"
#include "wzio.h"
#include "wzmisc.h"

#define RDX_MAGIC "YAMERDX1"
#define RDX_HDR 17   /* sig(8) + fmt(1) + n(8) in front of each stream */

/* One format-3 token at p, with at most `avail` bytes readable. Returns its
 * length in bytes (0 when it runs past avail) and the rows it covers, with
 * the M/U it holds (0/0 for a zero run). */
static int f3_token(const uint8_t *p, uint64_t avail, uint64_t *rows,
                    uint64_t *M, uint64_t *U) {
  if (avail < 1) return 0;
  switch (p[0] & 0x3) {
  case 0: {
    if (avail < 2) return 0;
    uint16_t v; memcpy(&v, p, 2);
    *rows = v >> 2; *M = 0; *U = 0;
    return 2; }
  case 1:
    *rows = 1; *M = p[0] >> 5; *U = (p[0] >> 2) & 0x7;
    return 1;
  case 2: {
    if (avail < 2) return 0;
    uint16_t v; memcpy(&v, p, 2);
    *rows = 1; *M = v >> 9; *U = (v >> 2) & 0x7f;
    return 2; }
  default: {
    if (avail < 8) return 0;
    uint64_t v; memcpy(&v, p, 8);
    *rows = 1; *M = v >> 33; *U = (v >> 2) & ((1ul << 31) - 1);
    return 8; }
  }
}

/* ---- builder ------------------------------------------------------------ */

static void blk_push(rdx_t *x, uint64_t *m, int64_t coff, uint64_t upos,
                     uint32_t usize) {
  if (x->nblk == *m) {
    *m = *m ? *m * 2 : 1024;
    x->blk_coff  = wzrealloc(x->blk_coff,  *m * sizeof(int64_t));
    x->blk_upos  = wzrealloc(x->blk_upos,  *m * sizeof(uint64_t));
    x->blk_usize = wzrealloc(x->blk_usize, *m * sizeof(uint32_t));
  }
  x->blk_coff[x->nblk] = coff;
  x->blk_upos[x->nblk] = upos;
  x->blk_usize[x->nblk] = usize;
  x->nblk++;
}

rdx_t *rdx_open_build(const char *fname) {
  if (!fname || strcmp(fname, "-") == 0) return NULL;
  int fd = open(fname, O_RDONLY);
  if (fd < 0) return NULL;
  struct stat st;
  if (fstat(fd, &st) != 0 || !S_ISREG(st.st_mode)) { close(fd); return NULL; }

  rdx_t *x = wzcalloc(1, sizeof(rdx_t));
  uint64_t m = 0, upos = 0;
  int64_t o = 0;
  while (o < st.st_size) {
    uint8_t h[18];
    uint8_t f[4];
    if (pread(fd, h, 18, o) != 18) goto fail;
    /* BGZF as htslib writes it: gzip with one 'BC' extra field, XLEN 6 */
    if (h[0] != 31 || h[1] != 139 || h[2] != 8 || !(h[3] & 4) ||
        h[10] != 6 || h[11] != 0 || h[12] != 'B' || h[13] != 'C' ||
        h[14] != 2 || h[15] != 0) goto fail;
    int64_t bsize = (int64_t) (h[16] | (h[17] << 8)) + 1;
    if (bsize < 28 || o + bsize > st.st_size) goto fail;
    if (pread(fd, f, 4, o + bsize - 4) != 4) goto fail;
    uint32_t isize = (uint32_t) f[0] | ((uint32_t) f[1] << 8) |
      ((uint32_t) f[2] << 16) | ((uint32_t) f[3] << 24);
    if (isize > 65536) goto fail;
    blk_push(x, &m, o, upos, isize);
    upos += isize;
    o += bsize;
  }
  close(fd);
  return x;

fail:
  close(fd);
  rdx_free(x);
  return NULL;
}

/* the block whose bytes include uncompressed position u, at or after b */
static uint64_t blk_of(const rdx_t *x, uint64_t b, uint64_t u) {
  while (b < x->nblk && u >= x->blk_upos[b] + x->blk_usize[b]) ++b;
  return b;
}

static rdx_rec_t *rec_push(rdx_t *x) {
  if (x->nrec == x->mrec) {
    x->mrec = x->mrec ? x->mrec * 2 : 16;
    x->rec = wzrealloc(x->rec, x->mrec * sizeof(rdx_rec_t));
  }
  rdx_rec_t *r = &x->rec[x->nrec++];
  memset(r, 0, sizeof(*r));
  return r;
}

void rdx_add(rdx_t *x, int64_t rec_voff, const cdata_t *c) {
  if (!x || c->fmt != '3' || !c->compressed || c->n == 0) return;

  /* the block the record's offset names, by exact compressed address */
  int64_t coff = rec_voff >> 16;
  uint64_t lo = 0, hi = x->nblk;
  while (lo < hi) {
    uint64_t mid = (lo + hi) / 2;
    if (x->blk_coff[mid] < coff) lo = mid + 1; else hi = mid;
  }
  if (lo == x->nblk || x->blk_coff[lo] != coff) return;
  uint64_t s0 = x->blk_upos[lo] + (uint64_t) (rec_voff & 0xFFFF) + RDX_HDR;

  rdx_ckpt_t *ck = NULL;
  uint64_t nck = 0, mck = 0, last_b = UINT64_MAX;
  uint64_t b = lo, spos = 0, row = 0;
  while (spos < c->n) {
    uint64_t u = s0 + spos;
    b = blk_of(x, b, u);
    if (b == x->nblk) { free(ck); return; }   /* the file ends early */
    if (b != last_b) {
      if (nck == mck) {
        mck = mck ? mck * 2 : 64;
        ck = wzrealloc(ck, mck * sizeof(rdx_ckpt_t));
      }
      ck[nck].row0 = row;
      ck[nck].spos = spos;
      ck[nck].voff = (x->blk_coff[b] << 16) | (int64_t) (u - x->blk_upos[b]);
      ck[nck].byte0 = c->s[spos];
      ++nck;
      last_b = b;
    }
    uint64_t rows, M, U;
    int len = f3_token(c->s + spos, c->n - spos, &rows, &M, &U);
    if (!len) { free(ck); return; }           /* a malformed stream */
    spos += len;
    row += rows;
  }

  /* where a reader stands after this record: inside the last block, or --
   * when the stream fills it exactly -- at the block after it, empty member
   * or not, which is what bgzf_tell() then reports */
  uint64_t ue = s0 + c->n;
  uint64_t bl = blk_of(x, b, ue - 1);
  int64_t end_voff = -1;
  if (bl < x->nblk) {
    if (ue < x->blk_upos[bl] + x->blk_usize[bl])
      end_voff = (x->blk_coff[bl] << 16) | (int64_t) (ue - x->blk_upos[bl]);
    else if (bl + 1 < x->nblk)
      end_voff = x->blk_coff[bl + 1] << 16;
  }

  rdx_rec_t *r = rec_push(x);
  r->rec_voff = rec_voff;
  r->end_voff = end_voff;
  r->nbytes = c->n;
  r->nrows = row;
  r->nckpt = nck;
  r->ckpt = ck;
}

static int64_t voff_shift(int64_t v, int64_t coff_delta) {
  return (((v >> 16) + coff_delta) << 16) | (v & 0xFFFF);
}

void rdx_add_moved(rdx_t *x, const rdx_rec_t *r, int64_t new_voff,
                   int64_t coff_delta, int64_t new_end_voff) {
  rdx_rec_t *o = rec_push(x);
  *o = *r;
  o->rec_voff = new_voff;
  o->end_voff = new_end_voff;
  o->ckpt = wzmalloc((r->nckpt ? r->nckpt : 1) * sizeof(rdx_ckpt_t));
  for (uint64_t i = 0; i < r->nckpt; ++i) {
    o->ckpt[i] = r->ckpt[i];
    o->ckpt[i].voff = voff_shift(r->ckpt[i].voff, coff_delta);
  }
}

/* ---- the .idx third column -------------------------------------------- */

/* Base 36, lower case: digits and letters only, so no reader that guesses a
 * separator (data.table::fread tries , ; : | and space) can find one here. */
static void put36(kstring_t *ks, uint64_t v) {
  char b[16]; int n = 0;
  do { int d = (int) (v % 36); b[n++] = (char) (d < 10 ? '0' + d : 'a' + d - 10);
       v /= 36; } while (v);
  while (n) kputc(b[--n], ks);
}

static int get36(const char **p, uint64_t *v) {
  const char *q = *p; uint64_t x = 0; int n = 0;
  for (;; ++q, ++n) {
    int d;
    if (*q >= '0' && *q <= '9') d = *q - '0';
    else if (*q >= 'a' && *q <= 'z') d = *q - 'a' + 10;
    else break;
    if (x > (UINT64_MAX - (uint64_t) d) / 36) return 0;
    x = x * 36 + (uint64_t) d;
  }
  if (!n) return 0;
  *p = q; *v = x;
  return 1;
}

static int want(const char **p, char c) {
  if (**p != c) return 0;
  ++*p; return 1;
}

/* r1.<nbytes>.<nrows>.<end>.<ckpt>.<ckpt>....k<crc>
 *   end   = x when unknown, else <coff - rec coff>_<uoff>
 *   ckpt  = <drow>_<dspos>_<dcoff>_<uoff>_<byte0>, each d a difference from
 *           the checkpoint before (the first from row 0, byte 0 and the
 *           record's own compressed address)
 *   crc   = CRC32 of everything before ".k". A digit edited by hand or by a
 *           tool would otherwise shift row numbers past every check the
 *           reader can make cheaply, and give wrong rows with exit 0. */
static void rdx_encode(kstring_t *ks, const rdx_rec_t *r) {
  int64_t c0 = r->rec_voff >> 16;
  size_t l0 = ks->l;
  kputs("r1.", ks);
  put36(ks, r->nbytes); kputc('.', ks);
  put36(ks, r->nrows); kputc('.', ks);
  if (r->end_voff < 0) kputc('x', ks);
  else {
    put36(ks, (uint64_t) ((r->end_voff >> 16) - c0)); kputc('_', ks);
    put36(ks, (uint64_t) (r->end_voff & 0xFFFF));
  }
  uint64_t prow = 0, pspos = 0; int64_t pcoff = c0;
  for (uint64_t i = 0; i < r->nckpt; ++i) {
    const rdx_ckpt_t *k = &r->ckpt[i];
    kputc('.', ks);
    put36(ks, k->row0 - prow); kputc('_', ks);
    put36(ks, k->spos - pspos); kputc('_', ks);
    put36(ks, (uint64_t) ((k->voff >> 16) - pcoff)); kputc('_', ks);
    put36(ks, (uint64_t) (k->voff & 0xFFFF)); kputc('_', ks);
    put36(ks, k->byte0);
    prow = k->row0; pspos = k->spos; pcoff = k->voff >> 16;
  }
  uint32_t crc = (uint32_t) crc32(0L, (const Bytef *) ks->s + l0,
                                  (uInt) (ks->l - l0));
  kputs(".k", ks);
  put36(ks, crc);
}

/* 1 when s is a column this build writes and r now holds it */
static int rdx_decode(const char *s, int64_t rec_voff, rdx_rec_t *r) {
  memset(r, 0, sizeof(*r));
  r->rec_voff = rec_voff;
  int64_t c0 = rec_voff >> 16;
  const char *p = s;
  uint64_t a, b, c, d, e;
  if (strncmp(p, "r1.", 3) != 0) return 0;
  /* the checksum first: nothing below is believed unless it matches */
  /* the LAST '.': a base-36 field may itself begin with k */
  const char *ck = strrchr(s, '.');
  if (!ck || ck[1] != 'k') return 0;
  { const char *q = ck + 2; uint64_t want_crc;
    if (!get36(&q, &want_crc) || *q != '\0') return 0;
    uint32_t got = (uint32_t) crc32(0L, (const Bytef *) s, (uInt) (ck - s));
    if ((uint64_t) got != want_crc) return 0; }
  p += 3;
  if (!get36(&p, &r->nbytes) || !want(&p, '.') ||
      !get36(&p, &r->nrows) || !want(&p, '.')) return 0;
  if (*p == 'x') { r->end_voff = -1; ++p; }
  else {
    if (!get36(&p, &a) || !want(&p, '_') || !get36(&p, &b) || b > 0xFFFF)
      return 0;
    r->end_voff = ((c0 + (int64_t) a) << 16) | (int64_t) b;
  }
  uint64_t m = 0, prow = 0, pspos = 0; int64_t pcoff = c0;
  while (p < ck && *p == '.') {
    ++p;
    if (!get36(&p, &a) || !want(&p, '_') || !get36(&p, &b) || !want(&p, '_') ||
        !get36(&p, &c) || !want(&p, '_') || !get36(&p, &d) || !want(&p, '_') ||
        !get36(&p, &e) || d > 0xFFFF || e > 0xFF) goto fail;
    if (r->nckpt == m) {
      m = m ? m * 2 : 64;
      r->ckpt = wzrealloc(r->ckpt, m * sizeof(rdx_ckpt_t));
    }
    rdx_ckpt_t *k = &r->ckpt[r->nckpt++];
    k->row0 = prow + a; k->spos = pspos + b;
    k->voff = ((pcoff + (int64_t) c) << 16) | (int64_t) d;
    k->byte0 = (uint8_t) e;
    prow = k->row0; pspos = k->spos; pcoff = k->voff >> 16;
  }
  if (p != ck) goto fail;
  /* checkpoints climb, and the first is the stream's first token, or the
   * searches the reader does are meaningless */
  if (r->nckpt == 0 || r->ckpt[0].row0 != 0 || r->ckpt[0].spos != 0)
    goto fail;
  for (uint64_t i = 1; i < r->nckpt; ++i)
    if (r->ckpt[i].spos <= r->ckpt[i-1].spos ||
        r->ckpt[i].row0 < r->ckpt[i-1].row0) goto fail;
  if (r->ckpt[r->nckpt - 1].spos >= r->nbytes ||
      r->ckpt[r->nckpt - 1].row0 >= r->nrows) goto fail;
  return 1;
fail:
  free(r->ckpt); r->ckpt = NULL; r->nckpt = 0;
  return 0;
}

static int cmp_rec(const void *a, const void *b) {
  int64_t x = ((const rdx_rec_t *) a)->rec_voff;
  int64_t y = ((const rdx_rec_t *) b)->rec_voff;
  return (x > y) - (x < y);
}

static void rdx_sort(rdx_t *x) {
  if (x) qsort(x->rec, x->nrec, sizeof(rdx_rec_t), cmp_rec);
}

void rdx_write_index(FILE *fp, index_t *idx, rdx_t *x) {
  if (!x || x->nrec == 0) { writeIndex(fp, idx); return; }
  rdx_sort(x);
  int n = 0;
  index_pair_t *pairs = index_pairs(idx, &n);
  kstring_t ks = {0};
  for (int i = 0; i < n; ++i) {
    const rdx_rec_t *r = rdx_find(x, pairs[i].value);
    ks.l = 0;
    if (r) rdx_encode(&ks, r);
    else kputc('.', &ks);       /* every line has the column, or read.table
                                 * stops at the first line without one */
    fprintf(fp, "%s\t%"PRId64"\t%s\n", pairs[i].key, pairs[i].value, ks.s);
  }
  free(ks.s);
  clean_index_pairs(pairs, n);
}

rdx_t *rdx_load_index(const char *fname_cx) {
  if (!fname_cx || strcmp(fname_cx, "-") == 0) return NULL;
  char *fn = get_fname_index(fname_cx);
  gzFile fh = fn ? gzopen(fn, "r") : NULL;
  free(fn);
  if (!fh) return NULL;
  rdx_t *x = NULL;
  char *line = NULL;
  int more = 1;
  while (more) {
    more = gzFile_read_line(fh, &line);
    char *t1 = strchr(line, '\t');
    if (!t1) continue;
    char *t2 = strchr(t1 + 1, '\t');
    if (!t2 || t2[1] == '.' || t2[1] == '\0') continue;
    char *end = NULL;
    int64_t voff = strtoll(t1 + 1, &end, 10);
    if (end != t2) continue;
    rdx_rec_t r;
    if (!rdx_decode(t2 + 1, voff, &r)) continue;
    if (!x) x = wzcalloc(1, sizeof(rdx_t));
    *rec_push(x) = r;
  }
  free(line);
  gzclose(fh);
  rdx_sort(x);
  return x;
}

const rdx_rec_t *rdx_find(const rdx_t *x, int64_t rec_voff) {
  if (!x) return NULL;
  uint64_t lo = 0, hi = x->nrec;
  while (lo < hi) {
    uint64_t mid = (lo + hi) / 2;
    if (x->rec[mid].rec_voff < rec_voff) lo = mid + 1; else hi = mid;
  }
  if (lo < x->nrec && x->rec[lo].rec_voff == rec_voff) return &x->rec[lo];
  return NULL;
}

void rdx_free(rdx_t *x) {
  if (!x) return;
  for (uint64_t i = 0; i < x->nrec; ++i) free(x->rec[i].ckpt);
  free(x->rec);
  free(x->blk_coff); free(x->blk_upos); free(x->blk_usize);
  free(x);
}

/* ---- reading rows ------------------------------------------------------- */

typedef struct {
  BGZF *fh;
  const rdx_rec_t *r;
  int positioned;
  uint64_t spos;         /* stream position after the current token */
  uint64_t row, rows;    /* current token covers [row, row + rows) */
  uint64_t M, U;
  uint64_t nextck;       /* the first checkpoint not yet passed */
} rdx_cur_t;

/* Read the token at the reader's position. Reading on across a checkpoint
 * is a free check of the index against the data: the token there must have
 * the checkpoint's first byte and start on its row. Any disagreement, or a
 * token running past the record's rows or bytes, means the file is not the
 * one indexed. */
static int cur_next(rdx_cur_t *k) {
  const rdx_rec_t *r = k->r;
  if (k->spos >= r->nbytes) return 0;
  uint8_t t[8];
  if (bgzf_read(k->fh, t, 1) != 1) return 0;
  uint64_t row = k->row + k->rows;  /* the row this token starts on */
  if (k->nextck < r->nckpt) {
    const rdx_ckpt_t *c = &r->ckpt[k->nextck];
    if (k->spos > c->spos) return 0;  /* stepped over it */
    if (k->spos == c->spos) {
      if (t[0] != c->byte0 || row != c->row0) return 0;
      ++k->nextck;
    }
  }
  int len = (t[0] & 0x3) == 1 ? 1 : (t[0] & 0x3) == 3 ? 8 : 2;
  if (k->spos + len > r->nbytes) return 0;
  if (len > 1 && bgzf_read(k->fh, t + 1, len - 1) != len - 1) return 0;
  uint64_t rows = 0;
  if (!f3_token(t, len, &rows, &k->M, &k->U)) return 0;
  if (row + rows > r->nrows) return 0;
  k->row = row;
  k->rows = rows;
  k->spos += len;
  return 1;
}

static int cur_seek(rdx_cur_t *k, uint64_t j) {
  const rdx_ckpt_t *c = &k->r->ckpt[j];
  if (bgzf_seek(k->fh, c->voff, SEEK_SET) < 0) return 0;
  k->positioned = 1;
  k->spos = c->spos;
  k->row = c->row0;
  k->rows = 0;
  k->nextck = j;                /* checked by the read that follows */
  return cur_next(k);
}

/* the last checkpoint with row0 <= w */
static uint64_t ckpt_for_row(const rdx_rec_t *r, uint64_t w) {
  uint64_t lo = 0, hi = r->nckpt;
  while (hi - lo > 1) {
    uint64_t mid = (lo + hi) / 2;
    if (r->ckpt[mid].row0 <= w) lo = mid; else hi = mid;
  }
  return lo;
}

/* the first checkpoint at or after stream position p */
static uint64_t ckpt_at_or_after(const rdx_rec_t *r, uint64_t p) {
  uint64_t lo = 0, hi = r->nckpt;
  while (lo < hi) {
    uint64_t mid = (lo + hi) / 2;
    if (r->ckpt[mid].spos < p) lo = mid + 1; else hi = mid;
  }
  return lo;
}

static const int64_t *sort_rows;
static int cmp_by_row(const void *a, const void *b) {
  int64_t x = sort_rows[*(const int64_t *) a];
  int64_t y = sort_rows[*(const int64_t *) b];
  return (x > y) - (x < y);
}

int rdx_read_rows(BGZF *fh, const rdx_rec_t *r, const int64_t *rows,
                  int64_t n, cdata_t *out) {

  /* the record header must still be the one the index describes */
  if (bgzf_seek(fh, r->rec_voff, SEEK_SET) < 0) return 0;
  if (fh->block_length == 0 && bgzf_read_block(fh) < 0) return 0;
  uint64_t sig, nb; char fmt;
  if (bgzf_read(fh, &sig, 8) != 8 || sig != CDSIG ||
      bgzf_read(fh, &fmt, 1) != 1 || fmt != '3' ||
      bgzf_read(fh, &nb, 8) != 8 || nb != r->nbytes) return 0;

  if (n < 0) return 0;
  size_t nn = n ? (size_t) n : 1;
  int64_t *ord = wzmalloc(nn * sizeof(int64_t));
  for (int64_t i = 0; i < n; ++i) ord[i] = i;
  sort_rows = rows;
  qsort(ord, n, sizeof(int64_t), cmp_by_row);

  cdata_t c = {0};
  c.fmt = '3'; c.unit = 8; c.n = n; c.compressed = 0;
  c.s = wzcalloc(nn, 8);

  rdx_cur_t k = { .fh = fh, .r = r };
  for (int64_t i = 0; i < n; ++i) {
    uint64_t w = (uint64_t) rows[ord[i]] - 1;
    while (!k.positioned || w >= k.row + k.rows) {
      uint64_t j = ckpt_for_row(r, w);
      /* Seek only to skip whole blocks: when the checkpoint is the next one
       * ahead, reading on reaches it without inflating anything twice. */
      int seek = !k.positioned ||
        (r->ckpt[j].spos > k.spos && j != ckpt_at_or_after(r, k.spos));
      if (!(seek ? cur_seek(&k, j) : cur_next(&k))) {
        free(c.s); free(ord); return 0;
      }
    }
    f3_set_mu(&c, ord[i], k.M, k.U);
  }
  free(ord);
  *out = c;
  return 1;
}
