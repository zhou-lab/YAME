/* summary_multi.c -- see summary_multi.h.
 *
 * SPDX-License-Identifier: AGPL-3.0-or-later
 */
#include <stdlib.h>
#include <string.h>
#include "summary_multi.h"
#include "wzmisc.h"

#define WBITS 64
static inline uint64_t nwords(uint64_t n) { return (n + WBITS - 1) / WBITS; }

/* popcount and count-trailing-zeros without assuming a builtin exists. gcc and
 * clang both fold these to one instruction; the fallbacks are for anything
 * else that can compile the rest of this tree. */
#if defined(__GNUC__) || defined(__clang__)
static inline int pc64(uint64_t x)  { return __builtin_popcountll(x); }
static inline int ctz64(uint64_t x) { return __builtin_ctzll(x); }
#else
static inline int pc64(uint64_t x) {
  int n = 0; while (x) { x &= x - 1; ++n; } return n;
}
static inline int ctz64(uint64_t x) {
  int n = 0; while (!(x & 1)) { x >>= 1; ++n; } return n;
}
#endif

/* Load 64 rows of a format 0 bitmap starting at row `base`.
 *
 * Byte-by-byte rather than a cast: FMT0_IN_SET puts row i at bit (i&7) of byte
 * (i>>3), so a uint64_t load only lines up with that on a little-endian
 * machine. Written this way the compiler emits a plain load where it can and
 * the result is the same everywhere.
 */
static inline uint64_t load0(const uint8_t *s, uint64_t nb, uint64_t byte) {
  uint64_t w = 0;
  uint64_t have = nb > byte ? nb - byte : 0;
  if (have > 8) have = 8;
  for (uint64_t k = 0; k < have; ++k) w |= (uint64_t) s[byte + k] << (k * 8);
  return w;
}

/* The same 64 rows of a format 6 record, reduced to "is a member".
 *
 * Two bits a row, set at the even bit and universe at the odd one, so a member
 * is `set & universe`. One byte holds four rows: (b & (b >> 1)) & 0x55 leaves a
 * member bit at each even position, which then packs down to four bits. Sixteen
 * bytes make 64 rows.
 */
static inline uint64_t load6(const uint8_t *s, uint64_t nb, uint64_t base_row) {
  uint64_t w = 0;
  for (int k = 0; k < 16; ++k) {
    uint64_t byte = (base_row >> 2) + (uint64_t) k;
    if (byte >= nb) break;
    uint8_t b = s[byte];
    uint8_t m = (uint8_t) ((b & (b >> 1)) & 0x55);   /* member at bits 0,2,4,6 */
    uint64_t nib = (uint64_t) ((m & 1) | ((m >> 1) & 2) | ((m >> 2) & 4) | ((m >> 3) & 8));
    w |= nib << (k * 4);
  }
  return w;
}

/* Bits above the record's last row, which must not be counted. */
static inline uint64_t tail_mask(uint64_t n, uint64_t word) {
  uint64_t base = word * WBITS;
  if (base + WBITS <= n) return ~(uint64_t) 0;
  if (base >= n) return 0;
  return (~(uint64_t) 0) >> (WBITS - (n - base));
}

int yame_qbits_build(const cdata_t *query, yame_qbits_t *out) {
  if (!query || !out) return -1;
  if (query->fmt != '3' && query->fmt != '6') return -1;

  memset(out, 0, sizeof *out);
  out->q = query;
  out->n_words = nwords(query->n);
  uint64_t alloc = out->n_words ? out->n_words : 1;
  out->cov = wzcalloc(alloc, sizeof(uint64_t));

  /* The one pass over the query. Row order is irrelevant here -- this marks and
   * counts, it does not sum -- so nothing about bit-identity depends on it. */
  if (query->fmt == '3') {
    for (uint64_t i = 0; i < query->n; ++i)
      if (f3_get_mu((cdata_t *) query, i)) {
        out->cov[i / WBITS] |= (uint64_t) 1 << (i % WBITS);
        out->n_q++;
      }
  } else {
    /* Both bitmaps: the universe decides which rows count at all, and the set
     * is what is being measured inside it. Neither follows from the other. */
    out->qset = wzcalloc(alloc, sizeof(uint64_t));
    for (uint64_t i = 0; i < query->n; ++i) {
      if (FMT6_IN_UNI(*query, i)) out->cov[i / WBITS]  |= (uint64_t) 1 << (i % WBITS);
      if (FMT6_IN_SET(*query, i)) out->qset[i / WBITS] |= (uint64_t) 1 << (i % WBITS);
    }
    /* n_q is per-mask here, since a format 6 mask narrows the universe. */
  }
  return 0;
}

void yame_qbits_free(yame_qbits_t *qb) {
  if (!qb) return;
  free(qb->cov);  qb->cov  = NULL;
  free(qb->qset); qb->qset = NULL;
  qb->n_words = 0;
}

/* The mask's universe, for a format 6 mask: the odd bits of each pair. A binary
 * mask has none, and every row it covers is in scope. */
static inline uint64_t load6_uni(const uint8_t *s, uint64_t nb, uint64_t base_row) {
  uint64_t w = 0;
  for (int k = 0; k < 16; ++k) {
    uint64_t byte = (base_row >> 2) + (uint64_t) k;
    if (byte >= nb) break;
    uint8_t b = (uint8_t) (s[byte] >> 1);            /* universe at bits 1,3,5,7 */
    uint8_t m = (uint8_t) (b & 0x55);
    uint64_t nib = (uint64_t) ((m & 1) | ((m >> 1) & 2) | ((m >> 2) & 4) | ((m >> 3) & 8));
    w |= nib << (k * 4);
  }
  return w;
}

/* The mask's SET bits alone, without ANDing in its universe. */
static inline uint64_t load6_set(const uint8_t *s, uint64_t nb, uint64_t base_row) {
  uint64_t w = 0;
  for (int k = 0; k < 16; ++k) {
    uint64_t byte = (base_row >> 2) + (uint64_t) k;
    if (byte >= nb) break;
    uint8_t m = (uint8_t) (s[byte] & 0x55);          /* set at bits 0,2,4,6 */
    uint64_t nib = (uint64_t) ((m & 1) | ((m >> 1) & 2) | ((m >> 2) & 4) | ((m >> 3) & 8));
    w |= nib << (k * 4);
  }
  return w;
}


/* The mask, however the caller holds it.
 *
 * A record is read 64 rows at a time by the loads above. A caller that came in
 * through the enumerator has already built a bitmap, and that bitmap is WORDS:
 * reading it as a format 0 record would only line up on a little-endian
 * machine, since FMT0_IN_SET addresses bytes. So both shapes go through this
 * instead of one pretending to be the other.
 */
typedef struct {
  const uint64_t *bm;       /* non-NULL: already a bitmap, one word per 64 rows */
  cdata_t        *cx;       /* else: a record */
  uint64_t        nb;       /* bytes of that record's payload */
  int             is6;      /* the record is format 6 */
} mview_t;

static inline uint64_t mv_member(const mview_t *v, uint64_t w, uint64_t base) {
  if (v->bm) return v->bm[w];
  return v->is6 ? load6(v->cx->s, v->nb, base) : load0(v->cx->s, v->nb, base / 8);
}
/* A bitmap and a binary record carry no universe of their own: every row is in
 * scope, and the query's universe is the only one that narrows anything. */
static inline uint64_t mv_uni(const mview_t *v, uint64_t w, uint64_t base) {
  (void) w;
  if (v->bm || !v->is6) return ~(uint64_t) 0;
  return load6_uni(v->cx->s, v->nb, base);
}
static inline uint64_t mv_set(const mview_t *v, uint64_t w, uint64_t base) {
  if (v->bm) return v->bm[w];
  return v->is6 ? load6_set(v->cx->s, v->nb, base) : load0(v->cx->s, v->nb, base / 8);
}

/* A format 3 query: the universe is every row, the mask count is unrestricted,
 * and beta is the MEAN of the per-row betas -- so the values have to be read,
 * at the overlapping rows only. */
static int one_q3(const yame_qbits_t *qb, const mview_t *mv, yame_acc_t *acc) {
  uint64_t n = qb->q->n;

  acc->n_u = n;
  acc->n_q = qb->n_q;

  for (uint64_t w = 0; w < qb->n_words; ++w) {
    uint64_t base = w * WBITS;
    uint64_t mw = mv_member(mv, w, base);
    mw &= tail_mask(n, w);
    if (!mw) continue;

    acc->n_m += (uint64_t) pc64(mw);

    uint64_t ov = mw & qb->cov[w];
    if (!ov) continue;
    acc->n_o += (uint64_t) pc64(ov);

    /* Increasing row order, which is what keeps the sum bit-identical to the
     * per-mask path: ctz takes the lowest set bit first, and clearing it with
     * ov &= ov - 1 moves up. */
    while (ov) {
      int b = ctz64(ov);
      uint64_t mu = f3_get_mu((cdata_t *) qb->q, base + (uint64_t) b);
      acc->sum_depth += MU2cov(mu);
      acc->sum_beta  += MU2beta(mu);
      ov &= ov - 1;
    }
  }
  acc->beta = acc->sum_beta / acc->n_o;     /* Inf at zero overlap, as before */
  return 0;
}

/* A format 6 query: every count is taken INSIDE a universe, and beta is the
 * overlap over the mask count. No value is read at all -- four popcounts a
 * word, which is why this is cheaper than the format 3 path rather than dearer.
 */
static int one_q6(const yame_qbits_t *qb, const mview_t *mv, yame_acc_t *acc) {
  uint64_t n = qb->q->n;

  for (uint64_t w = 0; w < qb->n_words; ++w) {
    uint64_t base = w * WBITS, keep = tail_mask(n, w);
    /* The universe: the query's, narrowed by the mask's when it has one. */
    uint64_t uni = qb->cov[w] & keep & mv_uni(mv, w, base);
    uint64_t mset = mv_set(mv, w, base);
    if (!uni) continue;

    uint64_t qs = qb->qset[w] & uni;
    uint64_t ms = mset & uni;
    acc->n_u += (uint64_t) pc64(uni);
    acc->n_q += (uint64_t) pc64(qs);
    acc->n_m += (uint64_t) pc64(ms);
    acc->n_o += (uint64_t) pc64(qs & ms);
  }
  acc->beta = (double) acc->n_o / acc->n_m; /* Inf at an empty mask, as before */
  return 0;
}

static int one_view(const yame_qbits_t *qb, const mview_t *mv, yame_acc_t *acc) {
  memset(acc, 0, sizeof *acc);
  if (qb->q->fmt == '3') return one_q3(qb, mv, acc);
  if (qb->q->fmt == '6' && qb->qset) return one_q6(qb, mv, acc);
  return -1;
}

int yame_summarize_one_cx(const yame_qbits_t *qb, cdata_t *mask, yame_acc_t *acc) {
  if (!qb || !qb->cov || !mask || !acc) return -1;
  if (mask->fmt != '0' && mask->fmt != '1' && mask->fmt != '6') return -1;
  if (mask->n != qb->q->n) return -1;

  mview_t mv;
  mv.bm  = NULL;
  mv.cx  = mask;
  mv.is6 = (mask->fmt == '6');
  mv.nb  = mv.is6 ? (qb->q->n + 3) / 4 : (qb->q->n + 7) / 8;
  return one_view(qb, &mv, acc);
}

/* ------------------------------------------------------------- record form -- */

int yame_summarize_multi_cx(const cdata_t *query, cdata_t *masks, uint32_t n_masks,
                            yame_acc_t *acc, uint64_t *n_q_out) {
  yame_qbits_t qb;
  if (yame_qbits_build(query, &qb) != 0) return -1;
  for (uint32_t j = 0; j < n_masks; ++j) {
    if (yame_summarize_one_cx(&qb, &masks[j], &acc[j]) != 0) {
      yame_qbits_free(&qb);
      return -1;
    }
  }
  if (n_q_out) *n_q_out = qb.n_q;
  yame_qbits_free(&qb);
  return 0;
}

/* --------------------------------------------------------- enumerator form -- */

/* Runs land in a per-slot bitmap, then the same word-wise accumulation runs
 * over it. A long run sets whole words at once, which is why a caller with runs
 * pays less here than one with records: there is nothing to scan. */
typedef struct {
  uint64_t **bm;            /* one bitmap per slot, allocated on first use */
  uint64_t   n_words, n;
  uint32_t   n_masks, n_acc;
  const uint32_t *state_base;
  int        bad;
} build_t;

static void build_emit(void *ctx, uint32_t mask,
                       uint64_t start, uint64_t len, uint32_t state) {
  build_t *b = ctx;
  if (!len || mask >= b->n_masks) { if (mask >= b->n_masks) b->bad = 1; return; }
  uint32_t slot = b->state_base ? b->state_base[mask] + state : mask;
  if (slot >= b->n_acc) { b->bad = 1; return; }

  if (start >= b->n) return;                     /* past the query: clipped */
  uint64_t end = start + len;
  if (end > b->n) end = b->n;

  if (!b->bm[slot]) b->bm[slot] = wzcalloc(b->n_words, sizeof(uint64_t));
  uint64_t *m = b->bm[slot];

  /* Whole words in the middle, partial words at the two ends. */
  uint64_t i = start;
  while (i < end) {
    uint64_t w = i / WBITS, off = i % WBITS;
    uint64_t room = WBITS - off, take = end - i < room ? end - i : room;
    uint64_t bits = (take == WBITS) ? ~(uint64_t) 0
                                    : (((uint64_t) 1 << take) - 1) << off;
    m[w] |= bits;
    i += take;
  }
}

int yame_summarize_multi(const cdata_t *query,
                         yame_runs_fn runs, void *runs_ctx,
                         uint32_t n_masks, const uint32_t *state_base,
                         yame_acc_t *acc, uint32_t n_acc,
                         uint64_t *n_q_out) {
  if (!query || !runs || !acc || !n_acc) return -1;

  yame_qbits_t qb;
  if (yame_qbits_build(query, &qb) != 0) return -1;

  build_t b;
  memset(&b, 0, sizeof b);
  b.n_words = qb.n_words; b.n = query->n;
  b.n_masks = n_masks; b.n_acc = n_acc; b.state_base = state_base;
  b.bm = wzcalloc(n_acc, sizeof(uint64_t *));

  int rc = runs(runs_ctx, build_emit, &b);
  if (rc != 0 || b.bad) {
    for (uint32_t s = 0; s < n_acc; ++s) free(b.bm[s]);
    free(b.bm); yame_qbits_free(&qb);
    return -1;
  }

  /* The same per-mask code as `summary -m` runs here, over the bitmap the runs
   * built: one definition of the arithmetic, not two. */
  for (uint32_t s = 0; s < n_acc; ++s) {
    if (!b.bm[s]) continue;
    mview_t mv; mv.bm = b.bm[s]; mv.cx = NULL; mv.nb = 0; mv.is6 = 0;
    if (one_view(&qb, &mv, &acc[s]) != 0) {
      for (uint32_t t = s; t < n_acc; ++t) free(b.bm[t]);
      free(b.bm); yame_qbits_free(&qb);
      return -1;
    }
    free(b.bm[s]);
  }
  free(b.bm);
  if (n_q_out) *n_q_out = qb.n_q;
  yame_qbits_free(&qb);
  return 0;
}
