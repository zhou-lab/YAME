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
  /* A COMPRESSED record counts bytes in `n`, not rows, while `unit` still
   * describes the inflated width -- so reading row i of a compressed format 3
   * record at unit 8 reaches eight times past the buffer. `summary` inflates
   * before it gets here, but this is a public entry point and the caller cannot
   * be assumed to have. */
  if (query->compressed) return -1;

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
  acc->beta = acc->sum_beta / acc->n_o;     /* NaN at zero overlap; the
                                             * per-mask path gives NaN too, and
                                             * both print NA */
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
  acc->beta = (double) acc->n_o / acc->n_m; /* NaN at an empty mask, as the
                                             * per-mask path also gives */
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
  /* NOT format 1. An inflated format 1 record is one BYTE per row, not one bit,
   * so reading it through load0() stays in bounds and answers zero -- silently.
   * `summary` never offers one, because prepare_mask() converts it to format 0
   * first, but a caller that has not done that gets a refusal rather than a
   * plausible wrong number. */
  if (mask->fmt != '0' && mask->fmt != '6') return -1;
  /* NO compressed test here, deliberately. For a QUERY the flag matters: a
   * compressed format 3 record counts bytes in `n`. For a format 0 or 6 MASK it
   * does not describe the layout at all -- convertToFmt0() hands back a format 0
   * record with compressed still set, because bit-packed IS its storage. A test
   * here rejected almost every mask and sent the whole knowledgebase down the
   * fallback: 6 s became 107 s, with the output still correct, so only a timing
   * check caught it. */
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

/*
 * Runs are consumed WHERE THEY ARRIVE. Nothing per-slot is allocated.
 *
 * The first version of this turned each slot's runs into a bitmap and reused
 * the word-wise code above. That is right for a handful of binary masks and
 * catastrophic for the case it was actually built for: a bitmap is n/8 bytes,
 * 3.7 MB on a whole-genome row space, and an MRMP bank is mostly STATE masks --
 * 178 sets with even five states apiece is 890 slots, so 3.3 GB of bitmaps
 * before a single count. Twenty states would be 13 GB.
 *
 * So each run is scored as it is emitted. The claimed count is its length. The
 * overlap is a popcount of the query's covered bits inside its range. The sums
 * read values only at those bits. Nothing scales with the slot count except the
 * accumulators, which are 56 bytes each.
 *
 * A consequence worth knowing: a row claimed by several slots has its value read
 * once per slot. For a bank where sets barely overlap that is the minimum work
 * anyway, and caching it would cost more than it saves.
 */
typedef struct {
  const yame_qbits_t *qb;
  yame_acc_t *acc;
  uint32_t    n_masks, n_acc;
  const uint32_t *state_base;
  int         is_q3;
  int         bad;
} scat_t;

/* Covered bits of [lo, hi) as a word and its in-range mask, for one word. */
static inline uint64_t range_word(const yame_qbits_t *qb, uint64_t w,
                                  uint64_t lo, uint64_t hi) {
  uint64_t base = w * WBITS;
  uint64_t keep = ~(uint64_t) 0;
  if (lo > base) keep &= ~(uint64_t) 0 << (lo - base);
  if (hi < base + WBITS) {
    uint64_t take = hi - base;
    keep &= (take >= WBITS) ? ~(uint64_t) 0 : (((uint64_t) 1 << take) - 1);
  }
  return qb->cov[w] & keep;
}

static void scat_emit(void *ctx, uint32_t mask,
                      uint64_t start, uint64_t len, uint32_t state) {
  scat_t *k = ctx;
  if (!len) return;
  if (mask >= k->n_masks) { k->bad = 1; return; }
  uint32_t slot = k->state_base ? k->state_base[mask] + state : mask;
  if (slot >= k->n_acc) { k->bad = 1; return; }

  uint64_t n = k->qb->q->n;
  if (start >= n) return;                        /* past the query: clipped */
  /* start + len can WRAP. It did: len == UINT64_MAX made end 0, the clamp below
   * never fired, n_m underflowed, and the sweep then indexed cov[] about 2^58
   * words past its allocation -- a guaranteed crash, not an overread. Clamping
   * before the addition is the only place that can catch it, because after the
   * wrap the value looks legitimate. */
  uint64_t end = (len > n - start) ? n : start + len;

  yame_acc_t *a = &k->acc[slot];
  a->n_m += end - start;

  uint64_t w0 = start / WBITS, w1 = (end - 1) / WBITS;
  for (uint64_t w = w0; w <= w1; ++w) {
    uint64_t ov = range_word(k->qb, w, start, end);
    if (!ov) continue;
    a->n_o += (uint64_t) pc64(ov);
    if (!k->is_q3) continue;                     /* fmt6: counts are the answer */
    uint64_t base = w * WBITS;
    /* Increasing row order within the run, and runs of one slot arrive in
     * increasing order from any sane enumerator, so the sum matches the
     * per-mask path. */
    while (ov) {
      int b = ctz64(ov);
      uint64_t mu = f3_get_mu((cdata_t *) k->qb->q, base + (uint64_t) b);
      a->sum_depth += MU2cov(mu);
      a->sum_beta  += MU2beta(mu);
      ov &= ov - 1;
    }
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

  /* A format 6 query needs its SET bits too, and a run says nothing about them,
   * so the enumerator form is format 3 only. A caller with a format 6 query and
   * records uses yame_summarize_one_cx(). */
  if (query->fmt != '3') { yame_qbits_free(&qb); return -1; }

  memset(acc, 0, (size_t) n_acc * sizeof(yame_acc_t));

  scat_t k;
  memset(&k, 0, sizeof k);
  k.qb = &qb; k.acc = acc;
  k.n_masks = n_masks; k.n_acc = n_acc; k.state_base = state_base;
  k.is_q3 = 1;

  int rc = runs(runs_ctx, scat_emit, &k);
  if (rc != 0 || k.bad) { yame_qbits_free(&qb); return -1; }

  for (uint32_t s = 0; s < n_acc; ++s) {
    acc[s].n_u = query->n;
    acc[s].n_q = qb.n_q;
    acc[s].beta = acc[s].sum_beta / acc[s].n_o;   /* NaN at zero overlap */
  }

  if (n_q_out) *n_q_out = qb.n_q;
  yame_qbits_free(&qb);
  return 0;
}
