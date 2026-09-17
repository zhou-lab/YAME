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
  if (query->fmt != '3') return -1;

  memset(out, 0, sizeof *out);
  out->q = query;
  out->n_words = nwords(query->n);
  out->cov = wzcalloc(out->n_words ? out->n_words : 1, sizeof(uint64_t));

  /* The one pass over the query. Row order is irrelevant here -- this counts
   * and marks, it does not sum -- so nothing about bit-identity depends on it. */
  for (uint64_t i = 0; i < query->n; ++i)
    if (f3_get_mu((cdata_t *) query, i)) {
      out->cov[i / WBITS] |= (uint64_t) 1 << (i % WBITS);
      out->n_q++;
    }
  return 0;
}

void yame_qbits_free(yame_qbits_t *qb) {
  if (!qb) return;
  free(qb->cov); qb->cov = NULL; qb->n_words = 0;
}

int yame_summarize_one_cx(const yame_qbits_t *qb, cdata_t *mask, yame_acc_t *acc) {
  if (!qb || !qb->cov || !mask || !acc) return -1;
  if (mask->fmt != '0' && mask->fmt != '1' && mask->fmt != '6') return -1;
  if (mask->n != qb->q->n) return -1;

  uint64_t n = qb->q->n;
  uint64_t nb = mask->fmt == '6' ? (n + 3) / 4 : (n + 7) / 8;

  for (uint64_t w = 0; w < qb->n_words; ++w) {
    uint64_t base = w * WBITS;
    uint64_t mw = (mask->fmt == '6') ? load6(mask->s, nb, base)
                                     : load0(mask->s, nb, base / 8);
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
  return 0;
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

  for (uint32_t s = 0; s < n_acc; ++s) {
    if (!b.bm[s]) continue;
    for (uint64_t w = 0; w < qb.n_words; ++w) {
      uint64_t mw = b.bm[s][w] & tail_mask(query->n, w);
      if (!mw) continue;
      acc[s].n_m += (uint64_t) pc64(mw);
      uint64_t ov = mw & qb.cov[w];
      if (!ov) continue;
      acc[s].n_o += (uint64_t) pc64(ov);
      uint64_t base = w * WBITS;
      while (ov) {
        int bit = ctz64(ov);
        uint64_t mu = f3_get_mu((cdata_t *) query, base + (uint64_t) bit);
        acc[s].sum_depth += MU2cov(mu);
        acc[s].sum_beta  += MU2beta(mu);
        ov &= ov - 1;
      }
    }
    free(b.bm[s]);
  }
  free(b.bm);
  if (n_q_out) *n_q_out = qb.n_q;
  yame_qbits_free(&qb);
  return 0;
}
