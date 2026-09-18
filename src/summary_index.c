/* summary_index.c -- see summary_index.h. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdarg.h>
#include "summary_index.h"
#include "cfile.h"

void prepare_mask(cdata_t *c);     /* summary.c */

static double now_s(void) {
  struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t);
  return t.tv_sec + 1e-9 * t.tv_nsec;
}

static void say(char *why, size_t why_len, const char *fmt, ...) {
  if (!why || !why_len) return;
  va_list ap; va_start(ap, fmt);
  vsnprintf(why, why_len, fmt, ap);
  va_end(ap);
}

/* Bytes as megabytes for a sentence: whole above 10 MB, two decimals below. */
static const char *mb(uint64_t bytes, char buf[32]) {
  double v = bytes / 1048576.0;
  snprintf(buf, 32, v >= 10 ? "%.0f" : "%.2f", v);
  return buf;
}

uint64_t yame_index_bytes(uint64_t n_rows, uint64_t n_memberships) {
  return 4 * (n_rows + 1) + 2 * n_memberships;
}

/* --------------------------------------------------- the shared skeleton -- */

/* Both builds fill the same three arrays in the same two passes; only where
 * the memberships come from differs. `cnt` is one uint32 per row: the count in
 * pass 1, the write cursor in pass 2. */
typedef struct {
  uint64_t  n_rows;
  uint32_t *cnt;
  uint32_t *off;
  uint16_t *ent;
  uint64_t *n_m;
  uint64_t  total;       /* memberships seen so far                         */
  uint64_t  budget;      /* 0: none                                         */
  int       over;        /* pass 1 crossed the budget: stop scattering      */
  int       bad;         /* a run outside the layout, or pass 2 overran     */
} build_t;

static int over_budget(const build_t *b) {
  return b->budget && yame_index_bytes(b->n_rows, b->total) > b->budget;
}

/* Turn the counts into offsets and reset cnt as the pass-2 cursor. */
static int prefix(build_t *b) {
  b->off = malloc((b->n_rows + 1) * sizeof(uint32_t));
  b->ent = malloc((b->total ? b->total : 1) * sizeof(uint16_t));
  if (!b->off || !b->ent) return -1;
  uint64_t acc = 0;
  for (uint64_t i = 0; i < b->n_rows; ++i) {
    b->off[i] = (uint32_t) acc; acc += b->cnt[i];
  }
  b->off[b->n_rows] = (uint32_t) acc;
  memset(b->cnt, 0, (b->n_rows + 1) * sizeof(uint32_t));
  return 0;
}

static yame_index_t *finish(build_t *b, uint32_t n_slots, double t0) {
  yame_index_t *ix = calloc(1, sizeof(yame_index_t));
  if (!ix) return NULL;
  ix->n_rows = b->n_rows; ix->n_slots = n_slots;
  ix->off = b->off; ix->ent = b->ent; ix->n_m = b->n_m;
  ix->bytes = yame_index_bytes(b->n_rows, b->total);
  ix->t_build = now_s() - t0;
  return ix;
}

static void scrap(build_t *b) {
  free(b->cnt); free(b->off); free(b->ent); free(b->n_m);
}

/* ------------------------------------------------------- from a mask file -- */

#define WBITS 64
#if defined(__GNUC__) || defined(__clang__)
#  define pc64(x)  ((unsigned) __builtin_popcountll(x))
#  define ctz64(x) ((int) __builtin_ctzll(x))
#else
static inline unsigned pc64(uint64_t x) { unsigned c = 0; for (; x; x &= x - 1) ++c; return c; }
static inline int ctz64(uint64_t x) { int n = 0; while (!(x & 1)) { x >>= 1; ++n; } return n; }
#endif

/* The word at bit offset w*64 of a bit-packed mask of n rows, without reading
 * past its last byte. */
static inline uint64_t mask_word(const uint8_t *s, uint64_t n, uint64_t w) {
  uint64_t byte = w * 8, nb = (n + 7) >> 3, x = 0;
  if (byte + 8 <= nb) memcpy(&x, s + byte, 8);
  else memcpy(&x, s + byte, nb - byte);
  uint64_t base = w * WBITS;
  if (n - base < WBITS) x &= (1ull << (n - base)) - 1;   /* past the last row */
  return x;
}

/* One fmt0 mask, 64 rows at a time. Pass 1 (fill == 0) counts per row unless
 * the build is already over budget, in which case only the popcount is kept,
 * which is what makes the decline message exact. Pass 2 writes the slot id at
 * each row's cursor. Returns the rows claimed. */
static uint64_t walk_bits(build_t *b, const uint8_t *s, int fill, uint32_t id) {
  uint64_t n = b->n_rows, nw = (n + WBITS - 1) / WBITS, claimed = 0;
  for (uint64_t w = 0; w < nw; ++w) {
    uint64_t x = mask_word(s, n, w);
    if (!x) continue;
    claimed += pc64(x);
    if (!fill && b->over) continue;
    uint64_t base = w * WBITS;
    while (x) {
      uint64_t i = base + (uint64_t) ctz64(x);
      if (fill) b->ent[b->off[i] + b->cnt[i]++] = (uint16_t) id;
      else ++b->cnt[i];
      x &= x - 1;
    }
  }
  return claimed;
}

yame_index_t *yame_index_build(const char *mask_path, uint64_t n_rows,
                               uint64_t budget_bytes, char *why, size_t why_len) {
  double t0 = now_s();
  char m1[32], m2[32];
  if (why && why_len) why[0] = 0;
  if (budget_bytes && yame_index_bytes(n_rows, 0) > budget_bytes) {
    say(why, why_len, "the offsets alone need %s MB over %"PRIu64" rows, "
        "budget %s MB", mb(yame_index_bytes(n_rows, 0), m1), n_rows,
        mb(budget_bytes, m2));
    return NULL;
  }
  build_t b = {0};
  b.n_rows = n_rows; b.budget = budget_bytes;
  b.cnt = calloc(n_rows + 1, sizeof(uint32_t));
  if (!b.cnt) { say(why, why_len, "out of memory for the row counts"); return NULL; }

  /* pass 1: how many slots claim each row, and how many rows each slot claims */
  uint64_t cap = 1024, n_masks = 0;
  b.n_m = malloc(cap * sizeof(uint64_t));
  cfile_t cf = open_cfile((char *) mask_path);
  for (;;) {
    cdata_t m = read_cdata1(&cf);
    if (m.n == 0) break;
    if (m.fmt > '1') {
      say(why, why_len, "record %"PRIu64" is format %c, not a binary mask; a "
          "state mask claims every row", n_masks + 1, m.fmt);
      free_cdata(&m); goto decline;
    }
    prepare_mask(&m);
    if (m.n != n_rows) {
      say(why, why_len, "record %"PRIu64" has %"PRIu64" rows, the query "
          "%"PRIu64, n_masks + 1, m.n, n_rows);
      free_cdata(&m); goto decline;
    }
    if (n_masks >= 65535) {
      say(why, why_len, "more than 65,535 records; an entry is 16 bits");
      free_cdata(&m); goto decline;
    }
    if (n_masks == cap) {
      cap *= 2; b.n_m = realloc(b.n_m, cap * sizeof(uint64_t));
    }
    b.n_m[n_masks] = walk_bits(&b, m.s, 0, 0);
    b.total += b.n_m[n_masks];
    ++n_masks;
    free_cdata(&m);
    if (!b.over && over_budget(&b)) b.over = 1;
  }
  bgzf_close(cf.fh);
  if (!n_masks) { say(why, why_len, "the mask file holds no record"); goto decline0; }
  if (b.over) {
    say(why, why_len, "index needs %s MB (%"PRIu64" masks, %"PRIu64" "
        "memberships), budget %s MB", mb(yame_index_bytes(n_rows, b.total), m1),
        n_masks, b.total, mb(budget_bytes, m2));
    goto decline0;
  }
  if (b.total > 0xFFFFFFFFull) {
    say(why, why_len, "%"PRIu64" memberships; an offset is 32 bits", b.total);
    goto decline0;
  }
  if (prefix(&b) != 0) { say(why, why_len, "out of memory for the index"); goto decline0; }

  /* pass 2: write the slot ids grouped by row */
  cf = open_cfile((char *) mask_path);
  for (uint64_t k = 0; k < n_masks; ++k) {
    cdata_t m = read_cdata1(&cf);
    if (m.n == 0) break;
    prepare_mask(&m);
    walk_bits(&b, m.s, 1, (uint32_t) k);
    free_cdata(&m);
  }
  bgzf_close(cf.fh);
  free(b.cnt); b.cnt = NULL;
  return finish(&b, (uint32_t) n_masks, t0);

 decline:
  bgzf_close(cf.fh);
 decline0:
  scrap(&b);
  return NULL;
}

/* ------------------------------------------------------------- from runs -- */

typedef struct {
  build_t  *b;
  uint32_t  n_masks;
  const uint32_t *state_base;
  uint32_t  n_slots;
  int       fill;
} runs_ctx_t;

static void index_emit(void *ctx, uint32_t mask,
                       uint64_t start, uint64_t len, uint32_t state) {
  runs_ctx_t *r = ctx;
  build_t *b = r->b;
  if (!len || b->bad) return;
  if (mask >= r->n_masks) { b->bad = 1; return; }
  uint32_t slot = r->state_base ? r->state_base[mask] + state : mask;
  if (slot >= r->n_slots) { b->bad = 1; return; }
  uint64_t n = b->n_rows;
  if (start >= n) return;
  uint64_t end = (len > n - start) ? n : start + len;   /* clamp before the add */
  if (!r->fill) {
    b->n_m[slot] += end - start;
    b->total += end - start;
    if (b->over) return;
    for (uint64_t i = start; i < end; ++i) ++b->cnt[i];
    if (over_budget(b)) b->over = 1;
  } else {
    for (uint64_t i = start; i < end; ++i) {
      if (b->off[i] + b->cnt[i] >= b->off[i + 1]) { b->bad = 1; return; }
      b->ent[b->off[i] + b->cnt[i]++] = (uint16_t) slot;
    }
  }
}

yame_index_t *yame_index_build_runs(yame_runs_fn runs, void *runs_ctx,
                                    uint64_t n_rows, uint32_t n_masks,
                                    const uint32_t *state_base, uint32_t n_slots,
                                    uint64_t budget_bytes,
                                    char *why, size_t why_len) {
  double t0 = now_s();
  char m1[32], m2[32];
  if (why && why_len) why[0] = 0;
  if (!n_slots || n_slots > 65535) {
    say(why, why_len, "%u slots; an entry is 16 bits, so at most 65,535", n_slots);
    return NULL;
  }
  if (budget_bytes && yame_index_bytes(n_rows, 0) > budget_bytes) {
    say(why, why_len, "the offsets alone need %s MB over %"PRIu64" rows, "
        "budget %s MB", mb(yame_index_bytes(n_rows, 0), m1), n_rows,
        mb(budget_bytes, m2));
    return NULL;
  }
  build_t b = {0};
  b.n_rows = n_rows; b.budget = budget_bytes;
  b.cnt = calloc(n_rows + 1, sizeof(uint32_t));
  b.n_m = calloc(n_slots, sizeof(uint64_t));
  if (!b.cnt || !b.n_m) { say(why, why_len, "out of memory for the row counts"); scrap(&b); return NULL; }

  runs_ctx_t r = { &b, n_masks, state_base, n_slots, 0 };
  if (runs(runs_ctx, index_emit, &r) != 0 || b.bad) {
    say(why, why_len, b.bad ? "a run named a mask or slot outside the layout"
                            : "the enumerator failed");
    scrap(&b); return NULL;
  }
  if (b.over) {
    say(why, why_len, "index needs %s MB (%u slots, %"PRIu64" memberships), "
        "budget %s MB", mb(yame_index_bytes(n_rows, b.total), m1),
        n_slots, b.total, mb(budget_bytes, m2));
    scrap(&b); return NULL;
  }
  if (b.total > 0xFFFFFFFFull) {
    say(why, why_len, "%"PRIu64" memberships; an offset is 32 bits", b.total);
    scrap(&b); return NULL;
  }
  if (prefix(&b) != 0) { say(why, why_len, "out of memory for the index"); scrap(&b); return NULL; }

  r.fill = 1;
  if (runs(runs_ctx, index_emit, &r) != 0 || b.bad) {
    say(why, why_len, b.bad ? "the second enumeration emitted more than the first"
                            : "the enumerator failed");
    scrap(&b); return NULL;
  }
  free(b.cnt); b.cnt = NULL;
  return finish(&b, n_slots, t0);
}

void yame_index_free(yame_index_t *ix) {
  if (!ix) return;
  free(ix->off); free(ix->ent); free(ix->n_m); free(ix);
}

/* ---------------------------------------------------------------- apply -- */

int yame_index_apply(const yame_index_t *ix, const cdata_t *query,
                     yame_acc_t *acc) {
  if (query->fmt != '3' || query->compressed || query->n != ix->n_rows) return -1;
  memset(acc, 0, (size_t) ix->n_slots * sizeof(yame_acc_t));

  uint64_t n_q = 0;
  const uint32_t *off = ix->off;
  const uint16_t *ent = ix->ent;
  for (uint64_t i = 0; i < ix->n_rows; ++i) {
    uint64_t mu = f3_get_mu((cdata_t *) query, i);
    if (!mu) continue;
    ++n_q;
    uint32_t e = off[i], e1 = off[i + 1];
    if (e == e1) continue;                    /* no slot claims this row */
    double b = MU2beta(mu);
    uint64_t d = MU2cov(mu);
    for (; e < e1; ++e) {
      yame_acc_t *a = &acc[ent[e]];
      ++a->n_o; a->sum_beta += b; a->sum_depth += d;
    }
  }
  for (uint32_t k = 0; k < ix->n_slots; ++k) {
    acc[k].n_u = ix->n_rows;
    acc[k].n_q = n_q;
    acc[k].n_m = ix->n_m[k];
    acc[k].beta = acc[k].sum_beta / acc[k].n_o;
  }
  return 0;
}
