/* The inverted index against the walk: bit-identical, from runs and from
 * records, and every way the build declines.
 *
 * summary_index.h is the shape methscope's featurizer uses (a bank of sets
 * against thousands of cells), fed from the same enumerator as the kernel. It
 * has to give the kernel's numbers exactly, refuse what it cannot hold, and say
 * what it would have needed.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cdata.h"
#include "summary.h"
#include "summary_multi.h"
#include "summary_index.h"

static int fails = 0;
#define CHECK(c, ...) do { if (!(c)) { printf("  FAIL %d: ", __LINE__); \
  printf(__VA_ARGS__); printf("\n"); fails++; } } while (0)

static cdata_t make_q3(uint64_t n, unsigned seed) {
  uint64_t *mu = malloc(n * sizeof(uint64_t));
  for (uint64_t i = 0; i < n; ++i) {
    unsigned r = (seed + (unsigned) i * 2654435761u);
    if ((r >> 3) % 5 == 0) { mu[i] = 0; continue; }
    mu[i] = ((uint64_t) ((r >> 7) % 9) << 32) | ((r >> 11) % 7 + 1);
  }
  cdata_t c = {0}; c.fmt='3'; c.n=n; c.unit=8; c.compressed=0; c.s=(uint8_t*)mu;
  return c;
}

/* Runs for NM masks x NS states, sparse the way a bank is. Slots of different
 * masks overlap; runs of one slot never do, which is the one obligation. */
#define NM 40u
#define NS 8u
static int emit_runs(void *ctx, yame_emit_fn emit, void *ec) {
  uint64_t n = *(uint64_t *) ctx;
  for (uint32_t m = 0; m < NM; ++m)
    for (uint32_t st = 0; st < NS; ++st)
      for (uint64_t start = (m * 37u + st * 11u) % 500; start + 4 < n; start += 500)
        emit(ec, m, start, 3, st);
  return 0;
}

/* The same runs, one slot's runs in DESCENDING order. The kernel needs
 * ascending; the index does not, and the header says so. */
static int emit_runs_desc(void *ctx, yame_emit_fn emit, void *ec) {
  uint64_t n = *(uint64_t *) ctx;
  for (uint32_t m = 0; m < NM; ++m)
    for (uint32_t st = 0; st < NS; ++st) {
      uint64_t first = (m * 37u + st * 11u) % 500, last = first;
      while (last + 500 + 4 < n) last += 500;
      for (uint64_t start = last; ; start -= 500) {
        if (start + 4 < n) emit(ec, m, start, 3, st);
        if (start == first) break;
      }
    }
  return 0;
}

/* A second enumeration that emits MORE than the first. */
static int calls = 0;
static int emit_grows(void *ctx, yame_emit_fn emit, void *ec) {
  uint64_t n = *(uint64_t *) ctx;
  emit(ec, 0, 10, 5, 0);
  if (calls++) emit(ec, 0, 100, 5, 0);
  (void) n;
  return 0;
}

static int emit_bad_mask(void *ctx, yame_emit_fn emit, void *ec) {
  (void) ctx; emit(ec, NM + 3, 0, 5, 0); return 0;
}

static int emit_wrap(void *ctx, yame_emit_fn emit, void *ec) {
  (void) ctx; emit(ec, 0, 1, UINT64_MAX, 0); return 0;
}

int main(void) {
  uint64_t n = 40000;
  cdata_t q = make_q3(n, 4242);
  uint32_t base[NM];
  for (uint32_t m = 0; m < NM; ++m) base[m] = m * NS;
  uint32_t n_acc = NM * NS;                     /* 320 slots */
  char why[256];

  /* ---- the kernel's answer, from the same runs ---------------------------- */
  yame_acc_t *ref = calloc(n_acc, sizeof(yame_acc_t));
  uint64_t n_q = 0;
  CHECK(yame_summarize_multi(&q, emit_runs, &n, NM, base, ref, n_acc, &n_q) == 0,
        "the kernel refused the runs");

  /* ---- the index, from runs, in the same order ---------------------------- */
  yame_index_t *ix = yame_index_build_runs(emit_runs, &n, n, NM, base, n_acc,
                                           0, why, sizeof why);
  CHECK(ix != NULL, "the runs-fed build declined: %s", why);
  if (!ix) return 1;
  CHECK(ix->n_slots == n_acc, "n_slots %u, want %u", ix->n_slots, n_acc);
  CHECK(ix->bytes == yame_index_bytes(n, ix->off[n]),
        "bytes %" PRIu64 " disagrees with yame_index_bytes", ix->bytes);

  yame_acc_t *acc = calloc(n_acc, sizeof(yame_acc_t));
  CHECK(yame_index_apply(ix, &q, acc) == 0, "apply refused a format 3 query");
  for (uint32_t s = 0; s < n_acc; ++s) {
    CHECK(acc[s].n_q == n_q, "slot %u n_q %" PRIu64 " vs %" PRIu64, s, acc[s].n_q, n_q);
    CHECK(acc[s].n_m == ref[s].n_m, "slot %u n_m %" PRIu64 " vs %" PRIu64, s, acc[s].n_m, ref[s].n_m);
    CHECK(acc[s].n_o == ref[s].n_o, "slot %u n_o %" PRIu64 " vs %" PRIu64, s, acc[s].n_o, ref[s].n_o);
    CHECK(acc[s].sum_depth == ref[s].sum_depth, "slot %u sum_depth differs", s);
    /* == not a tolerance: same doubles, same row order */
    CHECK(acc[s].sum_beta == ref[s].sum_beta, "slot %u sum_beta %.17g vs %.17g",
          s, acc[s].sum_beta, ref[s].sum_beta);
  }

  /* ---- driven by the coverage bitmap: the same numbers, covered rows only -- */
  {
    yame_qbits_t qb;
    CHECK(yame_qbits_build(&q, &qb) == 0, "qbits_build refused the query");
    yame_acc_t *a3 = calloc(n_acc, sizeof(yame_acc_t));
    CHECK(yame_index_apply_qb(ix, &qb, a3) == 0, "apply_qb refused");
    for (uint32_t s = 0; s < n_acc; ++s)
      CHECK(a3[s].n_q == acc[s].n_q && a3[s].n_o == acc[s].n_o &&
            a3[s].sum_depth == acc[s].sum_depth && a3[s].sum_beta == acc[s].sum_beta,
            "slot %u differs between apply and apply_qb", s);
    /* The covered-row LIST, the third route in: same numbers as apply, and a
     * list that overstates (every row) still gives them; a repeat or a wrong
     * order is refused rather than summed differently. */
    uint32_t *rows = malloc(n * sizeof(uint32_t)), *all = malloc(n * sizeof(uint32_t));
    uint64_t nr = 0;
    for (uint64_t i = 0; i < n; ++i) { all[i] = (uint32_t) i; if (f3_get_mu(&q, i)) rows[nr++] = (uint32_t) i; }
    yame_acc_t *a4 = calloc(n_acc, sizeof(yame_acc_t));
    CHECK(yame_index_apply_rows(ix, &q, rows, nr, a4) == 0, "apply_rows refused the list");
    for (uint32_t s = 0; s < n_acc; ++s)
      CHECK(a4[s].n_q == acc[s].n_q && a4[s].n_o == acc[s].n_o &&
            a4[s].sum_depth == acc[s].sum_depth && a4[s].sum_beta == acc[s].sum_beta,
            "slot %u differs between apply and apply_rows", s);
    CHECK(yame_index_apply_rows(ix, &q, all, n, a4) == 0, "apply_rows refused the full list");
    for (uint32_t s = 0; s < n_acc; ++s)
      CHECK(a4[s].n_q == acc[s].n_q && a4[s].sum_beta == acc[s].sum_beta,
            "slot %u differs for an overstated list", s);
    CHECK(yame_index_apply_rows(ix, &q, rows, 0, a4) == 0, "an empty list was refused");
    CHECK(a4[0].n_q == 0 && a4[0].n_o == 0, "an empty list counted something");
    { uint32_t bad[3] = { rows[0], rows[2], rows[1] };            /* out of order */
      CHECK(yame_index_apply_rows(ix, &q, bad, 3, a4) == -1, "an unordered list was accepted"); }
    { uint32_t bad[3] = { rows[0], rows[1], rows[1] };            /* a repeat */
      CHECK(yame_index_apply_rows(ix, &q, bad, 3, a4) == -1, "a repeated row was accepted"); }
    { uint32_t bad[1] = { (uint32_t) n };                          /* past the index */
      CHECK(yame_index_apply_rows(ix, &q, bad, 1, a4) == -1, "a row past the index was accepted"); }
    free(rows); free(all); free(a4);

    /* The query above is unit 8, the width a caller's own records have. The
     * same values at unit 1 (M and U fit four bits here) must read the same
     * through the accessor's other fast path. */
    cdata_t q1 = q; q1.unit = 1; q1.s = malloc(n);
    for (uint64_t i = 0; i < n; ++i) {
      uint64_t mu = f3_get_mu(&q, i);
      q1.s[i] = (uint8_t) (((mu >> 32) << 4) | (mu & 0xf));
    }
    yame_acc_t *a1 = calloc(n_acc, sizeof(yame_acc_t));
    CHECK(yame_index_apply(ix, &q1, a1) == 0, "apply refused unit 1");
    for (uint32_t s = 0; s < n_acc; ++s)
      CHECK(a1[s].sum_beta == acc[s].sum_beta && a1[s].n_o == acc[s].n_o &&
            a1[s].sum_depth == acc[s].sum_depth,
            "slot %u differs between unit 8 and unit 1", s);
    free(a3); free(a1); free(q1.s); yame_qbits_free(&qb);
  }

  /* ---- runs of one slot in DESCENDING order: still bit-identical ---------- */
  {
    yame_index_t *ix2 = yame_index_build_runs(emit_runs_desc, &n, n, NM, base,
                                              n_acc, 0, why, sizeof why);
    CHECK(ix2 != NULL, "descending runs declined: %s", why);
    if (ix2) {
      yame_acc_t *a2 = calloc(n_acc, sizeof(yame_acc_t));
      CHECK(yame_index_apply(ix2, &q, a2) == 0, "apply refused");
      for (uint32_t s = 0; s < n_acc; ++s)
        CHECK(a2[s].sum_beta == ref[s].sum_beta && a2[s].n_o == ref[s].n_o,
              "slot %u differs when its runs arrive descending", s);
      free(a2); yame_index_free(ix2);
    }
  }

  /* ---- what apply refuses --------------------------------------------------- */
  {
    cdata_t cq = q; cq.compressed = 1;
    CHECK(yame_index_apply(ix, &cq, acc) == -1, "a compressed query was accepted");
    cdata_t nq = q; nq.n = n - 1;
    CHECK(yame_index_apply(ix, &nq, acc) == -1, "a query of other rows was accepted");
    cdata_t q6 = q; q6.fmt = '6';
    CHECK(yame_index_apply(ix, &q6, acc) == -1, "a format 6 query was accepted");
  }

  /* ---- how the build declines, and what it says ---------------------------- */
  {
    /* over budget: the message names the size it wanted and the budget */
    uint64_t need = ix->bytes;
    yame_index_t *b = yame_index_build_runs(emit_runs, &n, n, NM, base, n_acc,
                                            need - 1, why, sizeof why);
    CHECK(b == NULL, "a build one byte over budget went ahead");
    CHECK(strstr(why, "index needs") && strstr(why, "budget"),
          "the decline did not say what it needed: %s", why);
    char want[64]; snprintf(want, sizeof want, "%" PRIu32 " memberships", ix->off[n]);
    CHECK(strstr(why, want) != NULL,
          "the decline's membership count is not the exact one: %s", why);
    /* exactly at budget: fits */
    b = yame_index_build_runs(emit_runs, &n, n, NM, base, n_acc, need, why, sizeof why);
    CHECK(b != NULL, "a build exactly at budget declined: %s", why);
    yame_index_free(b);
    /* the offsets alone over budget: declined before anything is enumerated */
    calls = 0;
    b = yame_index_build_runs(emit_grows, &n, n, 1, NULL, 1, 8, why, sizeof why);
    CHECK(b == NULL && calls == 0, "the offsets check did not run before the pass");
    CHECK(strstr(why, "offsets alone") != NULL, "unexpected reason: %s", why);
  }
  {
    /* too many slots for a 16-bit entry */
    yame_index_t *b = yame_index_build_runs(emit_runs, &n, n, NM, base, 65536,
                                            0, why, sizeof why);
    CHECK(b == NULL && strstr(why, "65,535"), "65,536 slots were accepted: %s", why);
  }
  {
    /* a run naming a mask outside the layout */
    yame_index_t *b = yame_index_build_runs(emit_bad_mask, &n, n, NM, base, n_acc,
                                            0, why, sizeof why);
    CHECK(b == NULL && strstr(why, "outside"), "a run outside the layout was indexed: %s", why);
  }
  {
    /* a second enumeration that emits more than the first must fail, not
     * write past the entries */
    calls = 0;
    yame_index_t *b = yame_index_build_runs(emit_grows, &n, n, 1, NULL, 1,
                                            0, why, sizeof why);
    CHECK(b == NULL && strstr(why, "second"), "a growing enumeration was accepted: %s", why);
  }
  {
    /* a run whose length wraps start + len is clipped, not indexed past n */
    yame_index_t *b = yame_index_build_runs(emit_wrap, &n, n, 1, NULL, 1,
                                            0, why, sizeof why);
    CHECK(b != NULL, "the wrapping run failed the build: %s", why);
    if (b) {
      CHECK(b->n_m[0] == n - 1, "a wrapping run claimed %" PRIu64 " rows, want %" PRIu64,
            b->n_m[0], n - 1);
      yame_index_free(b);
    }
  }

  yame_index_free(ix);
  free(acc); free(ref); free(q.s);
  if (fails) { printf("%d index assertion(s) failed\n", fails); return 1; }
  printf("ok: index from runs == kernel, bit for bit, over %u slots\n", n_acc);
  return 0;
}
