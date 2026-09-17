/* The multi-mask kernel against summarize1(), bit for bit.
 *
 * "Bit-identical, not approximately equal" is the requirement: within a mask
 * the rows are accumulated in increasing order either way, so the same doubles
 * are added in the same sequence, and a mismatch means the kernel changed the
 * arithmetic rather than merely the schedule.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cdata.h"
#include "summary.h"
#include "summary_multi.h"

stats_t* summarize1(cdata_t *c, cdata_t *c_mask, uint64_t *n_st,
                    char *sm, char *sq, config_t *config);

static int fails = 0;
#define CHECK(cond, ...) do { if (!(cond)) { \
  printf("  FAIL %s:%d: ", __func__, __LINE__); printf(__VA_ARGS__); \
  printf("\n"); fails++; } } while (0)

/* A format 3 query: M/U per row, some rows uncovered. */
static cdata_t make_q3(uint64_t n, unsigned seed) {
  uint64_t *mu = malloc(n * sizeof(uint64_t));
  for (uint64_t i = 0; i < n; ++i) {
    unsigned r = (seed + (unsigned)i * 2654435761u);
    if ((r >> 3) % 5 == 0) { mu[i] = 0; continue; }     /* uncovered */
    uint64_t m = (r >> 7) % 9, u = (r >> 11) % 7 + 1;
    mu[i] = (m << 32) | u;
  }
  cdata_t c = {0};
  c.fmt = '3'; c.n = n; c.unit = 8; c.compressed = 0;
  c.s = (uint8_t *) mu;
  return c;
}

/* A format 6 record: two bits a row, universe at the odd bit and set at the
 * even one, with some rows outside the universe so the restriction is exercised. */
static cdata_t make_q6(uint64_t n, unsigned seed) {
  uint64_t nb = (n + 3) / 4;
  uint8_t *s = calloc(nb, 1);
  for (uint64_t i = 0; i < n; ++i) {
    unsigned r = (seed + (unsigned) i * 2246822519u);
    int uni = ((r >> 5) % 4) != 0;            /* a quarter outside the universe */
    int set = uni && ((r >> 9) % 3) == 0;
    if (set) s[i >> 2] |= (uint8_t) (1u << ((i & 3) * 2));
    if (uni) s[i >> 2] |= (uint8_t) (1u << ((i & 3) * 2 + 1));
  }
  cdata_t c = {0};
  c.fmt = '6'; c.n = n; c.unit = 2; c.compressed = 0; c.s = s;
  return c;
}

/* A binary mask over the same rows, `every`-th row a member, offset by phase,
 * which gives runs rather than single rows. */
static cdata_t make_mask0(uint64_t n, int every, int phase, int runlen) {
  uint64_t nb = (n + 7) / 8;
  uint8_t *s = calloc(nb, 1);
  for (uint64_t i = 0; i < n; ++i)
    if ((int)((i + phase) % every) < runlen) s[i >> 3] |= (1u << (i & 7));
  cdata_t c = {0};
  c.fmt = '0'; c.n = n; c.unit = 1; c.compressed = 0; c.s = s;
  return c;
}

int main(void) {
  const uint64_t N = 5000;
  cdata_t q = make_q3(N, 12345);

  const uint32_t NM = 6;
  cdata_t m[6];
  int every[6]  = { 3, 7, 11, 2, 97, 5 };
  int phase[6]  = { 0, 1,  2, 0,  3, 4 };
  int runlen[6] = { 1, 3,  5, 1, 40, 2 };
  for (uint32_t j = 0; j < NM; ++j) m[j] = make_mask0(N, every[j], phase[j], runlen[j]);

  /* the reference: one call per mask */
  config_t cfg; memset(&cfg, 0, sizeof cfg);
  stats_t ref[6]; memset(ref, 0, sizeof ref);
  for (uint32_t j = 0; j < NM; ++j) {
    uint64_t n_st = 0;
    stats_t *st = summarize1(&q, &m[j], &n_st, (char *)"m", (char *)"q", &cfg);
    CHECK(n_st == 1, "mask %u gave %" PRIu64 " stats", j, n_st);
    ref[j] = st[0];
    free(st[0].sm); free(st[0].sq); free(st);
  }

  /* the kernel: one pass */
  yame_acc_t acc[6]; memset(acc, 0, sizeof acc);
  uint64_t n_q = 0;
  int rc = yame_summarize_multi_cx(&q, m, NM, acc, &n_q);
  CHECK(rc == 0, "the kernel refused a format 3 query and binary masks (rc %d)", rc);

  for (uint32_t j = 0; j < NM; ++j) {
    CHECK(acc[j].n_m == ref[j].n_m, "mask %u n_m %" PRIu64 " vs %" PRIu64,
          j, acc[j].n_m, ref[j].n_m);
    CHECK(acc[j].n_o == ref[j].n_o, "mask %u n_o %" PRIu64 " vs %" PRIu64,
          j, acc[j].n_o, ref[j].n_o);
    CHECK(acc[j].sum_depth == ref[j].sum_depth, "mask %u depth %" PRIu64 " vs %" PRIu64,
          j, acc[j].sum_depth, ref[j].sum_depth);
    /* the point of the exercise: the same double, not a near one */
    CHECK(acc[j].sum_beta == ref[j].sum_beta,
          "mask %u sum_beta %.17g vs %.17g (差 %.3g)", j,
          acc[j].sum_beta, ref[j].sum_beta, acc[j].sum_beta - ref[j].sum_beta);
    CHECK(n_q == ref[j].n_q, "n_q %" PRIu64 " vs mask %u's %" PRIu64,
          n_q, j, ref[j].n_q);
  }

  /* ---- a format 6 QUERY, which measures something else entirely ----------
   * Its universe is the query's own, narrowed further by a format 6 mask's;
   * every count is taken inside it; and beta is the overlap over the mask
   * count, not a mean of betas. Getting that wrong would print plausible
   * numbers, so it is checked against summarize1() the same way. */
  cdata_t q6 = make_q6(N, 999);   /* seeded so the mix of in and out of universe is fixed */
  yame_acc_t a6[6]; memset(a6, 0, sizeof a6);
  uint64_t nq6 = 0;
  CHECK(yame_summarize_multi_cx(&q6, m, NM, a6, &nq6) == 0,
        "the kernel refused a format 6 query");
  for (uint32_t j = 0; j < NM; ++j) {
    uint64_t n_st = 0;
    stats_t *st = summarize1(&q6, &m[j], &n_st, (char *)"m", (char *)"q", &cfg);
    CHECK(a6[j].n_u == st[0].n_u, "fmt6 mask %u n_u %" PRIu64 " vs %" PRIu64,
          j, a6[j].n_u, st[0].n_u);
    CHECK(a6[j].n_q == st[0].n_q, "fmt6 mask %u n_q %" PRIu64 " vs %" PRIu64,
          j, a6[j].n_q, st[0].n_q);
    CHECK(a6[j].n_m == st[0].n_m, "fmt6 mask %u n_m %" PRIu64 " vs %" PRIu64,
          j, a6[j].n_m, st[0].n_m);
    CHECK(a6[j].n_o == st[0].n_o, "fmt6 mask %u n_o %" PRIu64 " vs %" PRIu64,
          j, a6[j].n_o, st[0].n_o);
    CHECK(a6[j].beta == st[0].beta, "fmt6 mask %u beta %.17g vs %.17g",
          j, a6[j].beta, st[0].beta);
    free(st[0].sm); free(st[0].sq); free(st);
  }

  /* a format 6 MASK against that query narrows the universe again */
  cdata_t m6 = make_q6(N, 4242);
  yame_acc_t am6 = {0};
  CHECK(yame_summarize_multi_cx(&q6, &m6, 1, &am6, NULL) == 0,
        "the kernel refused a format 6 mask");
  {
    uint64_t n_st = 0;
    stats_t *st = summarize1(&q6, &m6, &n_st, (char *)"m", (char *)"q", &cfg);
    CHECK(am6.n_u == st[0].n_u, "fmt6 mask narrowing: n_u %" PRIu64 " vs %" PRIu64,
          am6.n_u, st[0].n_u);
    CHECK(am6.n_m == st[0].n_m, "fmt6 mask narrowing: n_m %" PRIu64 " vs %" PRIu64,
          am6.n_m, st[0].n_m);
    CHECK(am6.n_o == st[0].n_o, "fmt6 mask narrowing: n_o %" PRIu64 " vs %" PRIu64,
          am6.n_o, st[0].n_o);
    CHECK(am6.beta == st[0].beta, "fmt6 mask narrowing: beta %.17g vs %.17g",
          am6.beta, st[0].beta);
    free(st[0].sm); free(st[0].sq); free(st);
  }
  free(m6.s); free(q6.s);

  /* a format the kernel does not cover must be refused, not guessed at */
  cdata_t q4 = q; q4.fmt = '4';
  yame_acc_t a1[6]; memset(a1, 0, sizeof a1);
  CHECK(yame_summarize_multi_cx(&q4, m, NM, a1, NULL) == -1,
        "the kernel accepted a format 4 query");

  /* a mask of the wrong length must be refused, so the per-mask path reports it */
  cdata_t short_mask = make_mask0(N / 2, 3, 0, 1);
  yame_acc_t a2[1]; memset(a2, 0, sizeof a2);
  CHECK(yame_summarize_multi_cx(&q, &short_mask, 1, a2, NULL) == -1,
        "the kernel accepted a mask shorter than the query");
  free(short_mask.s);

  for (uint32_t j = 0; j < NM; ++j) free(m[j].s);
  free(q.s);

  if (fails) { printf("%d multi-mask assertion(s) failed\n", fails); return 1; }
  return 0;
}
