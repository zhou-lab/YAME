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
