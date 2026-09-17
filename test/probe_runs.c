/* The enumerator path, which is what an MRMP bank uses, against the record path.
 *
 * It matters that this is checked separately: the two share the arithmetic but
 * not the route in, and the enumerator form is the one that has to survive
 * thousands of slots without allocating per slot.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cdata.h"
#include "summary.h"
#include "summary_multi.h"

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

/* Runs for NM masks x NS states, sparse the way a bank is: most rows in no set. */
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

int main(void) {
  uint64_t n = 40000;
  cdata_t q = make_q3(n, 4242);

  uint32_t base[NM];
  for (uint32_t m = 0; m < NM; ++m) base[m] = m * NS;
  uint32_t n_acc = NM * NS;                     /* 320 slots */

  yame_acc_t *acc = calloc(n_acc, sizeof(yame_acc_t));
  uint64_t n_q = 0;
  CHECK(yame_summarize_multi(&q, emit_runs, &n, NM, base, acc, n_acc, &n_q) == 0,
        "the enumerator path refused sparse runs over %u slots", n_acc);

  /* The same membership as records, one mask per slot, through the other path. */
  cdata_t *m = calloc(n_acc, sizeof(cdata_t));
  uint64_t nb = (n + 7) / 8;
  for (uint32_t mm = 0; mm < NM; ++mm)
    for (uint32_t st = 0; st < NS; ++st) {
      uint32_t slot = mm * NS + st;
      m[slot].fmt='0'; m[slot].n=n; m[slot].unit=1; m[slot].s=calloc(nb,1);
      for (uint64_t s = (mm * 37u + st * 11u) % 500; s + 4 < n; s += 500)
        for (uint64_t k = 0; k < 3; ++k) m[slot].s[(s+k)>>3] |= 1u << ((s+k)&7);
    }
  yame_acc_t *ref = calloc(n_acc, sizeof(yame_acc_t));
  uint64_t n_q2 = 0;
  CHECK(yame_summarize_multi_cx(&q, m, n_acc, ref, &n_q2) == 0, "the record path refused");
  CHECK(n_q == n_q2, "covered counts differ: %" PRIu64 " vs %" PRIu64, n_q, n_q2);

  uint64_t nonzero = 0;
  for (uint32_t s = 0; s < n_acc; ++s) {
    CHECK(acc[s].n_m == ref[s].n_m, "slot %u n_m %" PRIu64 " vs %" PRIu64, s, acc[s].n_m, ref[s].n_m);
    CHECK(acc[s].n_o == ref[s].n_o, "slot %u n_o %" PRIu64 " vs %" PRIu64, s, acc[s].n_o, ref[s].n_o);
    CHECK(acc[s].sum_beta == ref[s].sum_beta, "slot %u sum_beta %.17g vs %.17g",
          s, acc[s].sum_beta, ref[s].sum_beta);
    if (acc[s].n_o) nonzero++;
  }
  CHECK(nonzero > NM, "only %" PRIu64 " of %u slots got any overlap", nonzero, n_acc);

  for (uint32_t s = 0; s < n_acc; ++s) free(m[s].s);
  free(m); free(acc); free(ref); free(q.s);
  if (fails) { printf("%d enumerator assertion(s) failed\n", fails); return 1; }
  /* What the first design would have cost here, for the record: it allocated a
   * bitmap per slot, n/8 bytes each. At a whole-genome row space and a bank's
   * slot count that is gigabytes, which is why runs are now scored where they
   * arrive. */
  printf("  enumerator path: %u slots, agrees with records"
         " (per-slot bitmaps would have been %.1f MB here,"
         " %.1f GB at 29.4M rows)\n", n_acc,
         (double) n_acc * ((n + 7) / 8) / 1048576.0,
         (double) n_acc * ((29401795 + 7) / 8) / 1073741824.0);
  return 0;
}
