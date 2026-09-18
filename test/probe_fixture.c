/* The .cm-exported-from-.mrmp fixture: runs against records, on the shape
 * methscope actually produces.
 *
 * methscope feeds the kernel membership RUNS walked out of a .mrmp; yame reads
 * the exported .cm as RECORDS. Nothing else enforces that the two views of one
 * membership give one answer. Here the .cm is decoded into runs (one slot per
 * state, the Pna state skipped as methscope skips it) and fed to the enumerator
 * path and to the inverted index, and each accumulator is compared bit for bit
 * with the record path over a binary mask built from the same codes.
 *
 * Usage: probe_fixture <mask.cm> <query.cg>
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cdata.h"
#include "cfile.h"
#include "summary.h"
#include "summary_multi.h"
#include "summary_index.h"

void prepare_mask(cdata_t *c);     /* summary.c */

static int fails = 0;
#define CHECK(c, ...) do { if (!(c)) { printf("  FAIL %d: ", __LINE__); \
  printf(__VA_ARGS__); printf("\n"); fails++; } } while (0)

typedef struct {
  cdata_t  *m;          /* the fmt2 records, inflated */
  uint32_t  n;
  uint32_t *base;       /* first slot of each record */
  uint32_t *na;         /* the Pna code of each record, or UINT32_MAX */
} bank_t;

/* Consecutive rows of one code are one run. Pna is not emitted: it marks rows
 * the set does not claim, and methscope's walker never emits it. */
static int emit_bank(void *ctx, yame_emit_fn emit, void *ec) {
  bank_t *b = ctx;
  for (uint32_t r = 0; r < b->n; ++r) {
    cdata_t *m = &b->m[r];
    uint64_t start = 0;
    uint64_t code = f2_get_uint64(m, 0);
    for (uint64_t i = 1; i <= m->n; ++i) {
      uint64_t c = (i < m->n) ? f2_get_uint64(m, i) : ~0ull;
      if (c == code) continue;
      if (code != b->na[r]) emit(ec, r, start, i - start, (uint32_t) code);
      start = i; code = c;
    }
  }
  return 0;
}

int main(int argc, char **argv) {
  if (argc < 3) { fprintf(stderr, "usage: probe_fixture mask.cm query.cg\n"); return 2; }

  /* the mask records */
  bank_t b = {0};
  cfile_t cf = open_cfile(argv[1]);
  uint32_t cap = 8; b.m = calloc(cap, sizeof(cdata_t));
  for (;;) {
    cdata_t m = read_cdata1(&cf);
    if (m.n == 0) break;
    if (b.n == cap) { cap *= 2; b.m = realloc(b.m, cap * sizeof(cdata_t)); }
    prepare_mask(&m);
    if (!m.aux) fmt2_set_aux(&m);
    b.m[b.n++] = m;
  }
  bgzf_close(cf.fh);
  CHECK(b.n > 0, "no record in %s", argv[1]);
  b.base = calloc(b.n + 1, sizeof(uint32_t));
  b.na = calloc(b.n, sizeof(uint32_t));
  for (uint32_t r = 0; r < b.n; ++r) {
    f2_aux_t *aux = b.m[r].aux;
    CHECK(b.m[r].fmt == '2', "record %u is format %c, want 2", r, b.m[r].fmt);
    b.base[r + 1] = b.base[r] + (uint32_t) aux->nk;
    b.na[r] = UINT32_MAX;
    for (uint64_t k = 0; k < aux->nk; ++k)
      if (strcmp(aux->keys[k], "Pna") == 0) b.na[r] = (uint32_t) k;
  }
  uint32_t n_slots = b.base[b.n];

  /* each slot as a binary mask, for the record path */
  uint64_t n = b.m[0].n, nb = (n + 7) / 8;
  cdata_t *bm = calloc(n_slots, sizeof(cdata_t));
  for (uint32_t r = 0; r < b.n; ++r) {
    for (uint32_t k = 0; k < b.base[r + 1] - b.base[r]; ++k) {
      cdata_t *x = &bm[b.base[r] + k];
      x->fmt = '0'; x->n = n; x->unit = 1; x->s = calloc(nb, 1);
    }
    for (uint64_t i = 0; i < n; ++i) {
      uint64_t c = f2_get_uint64(&b.m[r], i);
      if (c == b.na[r]) continue;
      uint8_t *s = bm[b.base[r] + c].s;
      s[i >> 3] |= (uint8_t) (1u << (i & 7));
    }
  }

  /* the index from the same runs, built once */
  char why[256];
  yame_index_t *ix = yame_index_build_runs(emit_bank, &b, n, b.n, b.base,
                                           n_slots, 0, why, sizeof why);
  CHECK(ix != NULL, "the index declined the fixture's runs: %s", why);

  /* every query record: enumerator, index and record path must agree */
  yame_acc_t *acc = calloc(n_slots, sizeof(yame_acc_t));
  yame_acc_t *aix = calloc(n_slots, sizeof(yame_acc_t));
  cf = open_cfile(argv[2]);
  uint32_t nq = 0, compared = 0;
  for (;; ++nq) {
    cdata_t q = read_cdata1(&cf);
    if (q.n == 0) break;
    prepare_mask(&q);
    CHECK(q.fmt == '3' && q.n == n, "query %u: format %c, %" PRIu64 " rows", nq, q.fmt, q.n);
    uint64_t n_q = 0;
    CHECK(yame_summarize_multi(&q, emit_bank, &b, b.n, b.base, acc, n_slots, &n_q) == 0,
          "query %u: the enumerator path refused", nq);
    if (ix) CHECK(yame_index_apply(ix, &q, aix) == 0, "query %u: the index refused", nq);
    yame_qbits_t qb;
    CHECK(yame_qbits_build(&q, &qb) == 0, "query %u: qbits refused", nq);
    for (uint32_t s = 0; s < n_slots; ++s) {
      uint32_t r = 0; while (b.base[r + 1] <= s) ++r;
      if (s - b.base[r] == b.na[r]) continue;            /* Pna: never claimed */
      yame_acc_t ref = {0};
      CHECK(yame_summarize_one_cx(&qb, &bm[s], &ref) == 0, "slot %u: record path refused", s);
      CHECK(acc[s].n_m == ref.n_m, "q%u slot %u n_m %" PRIu64 " vs %" PRIu64, nq, s, acc[s].n_m, ref.n_m);
      CHECK(acc[s].n_o == ref.n_o, "q%u slot %u n_o %" PRIu64 " vs %" PRIu64, nq, s, acc[s].n_o, ref.n_o);
      CHECK(acc[s].sum_depth == ref.sum_depth, "q%u slot %u sum_depth differs", nq, s);
      CHECK(acc[s].sum_beta == ref.sum_beta, "q%u slot %u sum_beta %.17g vs %.17g (runs vs records)",
            nq, s, acc[s].sum_beta, ref.sum_beta);
      if (ix) {
        CHECK(aix[s].n_m == ref.n_m && aix[s].n_o == ref.n_o && aix[s].sum_beta == ref.sum_beta,
              "q%u slot %u: the index differs from the record path", nq, s);
      }
      ++compared;
    }
    yame_qbits_free(&qb);
    free_cdata(&q);
  }
  bgzf_close(cf.fh);
  CHECK(nq > 0, "no record in %s", argv[2]);

  for (uint32_t s = 0; s < n_slots; ++s) free(bm[s].s);
  free(bm); free(acc); free(aix); free(b.base); free(b.na);
  for (uint32_t r = 0; r < b.n; ++r) free_cdata(&b.m[r]);
  free(b.m); yame_index_free(ix);
  if (fails) { printf("%d fixture assertion(s) failed\n", fails); return 1; }
  printf("ok: %u records x %u slots: runs, index and records agree bit for bit (%u pairs)\n",
         nq, n_slots, compared);
  return 0;
}
