// SPDX-License-Identifier: AGPL-3.0-or-later
/**
 * This file is part of YAME.
 *
 * Copyright (C) 2021-present Wanding Zhou
 *
 * YAME is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * YAME is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with YAME.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "cfile.h"

#include "yame_ui.h"
#include <math.h>
/**
 * yame pairwise
 * =============
 *
 * Goal
 * ----
 * Generate a differential-methylation "set" track (format 6) by comparing two
 * format-3 MU vectors site-by-site. The output can be used downstream for
 * overlap/enrichment or as a mask-like feature set.
 *
 * Inputs
 * ------
 * - Two cdata records with the same length N, each format 3 (M/U counts) or
 *   format 4 (betas, NA stored negative). A format-4 side is "covered" at a
 *   site when it holds a value; a format-3 side when M+U >= min_cov.
 * - If MU2.cx is not provided, the code reads the first two records from MU1.cx.
 * - -1/-2 name the record to take from each file through its index; without a
 *   name the first record is used.
 *
 * Summary mode (-S)
 * -----------------
 * The same row-aligned walk, but instead of writing a set it accumulates
 * agreement statistics over the sites BOTH sides cover and prints one TSV line
 * per pairing: n, MAE, RMSE, call agreement at a beta threshold (-t), Pearson r.
 * If file 2 has several records and none is named, every one of them is scored
 * against side 1 and the lines are printed in file order -- side 1 is
 * decompressed once, side 2 streams. Two multi-record files with neither side
 * named is an error, not a guess about which pairing was meant.
 *
 * Core definitions
 * ----------------
 * - mu = f3_get_mu(&c, i) returns packed (M,U) for site i.
 * - cov(mu)  = MU2cov(mu)  = M + U
 * - beta(mu) = MU2beta(mu) = M/(M+U)   (double)
 *
 * Universe rule (format 6)
 * ------------------------
 * A site is considered "valid/measured" only if BOTH samples meet minimum coverage:
 *
 *   if cov(mu1) >= min_coverage AND cov(mu2) >= min_coverage:
 *       site is in universe (FMT6 universe bit = 1)
 *       set bit is assigned by the direction/effect rule below
 *   else:
 *       site remains outside universe (universe bit stays 0; set bit irrelevant)
 *
 * Set rule (direction + effect size)
 * ----------------------------------
 * For sites passing the universe rule, the output set bit is decided by -H:
 *
 *   mode 1 (default): beta1 > beta2 AND (beta1 - beta2) > min_effect
 *   mode 2           : beta1 < beta2 AND (beta2 - beta1) > min_effect
 *   mode 3           : "different"
 *       - if min_effect <= 0: beta1 != beta2
 *       - else: |beta1 - beta2| > min_effect
 *
 * Output
 * ------
 * Writes a single compressed format-6 cdata record of length N to stdout or -o.
 * No index is generated (even with -o) in the current implementation.
 *
 * Notes / gotchas
 * ---------------
 * - min_effect defaults to 0 (after option parsing), meaning any non-zero beta
 *   difference can be flagged in mode 3, and only strict inequalities in modes 1/2.
 * - Comparisons use doubles; exact equality in mode 3 (min_effect <= 0) can be
 *   sensitive to floating representation (often fine here because beta derives
 *   from integer ratios, but still a consideration).
 */

static int usage(void) {
  yame_usage_head("yame pairwise [options] <MU1.cx> [MU2.cx] > out.cx");
  yame_usage_sec("Purpose:");
  yame_usage_text("Compare two samples site-by-site on the same row space. By default, write the");
  yame_usage_text("differential-methylation set as one format-6 track (set + universe). With -S,");
  yame_usage_text("print agreement statistics instead (one TSV line per pairing).");
  yame_usage_sec("Inputs:");
  yame_usage_text("<MU1.cx>   Format 3 (M/U) or format 4 (beta). Sample 1 is the record named by -1,");
  yame_usage_cont("or the first record.");
  yame_usage_text("[MU2.cx]   Optional second file. Sample 2 is the record named by -2, or the first");
  yame_usage_cont("record; with -S and no -2, EVERY record of MU2.cx is scored against");
  yame_usage_cont("sample 1. If omitted, sample 2 is the SECOND record of MU1.cx.");
  yame_usage_sec("Output:");
  yame_usage_text("Default: one format-6 record of length N (same as the inputs).");
  yame_usage_text("Universe: site i is in-universe only if BOTH samples cover it (-c for format 3;");
  yame_usage_cont("a value that is not NA for format 4).");
  yame_usage_text("Set:      site i is set if it passes the direction rule (-H) and effect threshold (-d).");
  yame_usage_text("-S:       TSV with a header: sample1 sample2 n mae rmse acc pearson, over the");
  yame_usage_cont("universe. acc = share of sites where both betas fall on the same side of -t.");
  yame_usage_sec("Options:");
  yame_usage_opt("-1 <name>", "Sample to take from MU1.cx (needs MU1.cx.idx).");
  yame_usage_opt("-2 <name>", "Sample to take from MU2.cx (needs MU2.cx.idx).");
  yame_usage_opt("-S", "Print agreement statistics instead of writing a set.");
  yame_usage_opt("-t <beta>", "Threshold for the acc column (default: 0.5). -S only.");
  yame_usage_opt("-o <out.cx>", "Write output to file (default: stdout).");
  yame_usage_opt("-c <cov>", "Minimum coverage (M+U) in a format-3 sample to count a site (default: 1).");
  yame_usage_opt("-d <delta>", "Minimum absolute beta difference required to call a site differential (default: 0).");
  yame_usage_opt("-H <mode>", "Direction mode (default: 1):");
  yame_usage_cont("1  beta1 > beta2  (hypermethylated in sample 1)");
  yame_usage_cont("2  beta1 < beta2  (hypomethylated  in sample 1)");
  yame_usage_cont("3  beta1 != beta2 (any difference; with -d uses |beta1-beta2|>delta)");
  yame_usage_opt("-h", "Show this help message.");
  yame_usage_sec("Notes:");
  yame_usage_text("* If you omit MU2.cx, MU1.cx must contain at least two records.");
  yame_usage_text("* The set output is binary; it does not store the delta magnitude.");
  yame_usage_text("* To score imputation, first blank the sites the model was shown:");
  yame_usage_cont("yame mask truth.cg query.cg | yame pairwise -S - prediction.cg");
  fprintf(stderr, "\n");

  return 1;
}

/* Beta at site i, or 0 when the side does not cover it: M+U below min_cov for
 * format 3, NA (stored negative) for format 4. The one place both formats are
 * read, so a third format fails loudly here rather than as numbers elsewhere. */
static int site_beta(cdata_t *c, uint64_t i, int64_t min_cov, double *b) {
  if (c->fmt == '3') {
    uint64_t mu = f3_get_mu(c, i);
    if (MU2cov(mu) < (uint64_t) min_cov) return 0;
    *b = MU2beta(mu);
    return 1;
  }
  float v = ((float*) c->s)[i];
  if (v < 0) return 0;
  *b = v;
  return 1;
}

/* Take one record from an open file: the one named through the index, or the
 * next in the stream. `label` is what the summary line calls it and `buf` is
 * the caller's storage for an ordinal label. A missing record is fatal when
 * `required`; otherwise it returns with n == 0, which is how a stream ends. */
static cdata_t take_record(cfile_t *cf, char *fname, char *name,
                           snames_t *names, int ordinal, int required,
                           char *buf, char **label) {
  cdata_t c = {0};
  if (name) {
    char *fname_index = get_fname_index(fname);
    index_t *idx = loadIndex(fname_index);
    if (!idx)
      wzfatal("[pairwise] -1/-2 need an index; none at %s (run yame index).\n",
              fname_index);
    free(fname_index);
    int64_t off = getIndex(idx, name);
    if (off < 0) wzfatal("[pairwise] %s is not in the index of %s.\n", name, fname);
    cleanIndex(idx);
    /* not an assert(): -DNDEBUG would delete the seek */
    if (bgzf_seek(cf->fh, off, SEEK_SET) != 0)
      wzfatal("[pairwise] Cannot seek to %s in %s.\n", name, fname);
    if (!read_cdata2(cf, &c) || c.n == 0)
      wzfatal("[pairwise] No record at the indexed position of %s in %s.\n",
              name, fname);
    *label = name;
  } else {
    if (!read_cdata2(cf, &c) || c.n == 0) {
      if (required) wzfatal("[pairwise] %s holds no record %d.\n", fname, ordinal + 1);
      free_cdata(&c);
      c.n = 0;
      return c;
    }
    if (names->n > ordinal) *label = names->s[ordinal];
    else { snprintf(buf, 32, "%d", ordinal + 1); *label = buf; }
  }
  if (c.fmt != '3' && c.fmt != '4')
    wzfatal("[pairwise] %s: need format 3 (M/U) or 4 (beta), got '%c'.\n",
            fname, c.fmt);
  decompress_in_situ(&c);
  return c;
}

/* Running sums for one pairing. Betas are in [0,1] and n is at most tens of
 * millions, so the raw-moment Pearson is well within double precision. */
typedef struct {
  uint64_t n, agree;
  double sae, sse, sx, sy, sxx, syy, sxy;
} pstat_t;

static void pstat_print(const pstat_t *s, const char *l1, const char *l2) {
  double n = (double) s->n;
  double den = (n * s->sxx - s->sx * s->sx) * (n * s->syy - s->sy * s->sy);
  double r = (s->n >= 2 && den > 0) ? (n * s->sxy - s->sx * s->sy) / sqrt(den)
                                     : NAN;
  fprintf(stdout, "%s\t%s\t%"PRIu64"\t%.6f\t%.6f\t%.6f\t%.6f\n", l1, l2, s->n,
          s->n ? s->sae / n : NAN, s->n ? sqrt(s->sse / n) : NAN,
          s->n ? (double) s->agree / n : NAN, r);
}

/* One pass over the two records. In set mode it fills c_out; in summary mode it
 * fills st. Both need the same universe, so they share the walk. */
static void walk(cdata_t *c1, cdata_t *c2, int64_t min_cov, int direc,
                 double min_effect, double thresh, cdata_t *c_out, pstat_t *st) {
  for (uint64_t i = 0; i < c1->n; ++i) {
    double b1, b2;
    if (!site_beta(c1, i, min_cov, &b1) || !site_beta(c2, i, min_cov, &b2))
      continue;
    if (st) {
      double d = b1 - b2;
      st->n++;
      st->sae += fabs(d);
      st->sse += d * d;
      st->agree += ((b1 >= thresh) == (b2 >= thresh));
      st->sx += b1; st->sy += b2;
      st->sxx += b1 * b1; st->syy += b2 * b2; st->sxy += b1 * b2;
      continue;
    }
    int hit;
    if (direc == 1)      hit = (b1 > b2 && b1 - b2 > min_effect);
    else if (direc == 2) hit = (b1 < b2 && b2 - b1 > min_effect);
    else                 hit = (min_effect <= 0) ? (b1 != b2)
                             : (b1 - b2 > min_effect || b2 - b1 > min_effect);
    if (hit) FMT6_SET1(*c_out, i); else FMT6_SET0(*c_out, i);
  }
}

int main_pairwise(int argc, char *argv[]) {

  int c; int direc = 1; double min_effect = -1.0;
  int64_t min_coverage = 1; char *fname_out = NULL;
  char *name1 = NULL, *name2 = NULL; int summary = 0; double thresh = 0.5;
  while ((c = getopt(argc, argv, "1:2:St:o:c:d:H:h"))>=0) {
    switch (c) {
    case '1': name1 = optarg; break;
    case '2': name2 = optarg; break;
    case 'S': summary = 1; break;
    case 't': thresh = atof(optarg); break;
    case 'o': fname_out = strdup(optarg); break;
    case 'c': min_coverage = atoi(optarg); break;
    case 'd': min_effect = atof(optarg); break;
    case 'H': direc = atoi(optarg); break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }
  if (min_coverage < 1) min_coverage = 1;
  if (!summary && (direc < 1 || direc > 3))
    wzfatal("[pairwise] -H must be 1, 2 or 3 (got %d).\n", direc);

  if (optind + 1 > argc) {
    usage();
    wzfatal("Please supply input file.\n");
  }

  char *fname1 = argv[optind++];
  char *fname2 = optind < argc ? argv[optind++] : NULL;

  cfile_t cf1 = open_cfile(fname1);
  snames_t names1 = loadSampleNamesFromIndex(fname1);
  char buf1[32], *label1;
  cdata_t c1 = take_record(&cf1, fname1, name1, &names1, 0, 1, buf1, &label1);

  /* Side 2 is the second record of the same file, a named record, or -- when
   * summarizing file 2 with no name -- every record of file 2 in turn. */
  cfile_t cf2 = cf1; char *fname2_eff = fname1; snames_t names2 = names1;
  int k = 1;                    /* ordinal of the next record on side 2 */
  if (fname2) {
    cf2 = open_cfile(fname2);
    fname2_eff = fname2;
    names2 = loadSampleNamesFromIndex(fname2);
    k = 0;
  }
  int broadcast = summary && fname2 && !name2;
  if (broadcast && !name1) {
    /* Two multi-record files with neither side named is a guess about which
     * pairing was meant. Peek for a second record in file 1 -- a compressed
     * read, no inflate -- and refuse rather than silently take the first. */
    cdata_t peek = {0};
    if (read_cdata2(&cf1, &peek) && peek.n > 0)
      wzfatal("[pairwise] %s has more than one record; say which with -1 "
              "(or which record of %s with -2).\n", fname1, fname2);
    free_cdata(&peek);
  }

  for (int first = 1; ; first = 0, ++k) {
    char buf2[32], *label2;
    /* the first pairing must exist; in a broadcast, running out is the end */
    cdata_t c2 = take_record(&cf2, fname2_eff, name2, &names2, k,
                             first || !broadcast, buf2, &label2);
    if (c2.n == 0) break;
    if (c1.n != c2.n)
      wzfatal("[pairwise] %s and %s differ in length: %"PRIu64" vs %"PRIu64
              " rows.\n", label1, label2, c1.n, c2.n);

    if (summary) {
      /* header only once a pairing has validated, so a refused run leaves
       * nothing on stdout for a sweep to append */
      if (first) fputs("sample1\tsample2\tn\tmae\trmse\tacc\tpearson\n", stdout);
      pstat_t st = {0};
      walk(&c1, &c2, min_coverage, direc, min_effect, thresh, NULL, &st);
      pstat_print(&st, label1, label2);
    } else {
      cdata_t c_out = {.fmt = '6', .n = c1.n };
      c_out.s = calloc((c_out.n+3)/4, sizeof(uint8_t));
      walk(&c1, &c2, min_coverage, direc, min_effect, thresh, &c_out, NULL);
      cdata_compress(&c_out);
      BGZF *fp_out;
      if (fname_out) fp_out = bgzf_open2(fname_out, "w");
      else fp_out = bgzf_dopen(fileno(stdout), "w");
      if (fp_out == NULL) {
        fprintf(stderr, "Error opening file for writing: %s\n", fname_out);
        exit(1);
      }
      cdata_write1(fp_out, &c_out);
      bgzf_close(fp_out);
      free_cdata(&c_out);
    }
    free_cdata(&c2);
    if (!broadcast) break;
  }

  free_cdata(&c1);
  if (fname2) { bgzf_close(cf2.fh); cleanSampleNames2(names2); }
  bgzf_close(cf1.fh);
  cleanSampleNames2(names1);
  return 0;
}
