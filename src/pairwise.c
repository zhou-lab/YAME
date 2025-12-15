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
 * - Two format-3 (M/U) cdata records with the same length N.
 * - If MU2.cx is not provided, the code reads the first two records from MU1.cx.
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
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "  yame pairwise [options] <MU1.cx> [MU2.cx] > out.cx\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Purpose:\n");
  fprintf(stderr, "  Compute a per-site differential-methylation set between two format-3 (M/U) samples,\n");
  fprintf(stderr, "  and output it as a single format-6 track (set + universe).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Inputs:\n");
  fprintf(stderr, "  <MU1.cx>   Format-3 input (M/U counts). The first record is used as sample 1.\n");
  fprintf(stderr, "  [MU2.cx]   Optional second format-3 input. If omitted, sample 2 is read as the\n");
  fprintf(stderr, "            SECOND record from MU1.cx (i.e., the top 2 samples in the same file).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Output:\n");
  fprintf(stderr, "  One format-6 record of length N (same as the inputs).\n");
  fprintf(stderr, "  Universe: site i is in-universe only if BOTH samples have coverage >= min_cov.\n");
  fprintf(stderr, "  Set:      site i is set if it passes the direction rule (-H) and effect threshold (-d).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "  -o <out.cx>  Write output to file (default: stdout).\n");
  fprintf(stderr, "  -c <cov>     Minimum coverage (M+U) in BOTH samples to include site in universe (default: 1).\n");
  fprintf(stderr, "  -d <delta>   Minimum absolute beta difference required to call a site differential (default: 0).\n");
  fprintf(stderr, "  -H <mode>    Direction mode (default: 1):\n");
  fprintf(stderr, "              1  beta1 > beta2  (hypermethylated in sample 1)\n");
  fprintf(stderr, "              2  beta1 < beta2  (hypomethylated  in sample 1)\n");
  fprintf(stderr, "              3  beta1 != beta2 (any difference; with -d uses |beta1-beta2|>delta)\n");
  fprintf(stderr, "  -h           Show this help message.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Notes:\n");
  fprintf(stderr, "  * If you omit MU2.cx, MU1.cx must contain at least two records.\n");
  fprintf(stderr, "  * The output is a binary set; it does not store the delta magnitude.\n");
  fprintf(stderr, "\n");
  return 1;
}

int main_pairwise(int argc, char *argv[]) {

  int c; int direc = 1; double min_effect = -1.0;
  int64_t min_coverage = 1; char *fname_out = NULL;
  while ((c = getopt(argc, argv, "o:c:d:H:h"))>=0) {
    switch (c) {
    case 'o': fname_out = strdup(optarg); break;
    case 'c': min_coverage = atoi(optarg); break;
    case 'd': min_effect = atof(optarg); break;
    case 'H': direc = atoi(optarg); break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }
  if (min_coverage < 1) min_coverage = 1;

  if (optind + 1 > argc) {
    usage(); 
    wzfatal("Please supply input file.\n");
  }

  cfile_t cf1 = open_cfile(argv[optind++]);
  cdata_t c1 = read_cdata1(&cf1);
  cdata_t c2 = {0};
  if (optind >= argc) {
    c2 = read_cdata1(&cf1);
  } else {
    cfile_t cf2 = open_cfile(argv[optind++]);
    c2 = read_cdata1(&cf2);
    bgzf_close(cf2.fh);
  }
  bgzf_close(cf1.fh);

  decompress_in_situ(&c1);
  decompress_in_situ(&c2);

  if (c1.n != c2.n) wzfatal("Two inputs have different dimensions: %"PRIu64" vs %"PRIu64"\n", c1.n, c2.n);

  cdata_t c_out = {.fmt = '6', .n = c1.n };
  c_out.s = calloc((c_out.n+3)/4, sizeof(uint8_t));
  for (uint64_t i=0; i<c1.n; ++i) {
    uint64_t mu1 = f3_get_mu(&c1, i);
    uint64_t mu2 = f3_get_mu(&c2, i);
    if (MU2cov(mu1) >= (uint64_t) min_coverage &&
        MU2cov(mu2) >= (uint64_t) min_coverage) {
      if (direc == 1) { // sample 1 is more methylated than sample 2
        if (MU2beta(mu1) > MU2beta(mu2) &&
            MU2beta(mu1) - MU2beta(mu2) > min_effect) {
          FMT6_SET1(c_out, i);
        } else {
          FMT6_SET0(c_out, i);
        }
      } else if (direc == 2) {  // sample 1 is less methylated than sample 2
        if (MU2beta(mu1) < MU2beta(mu2) &&
            MU2beta(mu2) - MU2beta(mu1) > min_effect) {
          FMT6_SET1(c_out, i);
        } else {
          FMT6_SET0(c_out, i);
        }
      } else if (direc == 3) {  // sample 1 is different from sample 2
        if ((min_effect <= 0 && MU2beta(mu1) != MU2beta(mu2)) ||
            (min_effect > 0 &&
             (MU2beta(mu1) - MU2beta(mu2) > min_effect ||
              MU2beta(mu2) - MU2beta(mu1) > min_effect))) {
          FMT6_SET1(c_out, i);
        } else {
          FMT6_SET0(c_out, i);
        }
      } else {
        fprintf(stderr, "[%s:%d] -H argument: %d unsupported.\n", __func__, __LINE__, direc);
        fflush(stderr);
        exit(1);
      }
    }
  }
  free_cdata(&c1);
  free_cdata(&c2);
  
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
  
  return 0;
}
