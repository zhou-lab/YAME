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

#include <sys/stat.h>
#include <sys/types.h>
#include "cfile.h"

/**
 * yame binarize
 * =============
 *
 * Goal
 * ----
 * Convert format-3 per-site methylation counts (M,U) into format-6 packed
 * set/universe representation suitable for fast overlap/enrichment operations.
 *
 * Format-6 semantics
 * ------------------
 * Each genomic position i is encoded by two bits (packed 4 sites per byte):
 *   - universe bit (U_i): whether this position is considered "measured/valid"
 *   - set bit      (S_i): whether this position is in the "positive" set
 *
 * In YAME macros (see cdata.h):
 *   FMT6_IN_UNI(c,i) tests universe membership
 *   FMT6_IN_SET(c,i) tests set membership
 * and FMT6_SET0/SET1/SET_NA write the corresponding states. :contentReference[oaicite:1]{index=1}
 *
 * Universe definition
 * -------------------
 * For each site, the input MU is fetched from format-3 (f3_get_mu), and
 * coverage is MU2cov(mu) = M+U.
 *
 *   If (M+U) >= min_cov:
 *       universe bit is set to 1 (site is included in universe)
 *   Else:
 *       universe bit remains 0 (site is NA/outside-universe)
 *
 * Set definition (two modes)
 * --------------------------
 * After passing universe filter, the set bit is assigned by one of:
 *
 *  (1) Beta-threshold mode (default):
 *      beta = MU2beta(mu) = M/(M+U)
 *      set=1 if beta > T, else set=0
 *
 *  (2) M-count mode (-m Mmin, overrides -t):
 *      set=1 if M >= Mmin, else set=0
 *
 * Multi-sample behavior
 * ---------------------
 * The input may contain multiple format-3 records (samples). For each record:
 *   - decompress in-place
 *   - build a new format-6 record of the same length
 *   - compress and write it
 *
 * Index propagation
 * -----------------
 * If an index exists for the input and -o is provided, the tool generates
 * a matching index for the output by replaying sample offsets. If output is
 * stdout, no index is written. :contentReference[oaicite:2]{index=2}
 */

static int usage(void) {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "  yame binarize [options] <mu.cx>\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Purpose:\n");
  fprintf(stderr, "  Convert per-site M/U counts (format 3) into a packed binary-with-universe\n");
  fprintf(stderr, "  track (format 6).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Input / Output:\n");
  fprintf(stderr, "  Input : format 3 (.cx) with per-site (M,U) stored as uint64.\n");
  fprintf(stderr, "  Output: format 6 (.cx), where each site stores two bits:\n");
  fprintf(stderr, "          - universe bit: 1 if depth>=min_cov, else 0 (NA/outside-universe)\n");
  fprintf(stderr, "          - set bit:      1 if methylated by rule, else 0\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Binarization rules:\n");
  fprintf(stderr, "  Default: set=1 if beta > T (beta = M/(M+U)), set=0 otherwise.\n");
  fprintf(stderr, "  If -m is provided (>0): set=1 if M >= Mmin, else 0 (overrides -t).\n");
  fprintf(stderr, "  Universe is always defined by coverage: (M+U) >= min_cov.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "  -t <Tmin>   Beta threshold (default: 0.5).\n");
  fprintf(stderr, "  -m <Mmin>   M-count threshold (default: 0; if >0 overrides -t).\n");
  fprintf(stderr, "  -c <cov>    Minimum coverage (M+U) to include a site in universe (default: 1).\n");
  fprintf(stderr, "  -o <out.cx> Write output to file (default: stdout).\n");
  fprintf(stderr, "  -h          Show this help message.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Notes:\n");
  fprintf(stderr, "  * Sites with depth < min_cov remain NA in format 6 (universe bit = 0).\n");
  fprintf(stderr, "  * If the input has a sample index and -o is used, an output index is written.\n");
  fprintf(stderr, "\n");
  return 1;
}


int main_binarize(int argc, char *argv[]) {

  int c; double Tmin = 0.5; uint64_t min_cov = 1;
  uint64_t Mmin = 0;
  char *fname_out = NULL;
  while ((c = getopt(argc, argv, "o:t:m:c:h"))>=0) {
    switch (c) {
    case 'o': fname_out = strdup(optarg); break;
    case 't': Tmin = atof(optarg); break;
    case 'm': Mmin = atoi(optarg); break;
    case 'c': min_cov = atoi(optarg); break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 1 > argc) {
    usage(); 
    wzfatal("Please supply input file.\n");
  }

  char *fname = argv[optind];

  BGZF *fp_out;
  if (fname_out) fp_out = bgzf_open2(fname_out, "w");
  else fp_out = bgzf_dopen(fileno(stdout), "w");
  if (fp_out == NULL) {
    fprintf(stderr, "Error opening file for writing: %s\n", fname_out);
    exit(1);
  }

  cfile_t cf = open_cfile(fname);
  char *fname_index = get_fname_index(fname);
  index_t *idx = loadIndex(fname_index);
  free(fname_index);

  while (1) {
    cdata_t c = read_cdata1(&cf);
    if (c.n == 0) break;
    decompress_in_situ(&c);
    if (c.fmt != '3') wzfatal("[%s:%d] Only format 3 files are supported (given %c).\n", __func__, __LINE__, c.fmt);

    cdata_t c6 = {.fmt = '6', .n = c.n};
    c6.s = calloc((c6.n+3)/4, sizeof(uint8_t));
    for (uint64_t i=0; i<c6.n; ++i) {
      uint64_t mu = f3_get_mu(&c, i);
      if (MU2cov(mu) >= min_cov) {
        if (Mmin>0) {           /* binarize by just M */
          if ((mu>>32) >= Mmin) FMT6_SET1(c6, i);
          else FMT6_SET0(c6, i);
        } else {                /* binarize by Beta */
          if (MU2beta(mu) >= Tmin) FMT6_SET1(c6, i);
          else FMT6_SET0(c6, i);
        }
      }
    }
    cdata_compress(&c6);
    cdata_write1(fp_out, &c6);
    free_cdata(&c6);
    free_cdata(&c);
  }

  if (idx && fname_out) {              // output index
    int npairs = 0;
    index_pair_t *pairs = index_pairs(idx, &npairs);
    
    cfile_t cf2 = open_cfile(fname_out);
    index_t *idx2 = kh_init(index);
    int64_t addr = bgzf_tell(cf2.fh);
    cdata_t c_tmp = {0};
    for (int i=0; i< npairs; ++i) {
      if (!read_cdata2(&cf2, &c_tmp))
        wzfatal("[Error] Data is shorter than the sample name list.\n");
      
      insert_index(idx2, pairs[i].key, addr);
      addr = bgzf_tell(cf2.fh);
    }
    free_cdata(&c_tmp);

    char *fname_index2 = get_fname_index(fname_out);
    FILE *out = fopen(fname_index2, "w");
    writeIndex(out, idx2);
    fclose(out);
    free(fname_index2);
    bgzf_close(cf2.fh);
    freeIndex(idx2);
    free(pairs);
    cleanIndex(idx);
  }

  if (fname_out) free(fname_out);
  bgzf_close(fp_out);
  bgzf_close(cf.fh);

  return 0;
}
