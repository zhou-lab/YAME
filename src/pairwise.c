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

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: yame pairwise [options] <MU1.cx> (<MU2.cx>)\n");
  fprintf(stderr, "Return a format 6 set that represent differential methylation between MU1 and MU2.\n");
  fprintf(stderr, "If MU2 is not given, use the top 2 samples in MU1.cx.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -o        output cx file name. if missing, output to stdout without index.\n");
  fprintf(stderr, "    -H        1: B1>B2 (default).\n");
  fprintf(stderr, "              2: B1<B2.\n");
  fprintf(stderr, "              3: B1!=B2, i.e., 1 and 2 combined.\n");
  fprintf(stderr, "    -c        minimum coverage (default: 1)\n");
  fprintf(stderr, "    -d        minimum delta meth level/effect size (default: 0)\n");
  fprintf(stderr, "    -h        This help\n");
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

  decompress2(&c1);
  decompress2(&c2);

  if (c1.n != c2.n) wzfatal("Two inputs have different dimensions: %"PRIu64" vs %"PRIu64"\n", c1.n, c2.n);

  cdata_t c_out = {.fmt = '6', .n = c1.n };
  c_out.s = calloc((c_out.n+3)/4, sizeof(uint8_t));
  for (uint64_t i=0; i<c1.n; ++i) {
    uint64_t mu1 = f3_get_mu(&c1, i);
    uint64_t mu2 = f3_get_mu(&c2, i);
    if (MU2cov(mu1) >= (uint64_t) min_coverage &&
        MU2cov(mu2) >= (uint64_t) min_coverage) {
      if (direc == 1) {
        if (MU2beta(mu1) > MU2beta(mu2) &&
            MU2beta(mu1) - MU2beta(mu2) > min_effect) FMT6_SET1(c_out, i);
        else FMT6_SET0(c_out, i);
      } else if (direc == 2) {
        if (MU2beta(mu1) < MU2beta(mu2) &&
            MU2beta(mu2) - MU2beta(mu1) > min_effect) FMT6_SET1(c_out, i);
        else FMT6_SET0(c_out, i);
      } else if (direc == 3) {
        if ((min_effect <= 0 && MU2beta(mu1) != MU2beta(mu2)) ||
            (min_effect > 0 &&
             (MU2beta(mu1) - MU2beta(mu2) > min_effect ||
              MU2beta(mu2) - MU2beta(mu1) > min_effect))) FMT6_SET1(c_out, i);
        else FMT6_SET0(c_out, i);
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
