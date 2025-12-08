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

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: yame binarize [options] <mu.cg>\n");
  fprintf(stderr, "If Beta>0, then 1, else 0. M+U>0 is used as universe.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -t [T]    1 if Beta>T else 0. (default: 0.5).\n");
  fprintf(stderr, "    -m [T]    1 if M>=T else 0. (default: 0, if >0 will override -t.).\n");
  fprintf(stderr, "    -c        Minimum depth (M+U) (default: 1).\n");
  fprintf(stderr, "    -o        output cx file name (format 6). if missing, output to stdout.\n");
  fprintf(stderr, "    -h        This help\n");
  fprintf(stderr, "\n");

  return 1;
}

int main_binarize(int argc, char *argv[]) {

  int c; double T = 0.5; uint64_t min_cov = 1;
  uint64_t Mmin = 0;
  char *fname_out = NULL;
  while ((c = getopt(argc, argv, "o:t:m:c:h"))>=0) {
    switch (c) {
    case 'o': fname_out = strdup(optarg); break;
    case 't': T = atof(optarg); break;
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
    decompress2(&c);
    if (c.fmt != '3') {
      fprintf(stderr, "[%s:%d] Only format %d files are supported.\n", __func__, __LINE__, c.fmt);
      fflush(stderr);
      exit(1);
    }

    cdata_t c6 = {.fmt = '6', .n = c.n};
    c6.s = calloc((c6.n+3)/4, sizeof(uint8_t));
    for (uint64_t i=0; i<c6.n; ++i) {
      uint64_t mu = f3_get_mu(&c, i);
      if (MU2cov(mu) >= min_cov) {
        if (Mmin>0) {           /* binarize by just M */
          if ((mu>>32) >= Mmin) FMT6_SET1(c6, i);
          else FMT6_SET0(c6, i);
        } else {                /* binarize by Beta */
          if (MU2beta(mu)>T) FMT6_SET1(c6, i);
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
      if (!read_cdata2(&cf2, &c_tmp)) {
        fprintf(stderr, "[Error] Data is shorter than the sample name list.\n");
        fflush(stderr);
        exit(1);
      }
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
