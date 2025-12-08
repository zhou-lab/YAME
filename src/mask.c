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
  fprintf(stderr, "Usage: yame mask [options] <in.cg> <mask.cx>\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -o        output cx file name. if missing, output to stdout without index.\n");
  fprintf(stderr, "    -c        contextualize to format 6 using '1's in mask.\n");
  fprintf(stderr, "              if format 3 is used as mask, then use M+U>0 (coverage).\n");
  fprintf(stderr, "    -v        reverse the mask (default is to mask '1's, if -v will mask '0's).\n");
  fprintf(stderr, "    -h        This help\n");
  fprintf(stderr, "\n");

  return 1;
}

void mask_fmt3(cdata_t *c, cdata_t c_mask, BGZF *fp_out) {
  for (uint64_t i=0; i<c->n; ++i) {
    if (FMT0_IN_SET(c_mask, i)) {
      f3_set_mu(c, i, 0, 0);
    }
  }
  cdata_compress(c);
  cdata_write1(fp_out, c);
}

void mask_fmt0(cdata_t *c, cdata_t c_mask, BGZF *fp_out) {
  for (uint64_t i=0; i<cdata_nbytes(c); ++i) {
    c->s[i] &= ~c_mask.s[i];
  }
  /* cdata_compress(&c); */
  cdata_write1(fp_out, c);
}

void fmt0ContextualizeFmt6(cdata_t *c, cdata_t c_mask, BGZF *fp_out) {
  cdata_t c6 = {.fmt = '6', .n = c->n};
  c6.s = calloc((c6.n+3)/4, sizeof(uint8_t));
  for (uint64_t i=0; i<c6.n; ++i) {
    if (FMT0_IN_SET(c_mask,i)) { // mask is used as universe, use -v to invert
      if (FMT0_IN_SET(*c, i)) FMT6_SET1(c6, i);
      else FMT6_SET0(c6, i);
    }
  }
  cdata_compress(&c6);
  cdata_write1(fp_out, &c6);
  free_cdata(&c6);
}

int main_mask(int argc, char *argv[]) {

  int c, reverse = 0, contextualize_to_fmt6 = 0;
  char *fname_out = NULL;
  while ((c = getopt(argc, argv, "o:cvh"))>=0) {
    switch (c) {
    case 'o': fname_out = strdup(optarg); break;
    case 'c': contextualize_to_fmt6 = 1; break;
    case 'v': reverse = 1; break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 2 > argc) {
    usage(); 
    wzfatal("Please supply input file.\n");
  }

  char *fname = argv[optind];
  char *fname_mask = argv[optind+1];

  cfile_t cf_mask = open_cfile(fname_mask);
  cdata_t c_mask = read_cdata1(&cf_mask);
  if (c_mask.fmt == '1') convertToFmt0(&c_mask);
  if (c_mask.fmt == '3') convertToFmt0(&c_mask);
  if (c_mask.fmt != '0') wzfatal("Mask format not supported.");
  if (reverse) {
    for (uint64_t i=0; i<cdata_nbytes(&c_mask); ++i) {
      c_mask.s[i] = ~(c_mask.s[i]);
    }
  }

  BGZF *fp_out;
  if (fname_out) fp_out = bgzf_open2(fname_out, "w");
  else fp_out = bgzf_dopen(fileno(stdout), "w");
  if (fp_out == NULL) {
    fprintf(stderr, "Error opening file for writing: %s\n", fname_out);
    exit(1);
  }

  cfile_t cf = open_cfile(fname);
  while (1) {
    cdata_t c = read_cdata1(&cf);
    if (c.n == 0) break;
    if (c.fmt == '1') convertToFmt0(&c);
    decompress2(&c);
    if (c.n != c_mask.n) {
      fprintf(stderr, "[%s:%d] mask (n=%"PRIu64") and query (N=%"PRIu64") are of different lengths.\n", __func__, __LINE__, c_mask.n, c.n);
      fflush(stderr);
      exit(1);
    }

    if (c.fmt == '3') {
      mask_fmt3(&c, c_mask, fp_out);
    } else if (c.fmt == '0') {
      if (contextualize_to_fmt6) {
        fmt0ContextualizeFmt6(&c, c_mask, fp_out);
      } else {
        mask_fmt0(&c, c_mask, fp_out);
      }
    } else {
      fprintf(stderr, "[%s:%d] Only format %d files are supported.\n", __func__, __LINE__, c.fmt);
      fflush(stderr);
      exit(1);
    }
    free_cdata(&c);
  }

  if (fname_out) free(fname_out);
  bgzf_close(fp_out);
  free_cdata(&c_mask);
  bgzf_close(cf.fh);
  bgzf_close(cf_mask.fh);
  return 0;
}
