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

#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdio.h>
#include "cfile.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: yame pack [options] <in.txt> <out.cx>\n");
  fprintf(stderr, "Pack tab-delimited text into a compressed cx file.\n");
  fprintf(stderr, "The input file must have one row per CpG and match the\n");
  fprintf(stderr, "dimension and order of the reference CpG BED file.\n\n");

  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -f [char] Format specification (one of b,c,s,m,d,n,r):\n");
  fprintf(stderr, "              (b) Binary data (format 0).\n");
  fprintf(stderr, "                  Each entry is 0 or 1.\n");
  fprintf(stderr, "                  Example (single-sample, one column):\n");
  fprintf(stderr, "                      0\n");
  fprintf(stderr, "                      1\n");
  fprintf(stderr, "                      1\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "              (c) Character / small integer data (format 1).\n");
  fprintf(stderr, "                  One byte per entry, typically 0–255.\n");
  fprintf(stderr, "                  Example:\n");
  fprintf(stderr, "                      0\n");
  fprintf(stderr, "                      5\n");
  fprintf(stderr, "                      9\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "              (s) State data (format 2).\n");
  fprintf(stderr, "                  Categorical strings compressed via an index + RLE.\n");
  fprintf(stderr, "                  Best for chromatin states or other labels.\n");
  fprintf(stderr, "                  Example:\n");
  fprintf(stderr, "                      quies\n");
  fprintf(stderr, "                      quies\n");
  fprintf(stderr, "                      enhA\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "              (m) Sequencing MU data (format 3).\n");
  fprintf(stderr, "                  Input is 2-column text: M and U counts per CpG.\n");
  fprintf(stderr, "                  M=U=0 is treated as missing.\n");
  fprintf(stderr, "                  Example (M U):\n");
  fprintf(stderr, "                      10\t5\n");
  fprintf(stderr, "                      20\t0\n");
  fprintf(stderr, "                      13\t17\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "              (d) Differential / mask data (format 6).\n");
  fprintf(stderr, "                  2-bit boolean for S (set) and U (universe).\n");
  fprintf(stderr, "                  Input is 2-column text: S and U, each 0 or 1.\n");
  fprintf(stderr, "                  Example (S U):\n");
  fprintf(stderr, "                      1\t1\n");
  fprintf(stderr, "                      0\t1\n");
  fprintf(stderr, "                      0\t0\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "              (n) Fraction / beta data (format 4).\n");
  fprintf(stderr, "                  Floating-point fraction in [0,1] or NA.\n");
  fprintf(stderr, "                  Example:\n");
  fprintf(stderr, "                      0.250\n");
  fprintf(stderr, "                      NA\n");
  fprintf(stderr, "                      1.000\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "              (r) Reference coordinates (format 7).\n");
  fprintf(stderr, "                  Compressed BED records for CpG coordinates.\n");
  fprintf(stderr, "                  Input is 4-column BED: chrom, start, end, name.\n");
  fprintf(stderr, "                  Example:\n");
  fprintf(stderr, "                      chr1\t100\t101\tCpG_1\n");
  fprintf(stderr, "                      chr1\t200\t201\tCpG_2\n");
  fprintf(stderr, "                      chr1\t300\t301\tCpG_3\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "              The examples above show single-sample input.\n");
  fprintf(stderr, "              Multi-sample input can be provided as additional\n");
  fprintf(stderr, "              columns per row, following the same conventions.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "    -u [int]  Number of bytes per unit when inflated (1-8).\n");
  fprintf(stderr, "              Lower values are more memory efficient but may be lossier.\n");
  fprintf(stderr, "              0 - infer from data.\n");
  fprintf(stderr, "    -v        Verbose mode.\n");
  fprintf(stderr, "    -h        Display this help message.\n\n");

  return 1;
}

cdata_t *fmt0_read_raw(char *fname, int verbose);
cdata_t *fmt1_read_raw(char *fname, int verbose);
cdata_t *fmt2_read_raw(char *fname, int verbose);
cdata_t *fmt3_read_raw(char *fname, uint8_t unit, int verbose);
cdata_t *fmt4_read_raw(char *fname, int verbose);
/* cdata_t *fmt5_read_raw(char *fname, int verbose); */
cdata_t *fmt6_read_raw(char *fname, int verbose);
cdata_t *fmt7_read_raw(char *fname, int verbose);
/* void fmta_tryBinary2byteRLE_ifsmaller(cdata_t *c); */

int main_pack(int argc, char *argv[]) {

  int c0; int verbose=0; char fmt='a'; uint8_t unit = 8;
  while ((c0 = getopt(argc, argv, "f:u:vh"))>=0) {
    switch (c0) {
    case 'f': fmt = optarg[0]; break;
    case 'u': unit = atoi(optarg); break;
    case 'v': verbose = 1; break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c0);
    }
  }

  if (optind + 1 > argc) { 
    usage(); 
    wzfatal("Please supply input file.\n"); 
  }

  char *fname_out = NULL;
  if (argc >= optind + 2)
    fname_out = strdup(argv[optind+1]);

  cdata_t *c = NULL;
  switch (fmt) {
  case 'b': {
    c = fmt0_read_raw(argv[optind], verbose);
    /* fmta_tryBinary2byteRLE_ifsmaller(c); */
    break;
  }
  case 'c': {
    c = fmt1_read_raw(argv[optind], verbose);
    break;
  }
  case 'd': {
    c = fmt6_read_raw(argv[optind], verbose);
    break;
  }
  case 's': {
    c = fmt2_read_raw(argv[optind], verbose);
    break;
  }
  case 'm': {
    c = fmt3_read_raw(argv[optind], unit, verbose);
    break;
  }
  case 'n': {
    c = fmt4_read_raw(argv[optind], verbose);
    break;
  }
  case 'r': {
    c = fmt7_read_raw(argv[optind], verbose);
    break;
  }
  case '0': {
    c = fmt0_read_raw(argv[optind], verbose);
    break;
  }
  case '1': {
    c = fmt1_read_raw(argv[optind], verbose);
    break;
  }
  case '2': {
    c = fmt2_read_raw(argv[optind], verbose);
    break;
  }
  case '3': {
    c = fmt3_read_raw(argv[optind], unit, verbose);
    break;
  }
  case '4': {
    c = fmt4_read_raw(argv[optind], verbose);
    break;
  }
  /* case '5': { */
  /*   c = fmt5_read_raw(argv[optind], verbose); */
  /*   break; */
  /* } */
  case '6': {
    c = fmt6_read_raw(argv[optind], verbose);
    break;
  }
  case '7': {
    c = fmt7_read_raw(argv[optind], verbose);
    break;
  }
  default: usage(); wzfatal("Unrecognized format: %c.\n", fmt);
  }
  cdata_write(fname_out, c, "w", verbose);
  free_cdata(c); free(c);

  if (fname_out) free(fname_out);
  return 0;
}

