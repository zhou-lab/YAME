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
#include "yame_ui.h"
#include <string.h>
#include <zlib.h>
#include <stdio.h>
#include "cfile.h"

static int usage() {
  yame_usage_head("yame pack [options] <in.txt> <out.cx>");
  yame_usage_text("Pack tab-delimited text into a compressed cx file.");
  yame_usage_text("The input file must have one row per CpG and match the");
  yame_usage_text("dimension and order of the reference CpG BED file.");
  yame_usage_sec("Options:");
  yame_usage_opt("-f [char]", "Format specification (one of b,c,s,m,d,n,r):");
  yame_usage_cont("(b) Binary data (format 0).");
  yame_usage_cont("    Each entry is 0 or 1.");
  yame_usage_cont("    Example (single-sample, one column):");
  yame_usage_cont("        0");
  yame_usage_cont("        1");
  yame_usage_cont("        1");
  yame_usage_cont("(c) Character / small integer data (format 1).");
  yame_usage_cont("    One byte per entry, typically 0–255.");
  yame_usage_cont("    Example:");
  yame_usage_cont("        0");
  yame_usage_cont("        5");
  yame_usage_cont("        9");
  yame_usage_cont("(s) State data (format 2).");
  yame_usage_cont("    Categorical strings compressed via an index + RLE.");
  yame_usage_cont("    Best for chromatin states or other labels.");
  yame_usage_cont("    Example:");
  yame_usage_cont("        quies");
  yame_usage_cont("        quies");
  yame_usage_cont("        enhA");
  yame_usage_cont("(m) Sequencing MU data (format 3).");
  yame_usage_cont("    Input is 2-column text: M and U counts per CpG.");
  yame_usage_cont("    M=U=0 is treated as missing.");
  yame_usage_cont("    Example (M U):");
  yame_usage_cont("        10	5");
  yame_usage_cont("        20	0");
  yame_usage_cont("        13	17");
  yame_usage_cont("(d) Differential / mask data (format 6).");
  yame_usage_cont("    2-bit boolean for S (set) and U (universe).");
  yame_usage_cont("    Input is 2-column text: S and U, each 0 or 1.");
  yame_usage_cont("    Example (S U):");
  yame_usage_cont("        1	1");
  yame_usage_cont("        0	1");
  yame_usage_cont("        0	0");
  yame_usage_cont("(n) Fraction / beta data (format 4).");
  yame_usage_cont("    Floating-point fraction in [0,1] or NA.");
  yame_usage_cont("    Example:");
  yame_usage_cont("        0.250");
  yame_usage_cont("        NA");
  yame_usage_cont("        1.000");
  yame_usage_cont("(r) Reference coordinates (format 7).");
  yame_usage_cont("    Compressed BED records for CpG coordinates.");
  yame_usage_cont("    Input is 4-column BED: chrom, start, end, name.");
  yame_usage_cont("    Example:");
  yame_usage_cont("        chr1	100	101	CpG_1");
  yame_usage_cont("        chr1	200	201	CpG_2");
  yame_usage_cont("        chr1	300	301	CpG_3");
  yame_usage_cont("The examples above show single-sample input.");
  yame_usage_cont("Multi-sample input can be provided as additional");
  yame_usage_cont("columns per row, following the same conventions.");
  yame_usage_opt("-u [int]", "Number of bytes per unit when inflated (1-8).");
  yame_usage_cont("Lower values are more memory efficient but may be lossier.");
  yame_usage_cont("0 - infer from data.");
  yame_usage_opt("-v", "Verbose mode.");
  yame_usage_opt("-h", "Display this help message.");
  fprintf(stderr, "\n");

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

