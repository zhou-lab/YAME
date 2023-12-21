#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdio.h>
#include "cfile.h"

static int usage() {
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: yame pack [options] <in.txt> <out.cx>\n");
    fprintf(stderr, "The input text file must match the dimension and order of the reference CpG bed file.\n\n");

    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -f [char] Format specification (choose one character or number):\n");
    fprintf(stderr, "              (b) Binary data. Format default to 0 or 1 depending on size:\n");
    fprintf(stderr, "                  0 - 1 byte for 8 binary CpGs\n");
    fprintf(stderr, "                  1 - Value (1 byte) + Run-Length Encoding (RLE) (2 bytes)\n");
    fprintf(stderr, "              (s) State data. Format default to 2:\n");
    fprintf(stderr, "                  2 - State text + Index RLE (Best for chromatin states).\n");
    fprintf(stderr, "                      Use format 0 for sequence context.\n");
    fprintf(stderr, "              (m) Sequencing MU data. Format default to 3:\n");
    fprintf(stderr, "                  3 - MU RLE + Ladder byte (Input: 2-column text, M and U).\n");
    fprintf(stderr, "              (d) Differential meth data. Format default to 6:\n");
    fprintf(stderr, "                  5 - 2-bits + NA-RLE (Input: only 0, 1, 2 values).\n");
    fprintf(stderr, "                  6 - 2-bits boolean for S (set) and U (universe).\n");
    fprintf(stderr, "                      (Input: 2-column text, S and U).\n");
    fprintf(stderr, "              (n) Fraction data. Format default to 4:\n");
    fprintf(stderr, "                  4 - Fraction / NA-RLE (32 bytes).\n");
    fprintf(stderr, "              (r) Reference coordinates. Format default to 7:\n");
    fprintf(stderr, "                  7 - Compressed BED format for CGs.\n");

    fprintf(stderr, "    -u [int]  Number of bytes per unit data when inflated (1-8).\n");
    fprintf(stderr, "              Lower values are more memory efficient but may be lossier.\n");
    fprintf(stderr, "              0 - Inferred from data.\n");
    fprintf(stderr, "    -v        Verbose mode\n");
    fprintf(stderr, "    -h        Display this help message\n\n");

    return 1;
}

cdata_t *fmt0_read_raw(char *fname, int verbose);
cdata_t *fmt1_read_raw(char *fname, int verbose);
cdata_t *fmt2_read_raw(char *fname, int verbose);
cdata_t *fmt3_read_raw(char *fname, uint8_t unit, int verbose);
cdata_t *fmt4_read_raw(char *fname, int verbose);
cdata_t *fmt5_read_raw(char *fname, int verbose);
cdata_t *fmt6_read_raw(char *fname, int verbose);
cdata_t *fmt7_read_raw(char *fname, int verbose);

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

  cdata_t *c;
  switch (fmt) {
  case 'b': {
    c = fmt0_read_raw(argv[optind], verbose);
    fmta_tryBinary2byteRLE_ifsmaller(c);
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
  case '5': {
    c = fmt5_read_raw(argv[optind], verbose);
    break;
  }
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

