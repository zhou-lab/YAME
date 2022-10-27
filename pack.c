#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdio.h>
#include "kycg.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: kycg pack [options] <in.bed> <out.cg>\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -f        format 0: 1 byte for 8 binary cpgs\n");
  fprintf(stderr, "                     1: value (1byte) + runlen (2bytes). Input is one ASCII\n");
  fprintf(stderr, "                     a: format 0 or format 1, whichever is smaller (default)\n");
  fprintf(stderr, "                     3: MU RLE + ladder byte \n");
  fprintf(stderr, "                     4: fraction / NA-RLE (32bytes)\n");
  fprintf(stderr, "                     5: 2-bit + NA-RLE. Input has only 0,1,2.\n");
  fprintf(stderr, "    -v        verbose\n");
  fprintf(stderr, "    -h        This help\n");
  fprintf(stderr, "\n");

  return 1;
}

/* static void process_2quaternaryVec(char *fname, char *fname_out) { */
/* } */

/* static void process_3um(char *fname, char *fname_out) { */
/* } */

/* The header design, 17 bytes
   uint64_t x 1: signature, used for validation
   uint8_t x 1: format (0=vec; 1=rle)
   uint64_t x 1: length (n_cgs or n_bytes for rle)
*/

int main_pack(int argc, char *argv[]) {

  int c; int verbose=0; char fmt='a';
  while ((c = getopt(argc, argv, "f:vh"))>=0) {
    switch (c) {
    case 'f': fmt = optarg[0]; break;
    case 'v': verbose = 1; break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 2 > argc) { 
    usage(); 
    wzfatal("Please supply input and output file.\n"); 
  }

  cgdata_t *cg;
  switch (fmt) {
  case 'a': {
    cg = fmt0_read_uncompressed(argv[optind], verbose);
    fmta_tryBinary2byteRLE_ifsmaller(cg);
    cgdata_write(argv[optind+1], cg, verbose);
    free_cgdata(cg);
    break;
  }
  case '0': {
    cg = fmt0_read_uncompressed(argv[optind], verbose);
    cgdata_write(argv[optind+1], cg, verbose);
    free_cgdata(cg);
    break;
  }
  case '1': {
    cg = fmt1_read_uncompressed(argv[optind], verbose);
    fmt1_compress(cg);
    cgdata_write(argv[optind+1], cg, verbose);
    free_cgdata(cg);
    break;
  }
    /* case '2': process_2quaternaryVec(argv[optind], argv[optind+1]); break; */
  case '3': {
    cg = fmt3_read_uncompressed(argv[optind], verbose);
    fmt3_compress(cg);
    cgdata_write(argv[optind+1], cg, verbose);
    free_cgdata(cg);
    break;
  }
  case '4': {
    cg = fmt4_read_uncompressed(argv[optind], verbose);
    fmt4_compress(cg);
    cgdata_write(argv[optind+1], cg, verbose);
    free_cgdata(cg);
    break;
  }
  case '5': {
    cg = fmt5_read_uncompressed(argv[optind], verbose);
    fmt5_compress(cg);
    cgdata_write(argv[optind+1], cg, verbose);
    free_cgdata(cg);
    break;
  }
  default: usage(); wzfatal("Unrecognized format: %c.\n", fmt);
  }
  return 0;
}

