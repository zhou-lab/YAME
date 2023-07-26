#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdio.h>
#include "cgfile.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: kycg pack [options] <in.bed> <out.cg>\n");
  fprintf(stderr, "Please only supply one of -b, -s, -m, -n, -f.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -b        binary data (.cb). Possible formats:\n");
  fprintf(stderr, "              (0) 1 byte for 8 binary cpgs\n"); 
  fprintf(stderr, "              (1) value (1 byte) + runlen (2 bytes). Input is one ASCII\n");
  fprintf(stderr, "              (5) 2-bit + NA-RLE. Input has only 0, 1, 2.\n");
  fprintf(stderr, "    -s        state data (.cs):\n");
  fprintf(stderr, "              (2) state text + index RLE.\n");
  fprintf(stderr, "    -m        sequencing MU data (.cm):\n");
  fprintf(stderr, "              (3) MU RLE + ladder byte.\n");
  fprintf(stderr, "    -n        fraction data (.cn):\n");
  fprintf(stderr, "              (4) fraction / NA-RLE (32 bytes)\n");
  fprintf(stderr, "    -f [int]  format. See options above. -b, -s, -m, -n provide the default based on the data size\n");
  fprintf(stderr, "              This option specifies the exact format.\n");
  fprintf(stderr, "    -v        verbose\n");
  fprintf(stderr, "    -h        This help\n");
  fprintf(stderr, "\n");
  return 1;
}

cgdata_t *fmt0_read_uncompressed(char *fname, int verbose);
cgdata_t *fmt1_read_uncompressed(char *fname, int verbose);
cgdata_t *fmt2_read_uncompressed(char *fname, int verbose);
cgdata_t *fmt3_read_uncompressed(char *fname, int verbose);
cgdata_t *fmt4_read_uncompressed(char *fname, int verbose);
cgdata_t *fmt5_read_uncompressed(char *fname, int verbose);
cgdata_t *fmt6_read_uncompressed(char *fname, int verbose);

int main_pack(int argc, char *argv[]) {

  int c; int verbose=0; char fmt='a';
  while ((c = getopt(argc, argv, "bsmnf:vh"))>=0) {
    switch (c) {
    case 'b': fmt = 'b'; break;
    case 's': fmt = 's'; break;
    case 'm': fmt = 'm'; break;
    case 'n': fmt = 'n'; break;
    case 'f': fmt = optarg[0]; break;
    case 'v': verbose = 1; break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 1 > argc) { 
    usage(); 
    wzfatal("Please supply input file.\n"); 
  }

  char *fname_out = NULL;
  if (argc >= optind + 2)
    fname_out = strdup(argv[optind+1]);

  cgdata_t *cg;
  switch (fmt) {
  case 'b': {
    cg = fmt0_read_uncompressed(argv[optind], verbose);
    fmta_tryBinary2byteRLE_ifsmaller(cg);
    break;
  }
  case 's': {
    cg = fmt2_read_uncompressed(argv[optind], verbose);
    break;
  }
  case 'm': {
    cg = fmt3_read_uncompressed(argv[optind], verbose);
    break;
  }
  case 'n': {
    cg = fmt4_read_uncompressed(argv[optind], verbose);
    break;
  }
  case '0': {
    cg = fmt0_read_uncompressed(argv[optind], verbose);
    break;
  }
  case '1': {
    cg = fmt1_read_uncompressed(argv[optind], verbose);
    fmt1_compress(cg);
    break;
  }
  case '2': {
    cg = fmt2_read_uncompressed(argv[optind], verbose);
    fmt2_compress(cg);
    break;
  }
  case '3': {
    cg = fmt3_read_uncompressed(argv[optind], verbose);
    fmt3_compress(cg);
    break;
  }
  case '4': {
    cg = fmt4_read_uncompressed(argv[optind], verbose);
    fmt4_compress(cg);
    break;
  }
  case '5': {
    cg = fmt5_read_uncompressed(argv[optind], verbose);
    fmt5_compress(cg);
    break;
  }
  case '6': {
    cg = fmt6_read_uncompressed(argv[optind], verbose);
    fmt6_compress(cg);
    break;
  }
  default: usage(); wzfatal("Unrecognized format: %c.\n", fmt);
  }
  cgdata_write(fname_out, cg, "w", verbose);
  free_cgdata(cg); free(cg);

  if (fname_out) free(fname_out);
  return 0;
}

