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
  fprintf(stderr, "    -f        format (0,1,2,3,4,a)\n");
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

void process_0binaryVec(char *fname, char *fname_out, int verbose);
void process_a(char *fname, char *fname_out, int verbose);
void process_1byteRLE(char *fname, char *fname_out, int verbose);
void process_4fractionNA(char *fname, char *fname_out, int verbose);
void process_5byteRLE(char *fname, char *fname_out, int verbose);


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

  switch (fmt) {
  case 'a': {
    cgdata_t *cg = fmt0_read_uncompressed(argv[optind], verbose);
    fmta_tryBinary2byteRLE_ifsmaller(cg);
    cgdata_write(argv[optind+1], cg, verbose);
    free_cgdata(cg);
    break;
  }
  case '0': {
    cgdata_t *cg = fmt0_read_uncompressed(argv[optind], verbose);
    cgdata_write(argv[optind+1], cg, verbose);
    free_cgdata(cg);
  }
  case '1': {
    cgdata_t *cg = fmt1_read_uncompressed(argv[optind], verbose);
    fmt1_compress(cg);
    cgdata_write(argv[optind+1], cg, verbose);
    free_cgdata(cg);
    break;
  }
  /* case '2': process_2quaternaryVec(argv[optind], argv[optind+1]); break; */
  case '3': {
    cgdata_t *cg = fmt3_read_uncompressed(argv[optind], verbose);
    fmt3_compress(cg);
    cgdata_write(argv[optind+1], cg, verbose);
    free_cgdata(cg);
    break;
  }
  case '4': {
    cgdata_t *cg = fmt4_read_uncompressed(argv[optind], verbose);
    fmt4_compress(cg);
    cgdata_write(argv[optind+1], cg, verbose);
    free_cgdata(cg);
    break;
  }
  case '5': {
    cgdata_t *cg = fmt5_read_uncompressed(argv[optind], verbose);
    fmt5_compress(cg);
    cgdata_write(argv[optind+1], cg, verbose);
    free_cgdata(cg);
    break;
  }
  default: usage(); wzfatal("Unrecognized format: %c.\n", fmt);
  }
  return 0;
}

static void unpack_5byteRLE(uint8_t *s, uint64_t n) {
  uint64_t i=0, j;
  for (i=0; i<n; ++i) {
    if (s[i] & (1<<7)) {
      int offset = 6;
      for (offset = 6; offset >= 0; offset -= 2) {
        if ((s[i]>>offset) & 0x2) {
          fputc(((s[i]>>offset) & 0x1)+'0', stdout);
          /* fputc('\n', stdout); */
          fprintf(stdout, "\t%u\n", s[i]);
        } else {
          break;
        }
      }
    } else {
      for (j=0; j < s[i]; ++j) {
        fputc('2', stdout);
        fputc('\n', stdout);
      }
    }
  }
}

int main_unpack(int argc, char *argv[]) {

  int c; int verbose = 0;
  while ((c = getopt(argc, argv, "vh"))>=0) {
    switch (c) {
    case 'v': verbose = 1; break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 1 > argc) { 
    usage(); 
    wzfatal("Please supply input file.\n"); 
  }

  FILE *fh = fopen(argv[optind],"rb");
  uint64_t sig;
  fread(&sig, sizeof(uint64_t), 1, fh);
  if (sig != CGSIG) {
    wzfatal("Unmatched signature. File corrupted.\n");
  }
  
  char fmt;
  fread(&fmt, sizeof(char), 1, fh);
  uint64_t n;
  fread(&n, sizeof(uint64_t), 1, fh);
  uint8_t *s = malloc(n);
  fread(s, 1, n, fh);
  
  switch (fmt) {
  /* case '0': unpack_0binaryVec(s, n, verbose); break; */
  /* case '1': unpack_1byteRLE(s, n, verbose); break; */
  /* case '4': unpack_4fractionNA(s, n, verbose); break; */
  case '5': unpack_5byteRLE(s, n); break;
  default: usage(); wzfatal("Unrecognized format: %c.\n", fmt);
  }
  return 0;
  
}
