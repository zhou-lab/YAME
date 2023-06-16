#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdio.h>
#include "kycg.h"

void cgdata_write(char *fname_out, cgdata_t *cg, const char *mode, int verbose) {

  if (!cg->compressed) recompress(cg);

  /* FILE *out; */
  /* if (fname_out) out = fopen(fname_out, mode); */
  /* else out = stdout; */
  /* uint64_t sig = CGSIG; fwrite(&sig, sizeof(uint64_t), 1, out); */
  /* fwrite(&cg->fmt, sizeof(uint8_t), 1, out); */
  /* fwrite(&cg->n, sizeof(uint64_t), 1, out); */
  /* fwrite(cg->s, 1, cgdata_nbytes(cg), out); */
  /* fclose(out); */

  BGZF* fp;
  if (fname_out) fp = bgzf_open2(fname_out, mode);
  else fp = bgzf_dopen(fileno(stdout), mode);
    
  if (fp == NULL) {
    fprintf(stderr, "Error opening file for writing: %s\n", fname_out);
    return;
  }

  // Write the signature
  uint64_t sig = CGSIG;
  if (bgzf_write(fp, &sig, sizeof(uint64_t)) < 0) {
    fprintf(stderr, "Error writing signature to file\n");
    bgzf_close(fp);
    return;
  }

  // Write the format
  if (bgzf_write(fp, &(cg->fmt), sizeof(uint8_t)) < 0) {
    fprintf(stderr, "Error writing format to file\n");
    bgzf_close(fp);
    return;
  }

  // Write the count
  if (bgzf_write(fp, &(cg->n), sizeof(uint64_t)) < 0) {
    fprintf(stderr, "Error writing count to file\n");
    bgzf_close(fp);
    return;
  }

  // Write the data
  if (bgzf_write(fp, cg->s, cgdata_nbytes(cg)) < 0) {
    fprintf(stderr, "Error writing data to file\n");
    bgzf_close(fp);
    return;
  }

  // Close the file
  bgzf_close(fp);

  if (verbose) {
    fprintf(stderr, "[%s:%d] Stored as Format %c\n", __func__, __LINE__, cg->fmt);
    fflush(stderr);
  }
}

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

  if (optind + 1 > argc) { 
    usage(); 
    wzfatal("Please supply input file.\n"); 
  }


  char *fname_out = NULL;
  if (argc >= optind + 2)
    fname_out = strdup(argv[optind+1]);
  
  cgdata_t *cg;
  switch (fmt) {
  case 'a': {
    cg = fmt0_read_uncompressed(argv[optind], verbose);
    fmta_tryBinary2byteRLE_ifsmaller(cg);
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
    /* case '2': process_2quaternaryVec(argv[optind], fname_out); break; */
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
  free_cgdata(cg);

  if (fname_out) free(fname_out);
  return 0;
}

