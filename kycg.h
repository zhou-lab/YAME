#ifndef _TBMATE_H
#define _TBMATE_H

#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <limits.h>
#include <inttypes.h>
#include <wordexp.h>
#include "wzmisc.h"
#include "wzmisc.h"
#include "wzbed.h"

#define PACKAGE_VERSION "0.1.20221020"
#define CGSIG 266563789635

typedef struct cgdata_t {
  uint8_t *s;
  uint64_t n; /* number of bytes, except for fmt 0, which is sub-byte you need the actual length*/
  int compressed;
  char fmt;
} cgdata_t;

DEFINE_VECTOR(cgdata_v, cgdata_t)

static inline uint64_t cgdata_nbytes(cgdata_t *cg) {
  uint64_t n = 0;
  switch(cg->fmt) {
  case '0': n = (cg->n>>3)+1; break;
  default: n = cg->n;
  }
  return n;
}

static inline void free_cgdata(cgdata_t *cg) {
  if(cg->s) free(cg->s);
  free(cg);
}

void fmta_tryBinary2byteRLE_ifsmaller(cgdata_t *cg);

cgdata_t *fmt0_read_uncompressed(char *fname, int verbose);
void fmt0_compress(cgdata_t *cg);
cgdata_t fmt0_decompress(cgdata_t *cg);

cgdata_t *fmt1_read_uncompressed(char *fname, int verbose);
void fmt1_compress(cgdata_t *cg);
cgdata_t fmt1_decompress(cgdata_t *cg);

cgdata_t *fmt3_read_uncompressed(char *fname, int verbose);
void fmt3_compress(cgdata_t *cg);
cgdata_t fmt3_decompress(cgdata_t *cg);

cgdata_t *fmt4_read_uncompressed(char *fname, int verbose);
void fmt4_compress(cgdata_t *cg);
cgdata_t fmt4_decompress(cgdata_t *cg);

cgdata_t *fmt5_read_uncompressed(char *fname, int verbose);
void fmt5_compress(cgdata_t *cg);
cgdata_t fmt5_decompress(cgdata_t *cg);

cgdata_t decompress(cgdata_t *cg);
void compress(cgdata_t *cg);

cgdata_t slice(cgdata_t *cg, uint64_t beg, uint64_t end);

static inline void cgdata_write(char *fname_out, cgdata_t *cg, int verbose) {

  FILE *out = fopen(fname_out, "wb");
  uint64_t sig = CGSIG; fwrite(&sig, sizeof(uint64_t), 1, out);
  fwrite(&cg->fmt, sizeof(uint8_t), 1, out);
  fwrite(&cg->n, sizeof(uint64_t), 1, out);
  fwrite(cg->s, 1, cgdata_nbytes(cg), out);
  fclose(out);

  if (verbose) {
    fprintf(stderr, "[%s:%d] Stored as Format %c\n", __func__, __LINE__, cg->fmt);
    fflush(stderr);
  }
}

/* static inline cgdata_t* read_cg(const char *fname) { */
/*   FILE *fh = fopen(fname,"rb"); */

/*   cgdata_t *cg = calloc(sizeof(cgdata_t),1); */

/*   uint64_t sig; */
/*   fread(&sig, sizeof(uint64_t), 1, fh); */
/*   if (sig != CGSIG) wzfatal("Unmatched signature. File corrupted.\n"); */
/*   fread(&cg->fmt, sizeof(char), 1, fh); */
/*   fread(&cg->n, sizeof(uint64_t), 1, fh); */
/*   cg->s = malloc(cgdata_nbytes(cg)); */
/*   fread(cg->s, 1, cgdata_nbytes(cg), fh); */
/*   /\* fprintf(stdout, "0s[0]: %u\n", cg->s[0]); *\/ */
/*   fclose(fh); */
/*   return cg; */
/* } */

typedef struct cgfile_t {
  gzFile fh;
  int n;                        /* number of samples read */
} cgfile_t;

static inline cgfile_t open_cgfile(const char *fname) {
  cgfile_t cgf = {0};
  cgf.fh = gzopen(fname, "rb");
  cgf.n = 0;
  return cgf;
}

static inline cgdata_t read_cg(cgfile_t *cgf) {
  
  cgdata_t cg = {0};

  uint64_t sig;
  if(!gzfread(&sig, sizeof(uint64_t), 1, cgf->fh)) return cg;
  
  if (sig != CGSIG) wzfatal("Unmatched signature. File corrupted.\n");
  gzfread(&cg.fmt, sizeof(char), 1, cgf->fh);
  gzfread(&cg.n, sizeof(uint64_t), 1, cgf->fh);
  cg.s = malloc(cgdata_nbytes(&cg));
  gzfread(cg.s, 1, cgdata_nbytes(&cg), cgf->fh);
  cgf->n++;
  /* fprintf(stdout, "0s[0]: %u\n", cg.s[0]); */
  return cg;
}

/* this function is memory intensive if there are many samples */
static inline cgdata_v* read_cg_all(cgfile_t *cgf) {

  cgdata_v *cgs = init_cgdata_v(10);
  while (1) {
    cgdata_t cg = read_cg(cgf);
    if (cg.n >0) push_cgdata_v(cgs, cg);
    else break;
  }
  return cgs;
}

#endif /* _TBMATE_H */
