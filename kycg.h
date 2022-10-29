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
void fmt0_decompress(cgdata_t *cg, cgdata_t *expanded);

cgdata_t *fmt1_read_uncompressed(char *fname, int verbose);
void fmt1_compress(cgdata_t *cg);
void fmt1_decompress(cgdata_t *cg, cgdata_t *expanded);

cgdata_t *fmt3_read_uncompressed(char *fname, int verbose);
void fmt3_compress(cgdata_t *cg);
void fmt3_decompress(cgdata_t *cg, cgdata_t *expanded);

cgdata_t *fmt4_read_uncompressed(char *fname, int verbose);
void fmt4_compress(cgdata_t *cg);
void fmt4_decompress(cgdata_t *cg, cgdata_t *expanded);

cgdata_t *fmt5_read_uncompressed(char *fname, int verbose);
void fmt5_compress(cgdata_t *cg);
void fmt5_decompress(cgdata_t *cg, cgdata_t *expanded);

cgdata_t *fmt6_read_uncompressed(char *fname, int verbose);
void fmt6_compress(cgdata_t *cg);
void fmt6_decompress(cgdata_t *cg, cgdata_t *expanded);

void decompress(cgdata_t *cg, cgdata_t *expanded);
void recompress(cgdata_t *cg);

static inline void slice(cgdata_t *cg, uint64_t beg, uint64_t end, cgdata_t *cg_sliced) {
  if (end > cg->n-1) end = cg->n-1;
  if (end < beg) wzfatal("Slicing negative span.");

  cg_sliced->s = realloc(cg_sliced->s, (end-beg+1)*sizeof(uint8_t));
  memcpy(cg_sliced->s, cg->s+beg, (end-beg+1)*sizeof(uint8_t));
  cg_sliced->n = end - beg + 1;
  cg_sliced->compressed = 0;
  cg_sliced->fmt = '5';
}

static inline void cgdata_write(char *fname_out, cgdata_t *cg, const char *mode, int verbose) {

  FILE *out;
  if (fname_out)
    out = fopen(fname_out, mode);
  else
    out = stdout;
  
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

static inline int read_cg_(cgfile_t *cgf, cgdata_t *cg) {
  uint64_t sig;
  if(!gzfread(&sig, sizeof(uint64_t), 1, cgf->fh)) return 0;
  if (sig != CGSIG) wzfatal("Unmatched signature. File corrupted.\n");
  gzfread(&(cg->fmt), sizeof(char), 1, cgf->fh);
  gzfread(&(cg->n), sizeof(uint64_t), 1, cgf->fh);
  cg->s = realloc(cg->s, cgdata_nbytes(cg));
  gzfread(cg->s, 1, cgdata_nbytes(cg), cgf->fh);
  cg->compressed = 1;
  cgf->n++;
  return 1;
}

static inline cgdata_t read_cg(cgfile_t *cgf) {
  cgdata_t cg = {0};
  if (!read_cg_(cgf, &cg)) return cg;
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

static inline cgdata_v* read_cgs(cgfile_t *cgf, int64_t beg, int64_t end) {

  if (beg < 0) beg = 0;
  if (end >= 0 && end < beg) wzfatal("End is smaller than beg");

  cgdata_v *cgs = init_cgdata_v(10);
  cgdata_t cg = {0};
  int64_t i=0;
  for (i=0; end<0 || i<=end; ++i) {
    read_cg_(cgf, &cg);
    if (i<beg) continue;
    if (cg.n>0) {
      (*next_ref_cgdata_v(cgs)) = cg;
      cg.s = NULL;
    } else {
      break;
    }
  }
  return cgs;
}

#endif /* _TBMATE_H */
