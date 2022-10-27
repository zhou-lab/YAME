#include <string.h>
#include <zlib.h>
#include <stdio.h>
#include "wzmisc.h"
#include "wzbed.h"
#include "kycg.h"

static int is_float(char *s) {
  size_t i;
  for (i=0; i<strlen(s); ++i) {
    if (!isdigit(s[i]) && s[i] != '.' && s[i] != '-')
      return 0;
  }
  return 1;
}

cgdata_t* fmt4_read_uncompressed(char *fname, int verbose) {

  gzFile fh = wzopen(fname);
  char *line = NULL;
  uint64_t n = 0, m=1<<22;
  float *s = calloc(m, sizeof(float));
  while (gzFile_read_line(fh, &line) > 0) {
    if (is_float(line)) {
      s[n++] = atof(line);
    } else {
      s[n++] = -1.0;
    }
    if (n+2>m) { m<<=1; s=realloc(s, m*sizeof(float)); }
  }
  free(line);
  wzclose(fh);
  if (verbose) {
    fprintf(stderr, "[%s:%d] Vector of length %zu loaded\n", __func__, __LINE__, n);
    fflush(stderr);
  }
  cgdata_t *cg = calloc(sizeof(cgdata_t),1);
  cg->s = (uint8_t*) s;
  cg->n = n;
  cg->compressed = 0;
  cg->fmt = '4';
  return cg;
}

/* 32 bit
   1 (1bit) + run length of NA (31 bits)
   0 (1bit) + floating number (always positive) (31bit, the sign bit is always 0)
 */
void fmt4_compress(cgdata_t *cg) {

  uint64_t n=0, m=1<<20;
  uint32_t *s = calloc(sizeof(uint32_t), m);
  uint64_t i = 0; uint32_t l = 0;
  uint32_t *s0 = (uint32_t*) cg->s;
  for (i=0, l=0; i<cg->n; ++i) {
    /* either not the same as before or reach block size max */
    if (!(s0[i] & (1ul<<31)) || l+2 >= (1ul<<31)) {
      if (l > 0) {
        if (n+2>m) { m<<=1; s = realloc(s, m*sizeof(uint32_t));}
        s[n++] = ((1<<31) | l);
        l = 0;
      }

      if (!(s0[i] & (1ul<<31))) {
        if (n+2>m) { m<<=1; s = realloc(s, m*sizeof(uint32_t));}
        memcpy(s+n, s0+i, sizeof(float_t));
        n++;
        l = 0;
      }
    } else {
      ++l;
    }
  }
  /* the last rle */
  if (l > 0) {
    if (n+2>m) { m<<=1; s = realloc(s, m*sizeof(uint32_t));}
    s[n++] = ((1<<31) | l);
  }
  fprintf(stdout, "%u\n", s[0]);
  fprintf(stdout, "%u\n", s[1]);
  
  free(cg->s);
  cg->s = (uint8_t*) s;
  cg->n = n*4;
  cg->compressed = 1;
}

cgdata_t fmt4_decompress(cgdata_t *cg) {

  uint64_t i=0, m = 1<<20,n = 0, j=0, l=0;
  uint32_t *s0 = (uint32_t*) cg->s;
  float_t *s = calloc(m, sizeof(float_t));

  for(i=0; i< cg->n>>2; ++i) {
    if (s0[i] >> 31) {
      l = s0[i]<<1>>1;
      if (n+l+10>m) {m=n+l+10; m<<=1; s = realloc(s, m*sizeof(float_t));}
      for (j=0; j<l; ++j) s[n++] = -1.0;
    } else {
      if (n+2>m) {m<<=1; s = realloc(s, m*sizeof(float_t));}
      memcpy(s+n, s0+i, sizeof(float_t));
      n++;
    }
  }

  cgdata_t cg2 = {0};
  cg2.s = (uint8_t*) s;
  cg2.n = n;
  cg2.compressed = 0;
  cg2.fmt = '4';
  return cg2;
}
