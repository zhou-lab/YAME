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
    if (n>m-2) { m<<=1; s=realloc(s, m*sizeof(float)); }
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

void fmt4_compress(cgdata_t *cg) {

  uint64_t n = 0;
  uint8_t *s = NULL;
  uint64_t i = 0; uint32_t l = 0;
  for (i=0, l=0; i<cg->n; ++i) {
    /* either not the same as before or reach block size max */
    if (cg->s[i] > 0 || l+2 >= 1<<15) {
      if (l > 0) {
        s = realloc(s, n+5);
        s[n] = 0;
        *((uint32_t*) (s+n+1)) = l;
        n += 5;
        l = 0;
      }

      if (cg->s[i] > 0) {
        s = realloc(s, n+5);
        s[n] = 1;
        ((float*)(s+n+1))[0] = cg->s[i];
        n += 5;
        l = 0;
      }
    } else {
      ++l;
    }
  }
  /* the last rle */
  if (l > 0) {
    s = realloc(s, n+5);
    s[n] = 0;
    *((uint32_t*) (s+n+1)) = l;
    n += 5;
  }
  
  free(cg->s);
  cg->s = s;
  cg->n = n;
  cg->compressed = 1;
}

void process_4fractionNA(char *fname, char *fname_out, int verbose) {
  float *s; uint64_t n = read_fractionVec(fname, &s, verbose);
  uint8_t *sr; uint64_t nr = fractionVec2RLE(s, n, &sr);
  write_bytevec(fname_out, sr, nr, '4', verbose, "fractionNA vector");
  free(s); free(sr);
  if (verbose) {
    fprintf(stderr, "[%s:%d] N: %zu; nr: %zu\n", __func__, __LINE__, n, nr);
    fflush(stderr);
  }
}
