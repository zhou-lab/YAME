#include "kycg.h"

cgdata_t* fmt1_read_uncompressed(char *fname, int verbose) {

  gzFile fh = wzopen(fname);
  char *line = NULL;
  uint64_t n = 0, m=1<<22;
  uint8_t *s = calloc(m, 1);
  while (gzFile_read_line(fh, &line) > 0) {
    s[n++] = line[0];
    if (n+2>m) { m<<=1; s=realloc(s,m); }
  }
  free(line);
  wzclose(fh);
  if (verbose) {
    fprintf(stderr, "[%s:%d] Vector of length %llu loaded\n", __func__, __LINE__, n);
    fflush(stderr);
  }
  cgdata_t *cg = calloc(sizeof(cgdata_t),1);
  cg->s = (uint8_t*) s;
  cg->n = n;
  cg->compressed = 0;
  cg->fmt = '1';
  return cg;
}

/* compressed:
   3byte --- value (1byte) + run len (2bytes)
   value is unrestricted ASCII
 */
void fmt1_compress(cgdata_t *cg) {
  uint64_t n = 0;
  uint8_t *s = NULL;
  uint64_t i = 0; uint16_t l = 0; uint8_t u0 = 0;
  for (i=0, l=0; i<cg->n; ++i) {
    /* either not the same as before or reach block size max */
    if ((l != 0 && cg->s[i] != u0) || l+2 >= 1<<15) {
      s = realloc(s, n+3);
      s[n] = u0;
      *((uint16_t*) (s+n+1)) = l;
      n += 3;
      l = 1;
    } else {
      ++l;
    }
    u0 = cg->s[i];
  }
  /* the last rle */
  s = realloc(s, n+3);
  s[n] = u0;
  *((uint16_t*) (s+n+1)) = l;
  n += 3;
  
  free(cg->s);
  cg->s = s;
  cg->n = n;
  cg->compressed = 1;
}

cgdata_t fmt1_decompress(cgdata_t *cg) {
  uint64_t i=0, j=0, n=0, m=1<<20;
  uint8_t *s = calloc(m, 1);
  for (i=0; i<cg->n; i+=3) {
    uint16_t l = ((uint16_t*) (cg->s+i+1))[0];
    if (n+l+2>m) {m=n+l+2; m<<=1; s = realloc(s, m);}
    for (j=0; j<l; ++j) s[n++] = cg->s[i];
  }
  cgdata_t cg2 = {0};
  cg2.s = (uint8_t*) s;
  cg2.n = n;
  cg2.compressed = 0;
  cg2.fmt = '1';

  return cg2;
}
