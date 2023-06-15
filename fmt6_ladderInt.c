#include "kycg.h"

/* uncompressed: [ M (uint32_t) | U (uint32_t) ] */
cgdata_t* fmt6_read_uncompressed(char *fname, int verbose) {
  gzFile fh = wzopen(fname, 1);
  char *line = NULL;
  uint64_t n = 0, m=1<<10;
  uint64_t *s = calloc(m, sizeof(uint64_t));
  while (gzFile_read_line(fh, &line) > 0) {
    if (n+2>m) { m<<=1; s=realloc(s,m*sizeof(uint64_t)); }
    s[n++] = atol(line);
  }
  free(line);
  wzclose(fh);
  if (verbose) {
    fprintf(stderr, "[%s:%d] Vector of length %lu loaded\n", __func__, __LINE__, n);
    fflush(stderr);
  }
  cgdata_t *cg = calloc(sizeof(cgdata_t),1);
  cg->s = (uint8_t*) s;
  cg->n = n;
  cg->compressed = 0;
  cg->fmt = '6';
  return cg;
}

/* compressed: assume little endian, TODO: use endian swap
   1byte | L in [0,64] ------ = L (6bit) + 1 (2bit)
   2byte | L in [64,2^14]---- = L (14bit) + 2 (2bit)
   8byte | L in [2^14,2^62]-- =  L (62bit) + 3 (2bit)
*/
void fmt6_compress(cgdata_t *cg) {
  uint8_t *s = NULL;
  uint64_t n = 0;
  uint64_t i = 0;
  uint64_t *s0 = (uint64_t*) cg->s;
  for (i=0; i<cg->n; ++i) {
    uint64_t L = s0[i];
    if (L < (1<<6)-2) {
      s = realloc(s, n+1);
      s[n] = (L<<2) | 0x1;
      n++;
    } else if (L < (1<<14)-2) {
      s = realloc(s, n+2);
      *((uint16_t*) (s+n)) = (L<<2) | 0x2;
      n += 2;
    } else {
      /* cap the counts, not the best solution */
      if (L > (1ul<<62)) L = (1ul<<62) - 1;
      s = realloc(s, n+8);
      *((uint64_t*) (s+n)) = (L<<62) | 3ul;
      n += 8;
    }
  }
  free(cg->s);
  cg->s = s;
  cg->n = n;
  cg->compressed = 1;
}

void fmt6_decompress(cgdata_t *cg, cgdata_t *expanded) {
  uint64_t i = 0, m = 1<<20,n = 0;
  uint64_t *s = realloc(expanded->s, m*sizeof(uint64_t));
  while (i < cg->n) {
    /* if ((cg->s[i] & 0x3) == 0) { */
    /*   l = (((uint16_t*) (cg->s+i))[0])>>2; */
    /*   if (n+l+10>m) {m=n+l+10; m<<=1; s = realloc(s, m*sizeof(uint64_t));} */
    /*   for (j=0; j<l; ++j) s[n++] = 0; */
    /*   i += 2; */
    if ((cg->s[i] & 0x3) == 1) {
      if (n+2>m) {m<<=1; s = realloc(s, m*sizeof(uint64_t));}
      s[n++] = ((cg->s[i])>>2);
      i++;
    } else if ((cg->s[i] & 0x3) == 2) {
      if (n+2>m) {m<<=1; s = realloc(s, m*sizeof(uint64_t));}
      s[n++] = (((uint16_t*) (cg->s+i))[0])>>2;
      i += 2;
    } else {
      if (n+2>m) {m<<=1; s = realloc(s, m*sizeof(uint64_t));}
      s[n++] = (((uint64_t*) (cg->s+i))[0])>>2;
      i += 8;
    }
  }
  expanded->s = (uint8_t*) s;
  expanded->n = n;
  expanded->compressed = 0;
  expanded->fmt = '6';
}
