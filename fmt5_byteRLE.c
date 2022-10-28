#include "kycg.h"


/* the input has only 0,1,2 */
cgdata_t* fmt5_read_uncompressed(char *fname, int verbose) {

  gzFile fh = wzopen(fname);
  char *line = NULL;
  uint64_t n = 0, m=1<<22;
  uint8_t *s = calloc(m, 1);
  while (gzFile_read_line(fh, &line) > 0) {
    s[n++] = line[0]-'0';
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
  cg->fmt = '5';
  return cg;
}

/*
  8 bits = 0 (1bit) | run length of NA (7bits)
  8 bits = 1 (1bit)|value (1bit) + 1 (1bit)|value (1bit) + ...
 */
void fmt5_compress(cgdata_t *cg) {
  uint64_t n = 0;
  uint8_t *s = NULL;
  uint64_t i = 0; uint16_t l = 0; int last = 0; uint8_t u = 0; int offset = 6;
  for (i=0, l=0; i<cg->n; ++i) {
    if (cg->s[i] == 0 || cg->s[i] == 1) { /* 0 or 1 */
      u |= (1<<(offset+1));
      u |= (cg->s[i]<<offset);
      offset -= 2;
      if (last <= 1) {            /* 0/1 > 0/1 */
        if (offset < 0) {
          s = realloc(s, n+1);
          s[n++] = u;
          u = 0; offset = 6;
        }
      } else if (l > 0) {       /* 2 > 0/1 */
        s = realloc(s, n+1);
        s[n++] = l;
        l = 0;
      }
      last = 1;
    } else {                    /* neither 0 nor 1, for missing value */
      if (last == 1 && u != 0) {               /* 0/1 > 2 */
        if (offset >= 0) u |= (0<<(offset+1)); /* add sentinel */
        s = realloc(s, n+1);
        s[n++] = u;
        u = 0; offset = 6;
      }
      l++;
      if (l+2 >= 1<<7) {        /* too many NA start a new count */
        s = realloc(s, n+1);
        s[n++] = l;
        l = 0;
      }
      last = 2;
    }
  }

  if (last == 1 && u != 0) {
    s = realloc(s, n+1);
    s[n++] = u;
  } else if (last == 2 && l > 0) {
    s = realloc(s, n+1);
    s[n++] = l;
  }

  free(cg->s);
  cg->s = s;
  cg->n = n;
  cg->compressed = 1;
}

cgdata_t fmt5_decompress(cgdata_t *cg) {
  uint64_t i = 0, m = 1<<20,n = 0, j=0;
  uint8_t *s = calloc(m, sizeof(uint8_t));

  for (i=0; i<cg->n; ++i) {
    if (cg->s[i] & (1<<7)) {
      int offset = 6;
      if (n+2>m) {m<<=1; s = realloc(s, m*sizeof(uint8_t));}
      for (offset = 6; offset >= 0; offset -= 2) {
        if ((cg->s[i]>>offset) & 0x2) {
          s[n++] = ((cg->s[i]>>offset) & 0x1);
        } else {
          break;
        }
      }
    } else {
      if (n+cg->s[i]+10>m) {m=n+cg->s[i]+10; m<<=1; s = realloc(s, m*sizeof(uint8_t));}
      for (j=0; j < cg->s[i]; ++j) s[n++] = 2;
    }
  }

  cgdata_t cg2 = {0};
  cg2.s = (uint8_t*) s;
  cg2.n = n;
  cg2.compressed = 0;
  cg2.fmt = '5';
  return cg2;
}
