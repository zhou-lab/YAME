#include "kycg.h"

void fmta_tryBinary2byteRLE_ifsmaller(cgdata_t *cg) {

  uint64_t n = 0;
  uint8_t *s = NULL;
  uint64_t i=0; uint16_t l=0; uint8_t u0=0;
  for (i=0, l=0; i<cg->n; ++i) {
    /* unsigned char u = s[i>>3] & (1<<(n&0x7)); */
    uint8_t u = (s[i>>3]>>(i&0x7))&0x1;
    /* either not the same as before or reach block size max */
    if ((l != 0 && u != u0) || l+2 >= 1<<15) {
      s = realloc(s, n+3);
      s[n] = u0;
      *((uint16_t*) (s+n+1)) = l;
      n += 3;
      l = 1;
    } else {
      ++l;
    }
    u0 = u;
  }
  /* the last rle */
  s = realloc(s, n+3);
  s[n] = u0;
  *((uint16_t*) (s+n+1)) = l;
  n += 3;
  
  if (cg->n>>3 > n) {
    free(cg->s);
    cg->s = s;
    cg->n = n;
    cg->fmt = '1';
    cg->compressed = 1;
  }
}


