#include "kycg.h"

/* https://developer.ibm.com/articles/au-endianc/ */
#define is_bigendian() ( (∗(char∗)&i) == 0 )

static int is_nonnegative_int(char *s) {
  size_t i;
  for (i=0; i<strlen(s); ++i) {
    if (!isdigit(s[i])) return 0;
  }
  return 1;
}


/* uncompressed: [ M (uint32_t) | U (uint32_t) ] */
cgdata_t* fmt3_read_uncompressed(char *fname, int verbose) {
  gzFile fh = wzopen(fname);
  char *line = NULL;
  uint64_t n = 0, m=1<<10;
  uint64_t *s = calloc(m, sizeof(uint64_t));
  char **fields; int nfields;
  uint64_t M,U;
  while (gzFile_read_line(fh, &line) > 0) {
    line_get_fields(line, "\t", &fields, &nfields);
    if (nfields < 2) wzfatal("Number of fields <2. Abort.");
    if (!is_nonnegative_int(fields[0]) || !is_nonnegative_int(fields[1]))
      wzfatal("Field 1 or 2 is not a nonnegative integer.");
    M = atol(fields[0]);
    U = atol(fields[1]);
    if (n+2>m) { m<<=1; s=realloc(s,m*sizeof(uint64_t)); }
    s[n++] = ((M<<32) | U);
    free_fields(fields, nfields);
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
  cg->fmt = '3';
  return cg;
}

/* compressed: assume little endian, TODO: use endian swap
   2byte | U=M=0 -------------- = run len (14 bit) + 0 2bit
   1byte | U,M in [0,7] ------- = M (3bit) | U (3bit) + 1 2bit
   2byte | U,M in [0,127]------ = M (7bit) | U (7bit) + 2 2bit
   8byte | M,U in [128,2**31]-- = M (31bit) | U (31bit) + 3 2bit
*/
void fmt3_compress(cgdata_t *cg) {
  uint8_t *s = NULL;
  uint64_t n = 0;
  uint64_t i = 0;
  uint64_t *s0 = (uint64_t*) cg->s;
  uint64_t M,U; uint32_t l=0;
  for (i=0; i<cg->n; ++i) {
    M = s0[i]>>32;
    U = s0[i]<<32>>32;
    if (M>0 || U>0 || l+2 >= 1<<14) {
      if (l>0) {
        s = realloc(s, n+2);
        *((uint16_t*) (s+n)) = (uint16_t) l<<2;
        n += 2;
        if (M>0 || U>0) l = 0;
        else l = 1;
      }
      if (M>0 || U>0) {
        if (M<7 && U<7) {
          s = realloc(s, n+1);
          s[n] = (M<<5) | (U<<2) | 0x1;
          n++;
        } else if (M<127 && U<127) {
          s = realloc(s, n+2);
          *((uint16_t*) (s+n)) = (M<<9) | (U<<2) | 0x2;
          n += 2;
        } else {
          /* cap the counts, not the best solution */
          if (M > (1ul<<31)) M = (1ul<<31) - 1;
          if (U > (1ul<<31)) U = (1ul<<31) - 1;
          s = realloc(s, n+8);
          *((uint64_t*) (s+n)) = (M<<33) | (U<<31) | 3ul;
          n += 8;
        }
      }
    } else {
      ++l;
    }
  }
  if (l>0) {
    s = realloc(s, n+2);
    *((uint16_t*) (s+n)) = (uint16_t) l<<2;
    n += 2;
  }
  free(cg->s);
  cg->s = s;
  cg->n = n;
  cg->compressed = 1;
}

cgdata_t fmt3_decompress(cgdata_t *cg) {
  uint64_t i = 0, m = 1<<20,n = 0, j=0, l=0;
  uint64_t *s = calloc(m, sizeof(uint64_t));
  uint64_t U=0,M=0;
  while (i < cg->n) {
    if ((cg->s[i] & 0x3) == 0) {
      l = (((uint16_t*) (cg->s+i))[0])>>2;
      if (n+l+10>m) {m=n+l+10; m<<=1; s = realloc(s, m*sizeof(uint64_t));}
      for (j=0; j<l; ++j) s[n++] = 0;
      i += 2;
    } else if ((cg->s[i] & 0x3) == 1) {
      M = (cg->s[i])>>5;
      U = ((cg->s[i])>>2) & 0x7;
      if (n+2>m) {m<<=1; s = realloc(s, m*sizeof(uint64_t));}
      s[n++] = (M<<32 | U);
      i++;
    } else if ((cg->s[i] & 0x3) == 2) {
      M = (((uint16_t*) (cg->s+i))[0])>>2;
      U = M & ((1<<7)-1);
      M = (M>>7) & ((1<<7)-1);
      if (n+2>m) {m<<=1; s = realloc(s, m*sizeof(uint64_t));}
      s[n++] = (M<<32 | U);
      i += 2;
    } else {
      M = (((uint64_t*) (cg->s+i))[0])>>2;
      U = M & ((1ul<<31)-1);
      M = (M>>7) & ((1ul<<31)-1);
      if (n+2>m) {m<<=1; s = realloc(s, m*sizeof(uint64_t));}
      s[n++] = (M<<32 | U);
      i += 8;
    }
  }
  cgdata_t cg2 = {0};
  cg2.s = (uint8_t*) s;
  cg2.n = n;
  cg2.compressed = 0;
  cg2.fmt = '3';
  return cg2;
}
