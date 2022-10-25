#include "kycg.h"

static int is_nonnegative_int(char *s) {
  size_t i;
  for (i=0; i<strlen(s); ++i) {
    if (!isdigit(s[i])) return 0;
  }
  return 1;
}


/* unpacked
   [ M (uint32_t) | U (uint32_t) ]
   packed
   [ block_type (2bit) |
   2byte | U=M=0 -------------- = 0 2bit + run len (14 bit)
   1byte | U,M in [0,7] ------- = 1 2bit + M (3bit) | U (3bit)
   2byte | U,M in [0,127]------ = 2 2bit + M (7bit) | U (7bit)
   8byte | M,U in [128,2**31]-- = 3 2bit + M (31bit) | U (31bit) */
cgdata_t* fmt3_read_uncompressed(char *fname, int verbose) {
  gzFile fh = wzopen(fname);
  char *line = NULL;
  uint64_t n = 0, m=1<<10;
  uint64_t *s = calloc(m, sizeof(uint64_t));
  char **fields; int nfields;
  int M,U;
  while (gzFile_read_line(fh, &line) > 0) {
    line_get_fields(line, '\t', &fields, &nfields);
    if (nfields < 2) wzfatal("Number of fields <2. Abort.");
    if (!is_nonnegative_int(fields[0]) || !is_nonnegative_int(fields[1]))
      wzfatal("Field 1 or 2 is not a nonnegative integer.");
    M = atol(fields[0]);
    U = atol(fields[1]);
    (*s)[n++] = M<<32 | U;
    if (n>m-2) { m<<=1; *s=realloc(*s,m); }
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

void fmt3_compress(cgdata_t *cg) {
  uint8_t *s = NULL;
  uint64_t n = 0;
  uint64_t i = 0;
  uint64_t *s0 = (uint64_t*) cg->s;
  int M,U; uint32_t l=0;
  for (i=0; i<cg->n; ++i) {
    M = s0[i]>>32;
    U = s0[i]<<32>>32;
    if (M>0 || U>0 || l+2 >= 1<<14) {
      if (l>0) {
        *s = realloc(*s, n+2);
        *((uint16_t*) (s+n)) = (uint16_t) l;
        n += 2;
        l = 0;
      }
      if (M>0 || U>0) {
        if (M<7 && U<7) {
          *s = realloc(*s, n+1);
          s[n] = (1<<6) | (M<<3) | U;
          n++;
        } else if (M<127 && U<127) {
          *s = realloc(*s, n+2);
          *((uint16_t*) (s+n)) = (2<<14) | (M<<7) | U;
          n += 2;
        } else {
          /* cap the counts, not the best solution */
          if (M>(1<<31)) M = (1<<31) - 1;
          if (U>(1<<31)) U = (1<<31) - 1;
          *s = realloc(*s, n+8);
          *((uint64_t*) (s+n)) = (3<<62) | (M<<31) | U;
          n += 8;
        }
      }
    } else {
      ++l;
    }
  }
  if (l>0) {
    *s = realloc(*s, n+2);
    *((uint16_t*) (s+n)) = (uint16_t) l;
  }
  free(cg->s);
  cg->s = s;
  cg->n = n;
  cg->compressed = 1;
}
