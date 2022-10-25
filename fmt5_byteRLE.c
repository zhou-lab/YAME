#include "kycg.h"


cgdata_t* fmt5_read_uncompressed(char *fname, int verbose) {

  gzFile fh = wzopen(fname);
  char *line = NULL;
  uint64_t n = 0, m=1<<22;
  uint8_t *s = calloc(m, 1);
  while (gzFile_read_line(fh, &line) > 0) {
    s[n++] = line[0]-'0';
    if (n>m-2) { m<<=1; s=realloc(s,m); }
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
  cg->fmt = '5';
  return cg;
}


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

void process_5byteRLE(char *fname, char *fname_out, int verbose) {
  uint8_t *s; uint64_t n = read_byteVector(fname, &s, verbose);
  uint8_t *sr; uint64_t nr = byteVec2RLE(s, n, &sr);
  write_bytevec(fname_out, sr, nr, '5', verbose, "byte RLE vector 5");
  free(s); free(sr);
  if (verbose) {
    fprintf(stderr, "[%s:%d] N_CpGs: %"PRId64"; N_bytes: %"PRId64"\n", __func__, __LINE__, n, nr);
    fflush(stderr);
  }
}
