#include "kycg.h"

/* 8 bit for 8 cpgs, each is binary */
cgdata_t* fmt0_read_uncompressed(char *fname, int verbose) {

  gzFile fh = wzopen(fname);
  char *line = NULL;
  uint64_t n = 0, m=1<<22;
  uint8_t *s = calloc(m, 1);
  while (gzFile_read_line(fh, &line) > 0) {
    if (line[0] != '0') {
      s[n>>3] |= (1<<(n&0x7));
    }
    n++;
    if (n+2>m) { m<<=1; s=realloc(s,m); }
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
  cg->fmt = '0';
  return cg;
}

void fmt0_decompress(cgdata_t *cg, cgdata_t *expanded) {
  expanded->s = realloc(expanded->s, cg->n>>3);
  memcpy(expanded->s, cg->s, cg->n>>3);
  expanded->n = cg->n;
  expanded->compressed = 0;
  expanded->fmt = '0';
}


