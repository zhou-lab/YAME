#include "kycg.h"

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
  cg->fmt = '0';
  return cg;
}

/* void process_0binaryVec(char *fname, char *fname_out, int verbose) { */
/*   uint8_t *s; uint64_t n = read_binary(fname, &s, verbose); */
/*   write_bytevec(fname_out, s, (n>>3)+1, '0', verbose, "bit vector"); */
/*   free(s); */
/* } */

