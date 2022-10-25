#ifndef _TBMATE_H
#define _TBMATE_H

#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <limits.h>
#include <inttypes.h>
#include <wordexp.h>
#include "wzmisc.h"
#include "wzmisc.h"
#include "wzbed.h"

#define PACKAGE_VERSION "0.1.20221020"
#define CGSIG 266563789635

typedef struct cgdata_t {
  uint8_t *s;
  uint64_t n; /* number of bytes, except for fmt 0, which is sub-byte you need the actual length*/
  int compressed;
  char fmt;
} cgdata_t;

static inline uint64_t cgdata_nbytes(cgdata_t *cg) {
  switch(cg->fmt) {
  case '0': n = (cg->n>>3)+1; break;
  default: n = cg->n;
  }
  return n;
}

static inline free_cgdata(cgdata_t *cg) {
  if(cg->s) free(cg->s);
  free(cg);
}

cgdata_t *fmt3_read_uncompressed(char *fname, int verbose);
void fmt3_compress(cgdata_t *cg);

static inline void cgdata_write(char *fname_out, cgdata_t *cg, int verbose) {

  FILE *out = fopen(fname_out, "wb");
  uint64_t sig = CGSIG; fwrite(&sig, sizeof(uint64_t), 1, out);
  fwrite(&fmt, sizeof(uint8_t), 1, out);
  fwrite(&cg->n, sizeof(uint64_t), 1, out);
  fwrite(cg->s, 1, cgdata_nbytes(cg), out);
  fclose(out);

  if (verbose) {
    fprintf(stderr, "[%s:%d] Stored as Format %c\n", __func__, __LINE__, cg->fmt);
    fflush(stderr);
  }
}

#endif /* _TBMATE_H */
