#ifndef _CDATA_H
#define _CDATA_H

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
#include "bgzf.h"

#define CDSIG 266563789635

typedef struct f2_aux_t {
  uint64_t nk;                  // num keys
  char **keys;                  // pointer to keys, doesn't own memory
  uint8_t *data;                // pointer to data, doesn't own memory
} f2_aux_t;

/** The header design, 17 bytes
    uint64_t: signature, used for validation
    uint8_t: format (0=vec; 1=rle)
    uint64_t: length (n_cs or n_bytes for rle)
**/
typedef struct cdata_t {
  uint8_t *s;
  uint64_t n; /* number of bytes, except for fmt 0, which is sub-byte you need the actual length */
  int compressed;
  char fmt;
  uint8_t unit; // how many bytes is needed for each decompressed data unit
  void *aux;
} cdata_t;

static inline uint64_t cdata_nbytes(cdata_t *c) {
  uint64_t n = 0;
  switch(c->fmt) {
  case '0': n = (c->n>>3)+1; break;
  default: n = c->n;
  }
  return n;
}

static inline void free_cdata(cdata_t *c) {
  if (c->s) free(c->s);
  if (c->fmt == '2' && c->aux) {
    free(((f2_aux_t*) c->aux)->keys);
    free(c->aux);
  }
  c->s = NULL;
}

void fmta_tryBinary2byteRLE_ifsmaller(cdata_t *c);

void decompress(cdata_t *c, cdata_t *expanded);
void decompress2(cdata_t *c);
void cdata_compress(cdata_t *c);

static inline uint64_t cdata_n(cdata_t *c) {
  cdata_t c2 = {0};
  decompress(c, &c2);
  uint64_t n = c2.n;
  free(c2.s);
  return n;
}

void fmt0_decompress(cdata_t *c, cdata_t *inflated);

void fmt1_compress(cdata_t *c);
void fmt1_decompress(cdata_t *c, cdata_t *inflated);

// ----- format 2 (state data) ----
// key section + data section
// The key section and data section are separated by an extra '\0'.
// The key section is made of multiple c-strings concatenated by '\0'.
// The data section is either an RLE (compressed) or a integer vector (inflated).
// When compressed, the RLE is made of a value part and a length part.
// The value part size is defined by a uint8_t that leads the data section.
// The length part is always 2 bytes in size.
void fmt2_compress(cdata_t *c);
void fmt2_decompress(cdata_t *c, cdata_t *inflated);
void fmt2_set_aux(cdata_t *c);
uint8_t* fmt2_get_data(cdata_t *c);
uint64_t fmt2_get_keys_n(cdata_t *c);

void fmt3_compress(cdata_t *c);
void fmt3_decompress(cdata_t *c, cdata_t *inflated);
uint64_t f3_unpack_mu(cdata_t *c, uint64_t i);

void fmt4_compress(cdata_t *c);
void fmt4_decompress(cdata_t *c, cdata_t *inflated);

void fmt5_compress(cdata_t *c);
void fmt5_decompress(cdata_t *c, cdata_t *inflated);

void fmt6_compress(cdata_t *c);
void fmt6_decompress(cdata_t *c, cdata_t *inflated);

void convertToFmt0(cdata_t *c);

static inline void slice(cdata_t *c, uint64_t beg, uint64_t end, cdata_t *c_sliced) {

  if (c->compressed) {
    fprintf(stderr, "[%s:%d] Cannot slice compressed data.\n", __func__, __LINE__);
    fflush(stderr);
    exit(1);
  }
  if (end > c->n-1) end = c->n-1;
  if (end < beg) wzfatal("Slicing negative span.");

  c_sliced->s = realloc(c_sliced->s, (end-beg+1)*c->unit);
  memcpy(c_sliced->s, c->s+beg*c->unit, (end-beg+1)*c->unit);
  c_sliced->n = end - beg + 1;
  c_sliced->compressed = 0;
  c_sliced->fmt = c->fmt;
}

static inline uint64_t sumMUpair(uint64_t MU1, uint64_t MU2) {
  uint64_t M = (MU1>>32) + (MU2>>32);
  uint64_t U = (MU1&0xffffffff) + (MU2&0xffffffff);
  if (M > 0xffffffff || U > 0xffffffff) {
    uint64_t tmp;
    int im = 0; tmp=M; while(tmp>>32) { tmp>>=1; ++im; }
    int iu = 0; tmp=U; while(tmp>>32) { tmp>>=1; ++iu; }
    im = (im>iu ? im : iu);
    M>>=im; U>>=im;
  }
  return (M<<32|U);
}

static inline uint64_t MUbinarize(uint64_t MU) {
  uint64_t M = MU>>32;
  uint64_t U = MU&0xffffffff;
  if (M==0 && U==0) {
    return 0;
  } else if (M>=U) {
    return (1ul<<32);
  } else {
    return 1ul;
  }
}

static inline uint64_t f2_unpack_uint64(cdata_t *c, uint64_t i) {
  if (!c->aux) fmt2_set_aux(c);
  f2_aux_t *aux = (f2_aux_t*) c->aux;
  uint8_t *d = aux->data + c->unit*i;
  uint64_t value = 0;
  for (uint8_t j=0; j<c->unit; ++j) value |= (d[j] << (8*j));
  return value;
}

static inline char* f2_unpack_string(cdata_t *c, uint64_t i) {
  if (!c->aux) fmt2_set_aux(c);
  f2_aux_t *aux = (f2_aux_t*) c->aux;
  uint64_t val = f2_unpack_uint64(c, i);
  if (val >= aux->nk) {
    fprintf(stderr, "[%s:%d] State data is corrupted.\n", __func__, __LINE__);
    fflush(stderr);
    exit(1);
  }
  return aux->keys[val];
}

#endif /* _CDATA_H */
