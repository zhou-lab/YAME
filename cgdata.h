#ifndef _CGDATA_H
#define _CGDATA_H

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

#define CGSIG 266563789635

typedef struct f2_aux_t {
  uint64_t nk;                  // num keys
  char **keys;                  // pointer to keys
  uint8_t *data;                // pointer to data
} f2_aux_t;

/** The header design, 17 bytes
    uint64_t: signature, used for validation
    uint8_t: format (0=vec; 1=rle)
    uint64_t: length (n_cgs or n_bytes for rle)
**/
typedef struct cgdata_t {
  uint8_t *s;
  uint64_t n; /* number of bytes, except for fmt 0, which is sub-byte you need the actual length */
  int compressed;
  char fmt;
  uint8_t unit; // how many bytes is needed for each decompressed data unit
  void *aux;
} cgdata_t;

static inline uint64_t cgdata_nbytes(cgdata_t *cg) {
  uint64_t n = 0;
  switch(cg->fmt) {
  case '0': n = (cg->n>>3)+1; break;
  default: n = cg->n;
  }
  return n;
}

/* /\* unit size of uncompressed data *\/ */
/* static inline uint64_t cgdata_unit_size(cgdata_t *cg) { */
/*   switch(cg->fmt) { */
/*   case '2': return 8; break; // TODO: should fix */
/*   case '3': return 8; break; */
/*   case '4': return 4; break; */
/*   case '5': return 1; break; */
/*   case '6': return 8; break; */
/*   default: return 1; */
/*   } */
/*   return 1; */
/* } */

static inline void free_cgdata(cgdata_t *cg) {
  if(cg->s) free(cg->s);
  if (cg->fmt == '2' && cg->aux) {
    free(((f2_aux_t*) cg->aux)->keys);
    free(cg->aux);
  }
  cg->s = NULL;
}

void fmta_tryBinary2byteRLE_ifsmaller(cgdata_t *cg);

void decompress(cgdata_t *cg, cgdata_t *expanded);
void decompress2(cgdata_t *cg);
void cdata_compress(cgdata_t *cg);

void fmt0_decompress(cgdata_t *cg, cgdata_t *inflated);

void fmt1_compress(cgdata_t *cg);
void fmt1_decompress(cgdata_t *cg, cgdata_t *inflated);

// ----- format 2 (state data) ----
// key section + data section
// The key section and data section are separated by an extra '\0'.
// The key section is made of multiple c-strings concatenated by '\0'.
// The data section is either an RLE (compressed) or a integer vector (inflated).
// When compressed, the RLE is made of a value part and a length part.
// The value part size is defined by a uint8_t that leads the data section.
// The length part is always 2 bytes in size.
void fmt2_compress(cgdata_t *cg);
void fmt2_decompress(cgdata_t *cg, cgdata_t *inflated);
void fmt2_set_aux(cgdata_t *cg);
uint8_t* fmt2_get_data(cgdata_t *cg);
uint64_t fmt2_get_keys_n(cgdata_t *cg);

void fmt3_compress(cgdata_t *cg);
void fmt3_decompress(cgdata_t *cg, cgdata_t *inflated);
uint64_t f3_unpack_mu(cgdata_t *cg, uint64_t i);

void fmt4_compress(cgdata_t *cg);
void fmt4_decompress(cgdata_t *cg, cgdata_t *inflated);

void fmt5_compress(cgdata_t *cg);
void fmt5_decompress(cgdata_t *cg, cgdata_t *inflated);

void fmt6_compress(cgdata_t *cg);
void fmt6_decompress(cgdata_t *cg, cgdata_t *inflated);

void convertToFmt0(cgdata_t *cg);

static inline void slice(cgdata_t *cg, uint64_t beg, uint64_t end, cgdata_t *cg_sliced) {

  if (cg->compressed) {
    fprintf(stderr, "[%s:%d] Cannot slice compressed data.\n", __func__, __LINE__);
    fflush(stderr);
    exit(1);
  }
  if (end > cg->n-1) end = cg->n-1;
  if (end < beg) wzfatal("Slicing negative span.");

  cg_sliced->s = realloc(cg_sliced->s, (end-beg+1)*cg->unit);
  memcpy(cg_sliced->s, cg->s+beg*cg->unit, (end-beg+1)*cg->unit);
  cg_sliced->n = end - beg + 1;
  cg_sliced->compressed = 0;
  cg_sliced->fmt = cg->fmt;
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

/* static inline void set_data_uint64(uint8_t *data, uint64_t value, uint8_t value_nbytes) { */
/*   for (uint8_t i=0; i<value_nbytes; ++i) { */
/*     data[i] = (value & 0xff); */
/*     value >>= 8; */
/*   } */
/* } */

static inline uint64_t f2_unpack_uint64(cgdata_t *cg, uint64_t i) {
  if (!cg->aux) fmt2_set_aux(cg);
  f2_aux_t *aux = (f2_aux_t*) cg->aux;
  uint8_t *d = aux->data + cg->unit*i;
  uint64_t value = 0;
  for (uint8_t j=0; j<cg->unit; ++j) value |= (d[j] << (8*j));
  return value;
}

static inline char* f2_unpack_string(cgdata_t *cg, uint64_t i) {
  if (!cg->aux) fmt2_set_aux(cg);
  f2_aux_t *aux = (f2_aux_t*) cg->aux;
  uint64_t val = f2_unpack_uint64(cg, i);
  if (val >= aux->nk) {
    fprintf(stderr, "[%s:%d] State data is corrupted.\n", __func__, __LINE__);
    fflush(stderr);
    exit(1);
  }
  return aux->keys[val];
}

#endif /* _CGDATA_H */
