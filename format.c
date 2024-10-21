#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "cdata.h"

static int is_float(char *s) {
  size_t i;
  for (i=0; i<strlen(s); ++i) {
    if (!isdigit(s[i]) && s[i] != '.' && s[i] != '-')
      return 0;
  }
  return 1;
}

/* 8 bit for 8 cpgs, each is binary */
cdata_t* fmt0_read_raw(char *fname, int verbose) {

  gzFile fh = wzopen(fname, 1);
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
  cdata_t *c = calloc(sizeof(cdata_t),1);
  c->s = (uint8_t*) s;
  c->n = n;
  c->compressed = 0;
  c->fmt = '0';
  c->unit = 1;
  return c;
}

/* just copy, nothing done */
void fmt0_decompress(cdata_t *c, cdata_t *expanded) {
  expanded->s = realloc(expanded->s, cdata_nbytes(c));
  memcpy(expanded->s, c->s, cdata_nbytes(c));
  expanded->unit = 1;
  expanded->n = c->n;
  expanded->compressed = 0;
  expanded->fmt = '0';
}

/**
 * - Input format is '0': no further operation is performed.
 * - Input format is '1': if the value is 0, return 0 else 1.
 * - Input format is '3': if the M+U is 0, return 0 else 1
 * - Other input formats are not allowed.
 */
void convertToFmt0(cdata_t *c) {
  cdata_t c_out = {0};
  switch (c->fmt) {
  case '0': return;
  case '1': {
    c_out.fmt = '0';
    c_out.compressed = 1;
    c_out.n=0;
    uint64_t i;
    for (i=0; i<c->n/3; ++i) {
      c_out.n += *((uint16_t*) (c->s+i*3+1));
    }
    c_out.s = calloc((c_out.n>>3)+1, 1);
    size_t sum; uint16_t l=0;
    for (i=0, sum=0; i<c->n/3; ++i, sum+=l) {
      l = *((uint16_t*) (c->s+i*3+1));
      if (c->s[i*3] > '0') {
        for(size_t j=sum; j<sum+l; ++j) {
          FMT0_SET(c_out, j);
        }
      }
    }
    break;
  }
  case '3': {
    cdata_t expanded = {0};
    fmt3_decompress(c, &expanded);

    c_out.fmt = '0';
    c_out.compressed = 1;
    c_out.n = expanded.n;
    c_out.s = calloc((c_out.n>>3)+1,1);
    for (uint64_t i=0; i<expanded.n; ++i) {
      uint64_t mu = f3_get_mu(&expanded, i);
      if (mu > 0) { /* equivalent to: 1 if M+U > 0 else 0 */
        c_out.s[i>>3] |= (1<<(i&0x7));
      }
    }
    free(expanded.s);
    break;
  }
  default: wzfatal("Format %c unsupported.\n", c->fmt);
  }
  free(c->s);
  *c = c_out;
}

cdata_t* fmt1_read_raw(char *fname, int verbose) {

  gzFile fh = wzopen(fname, 1);
  char *line = NULL;
  uint64_t n = 0, m=1<<22;
  uint8_t *s = calloc(m, 1);
  while (gzFile_read_line(fh, &line) > 0) {
    s[n++] = line[0];
    if (n+2>m) { m<<=1; s=realloc(s,m); }
  }
  free(line);
  wzclose(fh);
  if (verbose) {
    fprintf(stderr, "[%s:%d] Vector of length %lu loaded\n", __func__, __LINE__, n);
    fflush(stderr);
  }
  cdata_t *c = calloc(sizeof(cdata_t),1);
  c->s = (uint8_t*) s;
  c->n = n;
  c->compressed = 0;
  c->fmt = '1';
  return c;
}

/* compressed:
   3byte --- value (1byte) + run len (2bytes)
   value is unrestricted ASCII
 */
void fmt1_compress(cdata_t *c) {
  uint64_t n = 0;
  uint8_t *s = NULL;
  uint64_t i = 0; uint16_t l = 0; uint8_t u0 = 0;
  for (i=0, l=0; i<c->n; ++i) {
    /* either not the same as before or reach block size max */
    if ((l != 0 && c->s[i] != u0) || l+2 >= 1<<15) {
      s = realloc(s, n+3);
      s[n] = u0;
      *((uint16_t*) (s+n+1)) = l;
      n += 3;
      l = 1;
    } else {
      ++l;
    }
    u0 = c->s[i];
  }
  /* the last rle */
  s = realloc(s, n+3);
  s[n] = u0;
  *((uint16_t*) (s+n+1)) = l;
  n += 3;
  
  free(c->s);
  c->s = s;
  c->n = n;
  c->compressed = 1;
}

void fmt1_decompress(cdata_t *c, cdata_t *expanded) {
  uint64_t i=0, j=0, n=0, m=1<<20;
  uint8_t *s = realloc(expanded->s, m);
  for (i=0; i<c->n; i+=3) {
    uint16_t l = ((uint16_t*) (c->s+i+1))[0];
    if (n+l+2>m) {m=n+l+2; m<<=1; s = realloc(s, m);}
    for (j=0; j<l; ++j) s[n++] = c->s[i];
  }
  expanded->s = (uint8_t*) s;
  expanded->n = n;
  expanded->compressed = 0;
  expanded->fmt = '1';
  expanded->unit = 1;
}

cdata_t* fmt4_read_raw(char *fname, int verbose) {

  gzFile fh = wzopen(fname, 1);
  char *line = NULL;
  uint64_t n = 0, m=1<<22;
  float *s = calloc(m, sizeof(float));
  while (gzFile_read_line(fh, &line) > 0) {
    if (is_float(line)) {
      s[n++] = atof(line);
    } else {
      s[n++] = -1.0;
    }
    if (n+2>m) { m<<=1; s=realloc(s, m*sizeof(float)); }
  }
  free(line);
  wzclose(fh);
  if (verbose) {
    fprintf(stderr, "[%s:%d] Vector of length %lu loaded\n", __func__, __LINE__, n);
    fflush(stderr);
  }
  cdata_t *c = calloc(sizeof(cdata_t),1);
  c->s = (uint8_t*) s;
  c->n = n;
  c->compressed = 0;
  c->fmt = '4';
  return c;
}

/* 32 bit
   1 (1bit) + run length of NA (31 bits)
   0 (1bit) + floating number (always positive) (31bit, the sign bit is always 0)
 */
void fmt4_compress(cdata_t *c) {

  uint64_t n=0, m=1<<20;
  uint32_t *s = calloc(sizeof(uint32_t), m);
  uint64_t i = 0; uint32_t l = 0;
  uint32_t *s0 = (uint32_t*) c->s;
  for (i=0, l=0; i<c->n; ++i) {
    /* either not the same as before or reach block size max */
    if (!(s0[i] & (1ul<<31)) || l+2 >= (1ul<<31)) {
      if (l > 0) {
        if (n+2>m) { m<<=1; s = realloc(s, m*sizeof(uint32_t));}
        s[n++] = ((1<<31) | l);
        l = 0;
      }

      if (!(s0[i] & (1ul<<31))) {
        if (n+2>m) { m<<=1; s = realloc(s, m*sizeof(uint32_t));}
        memcpy(s+n, s0+i, sizeof(float_t));
        n++;
        l = 0;
      }
    } else {
      ++l;
    }
  }
  /* the last rle */
  if (l > 0) {
    if (n+2>m) { m<<=1; s = realloc(s, m*sizeof(uint32_t));}
    s[n++] = ((1<<31) | l);
  }
  
  free(c->s);
  c->s = (uint8_t*) s;
  c->n = n*4;
  c->compressed = 1;
}

void fmt4_decompress(cdata_t *c, cdata_t *expanded) {

  uint64_t i=0, m = 1<<20,n = 0, j=0, l=0;
  uint32_t *s0 = (uint32_t*) c->s;
  float_t *s = realloc(expanded->s, m*sizeof(float_t));

  for(i=0; i< c->n>>2; ++i) {
    if (s0[i] >> 31) {
      l = s0[i]<<1>>1;
      if (n+l+10>m) {m=n+l+10; m<<=1; s = realloc(s, m*sizeof(float_t));}
      for (j=0; j<l; ++j) s[n++] = -1.0;
    } else {
      if (n+2>m) {m<<=1; s = realloc(s, m*sizeof(float_t));}
      memcpy(s+n, s0+i, sizeof(float_t));
      n++;
    }
  }

  expanded->s = (uint8_t*) s;
  expanded->n = n;
  expanded->compressed = 0;
  expanded->fmt = '4';
  expanded->unit = 4;
}

/* the input has only 0,1,2 */
cdata_t* fmt5_read_raw(char *fname, int verbose) {

  gzFile fh = wzopen(fname, 1);
  char *line = NULL;
  uint64_t n = 0, m=1<<22;
  uint8_t *s = calloc(m, 1);
  while (gzFile_read_line(fh, &line) > 0) {
    if (line[0] == '0' || line[0] == '1') s[n++] = line[0]-'0';
    else s[n++] = 2;
    if (n+2>m) { m<<=1; s=realloc(s,m); }
  }
  free(line);
  wzclose(fh);
  if (verbose) {
    fprintf(stderr, "[%s:%d] Vector of length %lu loaded\n", __func__, __LINE__, n);
    fflush(stderr);
  }
  cdata_t *c = calloc(sizeof(cdata_t),1);
  c->s = (uint8_t*) s;
  c->n = n;
  c->compressed = 0;
  c->fmt = '5';
  return c;
}

/*
  8 bits = 0 (1bit) | run length of NA (7bits)
  8 bits = 1 (1bit)|value (1bit) + 1 (1bit)|value (1bit) + ...
 */
void fmt5_compress(cdata_t *c) {
  uint64_t n = 0;
  uint8_t *s = NULL;
  uint64_t i = 0; uint16_t l = 0; int last = 0; uint8_t u = 0; int offset = 6;
  for (i=0, l=0; i<c->n; ++i) {
    if (c->s[i] == 0 || c->s[i] == 1) { /* 0 or 1 */
      u |= (1<<(offset+1));
      u |= (c->s[i]<<offset);
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

  free(c->s);
  c->s = s;
  c->n = n;
  c->compressed = 1;
}

void fmt5_decompress(cdata_t *c, cdata_t *expanded) {
  uint64_t i = 0, m = 1<<20,n = 0, j=0;
  uint8_t *s = realloc(expanded->s, m*sizeof(uint8_t));

  for (i=0; i<c->n; ++i) {
    if (c->s[i] & (1<<7)) {
      int offset = 6;
      if (n+2>m) {m<<=1; s = realloc(s, m*sizeof(uint8_t));}
      for (offset = 6; offset >= 0; offset -= 2) {
        if ((c->s[i]>>offset) & 0x2) {
          s[n++] = ((c->s[i]>>offset) & 0x1);
        } else {
          break;
        }
      }
    } else {
      if (n+c->s[i]+10>m) {m=n+c->s[i]+10; m<<=1; s = realloc(s, m*sizeof(uint8_t));}
      for (j=0; j < c->s[i]; ++j) s[n++] = 2;
    }
  }

  expanded->s = (uint8_t*) s;
  expanded->n = n;
  expanded->compressed = 0;
  expanded->fmt = '5';
  expanded->unit = 1;
}

void fmta_tryBinary2byteRLE_ifsmaller(cdata_t *c) {
  uint64_t n = 0;
  uint8_t *s = NULL;
  uint64_t i=0; uint16_t l=0; uint8_t u0=0;
  for (i=0, l=0; i<c->n; ++i) {
    uint8_t u = (c->s[i>>3]>>(i&0x7))&0x1;
    /* either not the same as before or reach block size max */
    if ((l != 0 && u != u0) || l+2 >= 1<<15) {
      s = realloc(s, n+3);
      s[n] = u0+'0';
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
  s[n] = u0+'0';
  *((uint16_t*) (s+n+1)) = l;
  n += 3;
  if (c->n>>3 > n) {
    free(c->s);
    c->s = s;
    c->n = n;
    c->fmt = '1';
    c->compressed = 1;
  }
}


