#include "cgdata.h"

/* https://developer.ibm.com/articles/au-endianc/ */
#define is_bigendian() ( (∗(char∗)&i) == 0 )

static int is_nonnegative_int(char *s) {
  size_t i;
  for (i=0; i<strlen(s); ++i) {
    if (!isdigit(s[i])) return 0;
  }
  return 1;
}

static int is_float(char *s) {
  size_t i;
  for (i=0; i<strlen(s); ++i) {
    if (!isdigit(s[i]) && s[i] != '.' && s[i] != '-')
      return 0;
  }
  return 1;
}

/* 8 bit for 8 cpgs, each is binary */
cgdata_t* fmt0_read_uncompressed(char *fname, int verbose) {

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
  cgdata_t *cg = calloc(sizeof(cgdata_t),1);
  cg->s = (uint8_t*) s;
  cg->n = n;
  cg->compressed = 0;
  cg->fmt = '0';
  return cg;
}

/* just copy, nothing done */
void fmt0_decompress(cgdata_t *cg, cgdata_t *expanded) {
  expanded->s = realloc(expanded->s, cgdata_nbytes(cg));
  memcpy(expanded->s, cg->s, cgdata_nbytes(cg));
  expanded->n = cg->n;
  expanded->compressed = 0;
  expanded->fmt = '0';
}

/**
 * Reads data from a given cgfile, converts it into format 0 if needed, and returns it.
 * Format 0 refers to a bit-packed format where 1 byte represents an 8-bit binary vector.
 *
 * @param cgf: A pointer to a cgfile_t structure from which data is to be read.
 * @return: A cgdata_t structure which contains the read and converted data.
 * 
 * The function performs the following operations:
 * 
 * 1. It begins by reading data from the given cgfile into a cgdata_t structure.
 * 2. Input format is '0': no further operation is performed.
 * 3. Input format is '1': if the value is 0, return 0 else 1.
 * 4. Input format is '3': if the M+U is 0, return 0 else 1
 * 5. Other input formats are not allowed.
 */
void convertToFmt0(cgdata_t *cg) {
  cgdata_t cg_out = {0};
  switch (cg->fmt) {
  case '0': return;
  case '1': {
    cg_out.fmt = '0';
    cg_out.compressed = 1;
    cg_out.n=0;
    uint64_t i;
    for (i=0; i<cg->n/3; ++i) {
      cg_out.n += *((uint16_t*) (cg->s+i*3+1));
    }
    cg_out.s = calloc((cg_out.n>>3)+1, 1);
    size_t sum; uint16_t l=0;
    for (i=0, sum=0; i<cg->n/3; ++i, sum+=l) {
      l = *((uint16_t*) (cg->s+i*3+1));
      if (cg->s[i*3] > '0') {
        for(size_t j=sum; j<sum+l; ++j) {
          cg_out.s[j>>3] |= (1<<(j&0x7));
        }
      }
    }
    break;
  }
  case '3': {
    cgdata_t expanded = {0};
    fmt3_decompress(cg, &expanded);

    cg_out.fmt = '0';
    cg_out.compressed = 1;
    cg_out.n = expanded.n;
    cg_out.s = calloc((cg_out.n>>3)+1,1);
    uint64_t *s = (uint64_t*) expanded.s;
    for (uint64_t i=0; i<expanded.n; ++i) {
      if (s[i] > 0) { /* equivalent to: 1 if M+U > 0 else 0 */
        cg_out.s[i>>3] |= (1<<(i&0x7));
      }
    }
    free(expanded.s);
    break;
  }
  default: wzfatal("Format %c unsupported.\n", cg->fmt);
  }
  free(cg->s);
  *cg = cg_out;
}

cgdata_t* fmt1_read_uncompressed(char *fname, int verbose) {

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
  cgdata_t *cg = calloc(sizeof(cgdata_t),1);
  cg->s = (uint8_t*) s;
  cg->n = n;
  cg->compressed = 0;
  cg->fmt = '1';
  return cg;
}

/* compressed:
   3byte --- value (1byte) + run len (2bytes)
   value is unrestricted ASCII
 */
void fmt1_compress(cgdata_t *cg) {
  uint64_t n = 0;
  uint8_t *s = NULL;
  uint64_t i = 0; uint16_t l = 0; uint8_t u0 = 0;
  for (i=0, l=0; i<cg->n; ++i) {
    /* either not the same as before or reach block size max */
    if ((l != 0 && cg->s[i] != u0) || l+2 >= 1<<15) {
      s = realloc(s, n+3);
      s[n] = u0;
      *((uint16_t*) (s+n+1)) = l;
      n += 3;
      l = 1;
    } else {
      ++l;
    }
    u0 = cg->s[i];
  }
  /* the last rle */
  s = realloc(s, n+3);
  s[n] = u0;
  *((uint16_t*) (s+n+1)) = l;
  n += 3;
  
  free(cg->s);
  cg->s = s;
  cg->n = n;
  cg->compressed = 1;
}

void fmt1_decompress(cgdata_t *cg, cgdata_t *expanded) {
  uint64_t i=0, j=0, n=0, m=1<<20;
  uint8_t *s = realloc(expanded->s, m);
  for (i=0; i<cg->n; i+=3) {
    uint16_t l = ((uint16_t*) (cg->s+i+1))[0];
    if (n+l+2>m) {m=n+l+2; m<<=1; s = realloc(s, m);}
    for (j=0; j<l; ++j) s[n++] = cg->s[i];
  }
  expanded->s = (uint8_t*) s;
  expanded->n = n;
  expanded->compressed = 0;
  expanded->fmt = '1';
}

/* uncompressed: [ M (uint32_t) | U (uint32_t) ] */
cgdata_t* fmt3_read_uncompressed(char *fname, int verbose) {
  gzFile fh = wzopen(fname, 1);
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
    fprintf(stderr, "[%s:%d] Vector of length %lu loaded\n", __func__, __LINE__, n);
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
   2byte | U=M=0 -------------- = run len (14 bit) + 0 (2bit)
   1byte | U,M in [0,7] ------- = M (3bit) | U (3bit) + 1 (2bit)
   2byte | U,M in [0,127]------ = M (7bit) | U (7bit) + 2 (2bit)
   8byte | M,U in [128,2**31]-- = M (31bit) | U (31bit) + 3 (2bit)
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
          *((uint64_t*) (s+n)) = (M<<33) | (U<<2) | 3ul;
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

void fmt3_decompress(cgdata_t *cg, cgdata_t *expanded) {
  uint64_t i = 0, m = 1<<20,n = 0, j=0, l=0;
  uint64_t *s = realloc(expanded->s, m*sizeof(uint64_t));
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
      M >>= 7;
      if (n+2>m) {m<<=1; s = realloc(s, m*sizeof(uint64_t));}
      s[n++] = (M<<32 | U);
      i += 2;
    } else {
      M = (((uint64_t*) (cg->s+i))[0])>>2;
      U = M & ((1ul<<31)-1);
      M >>= 31;
      if (n+2>m) {m<<=1; s = realloc(s, m*sizeof(uint64_t));}
      s[n++] = (M<<32 | U);
      i += 8;
    }
  }
  expanded->s = (uint8_t*) s;
  expanded->n = n;
  expanded->compressed = 0;
  expanded->fmt = '3';
}

cgdata_t* fmt4_read_uncompressed(char *fname, int verbose) {

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
  cgdata_t *cg = calloc(sizeof(cgdata_t),1);
  cg->s = (uint8_t*) s;
  cg->n = n;
  cg->compressed = 0;
  cg->fmt = '4';
  return cg;
}

/* 32 bit
   1 (1bit) + run length of NA (31 bits)
   0 (1bit) + floating number (always positive) (31bit, the sign bit is always 0)
 */
void fmt4_compress(cgdata_t *cg) {

  uint64_t n=0, m=1<<20;
  uint32_t *s = calloc(sizeof(uint32_t), m);
  uint64_t i = 0; uint32_t l = 0;
  uint32_t *s0 = (uint32_t*) cg->s;
  for (i=0, l=0; i<cg->n; ++i) {
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
  
  free(cg->s);
  cg->s = (uint8_t*) s;
  cg->n = n*4;
  cg->compressed = 1;
}

void fmt4_decompress(cgdata_t *cg, cgdata_t *expanded) {

  uint64_t i=0, m = 1<<20,n = 0, j=0, l=0;
  uint32_t *s0 = (uint32_t*) cg->s;
  float_t *s = realloc(expanded->s, m*sizeof(float_t));

  for(i=0; i< cg->n>>2; ++i) {
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
}

/* the input has only 0,1,2 */
cgdata_t* fmt5_read_uncompressed(char *fname, int verbose) {

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
  cgdata_t *cg = calloc(sizeof(cgdata_t),1);
  cg->s = (uint8_t*) s;
  cg->n = n;
  cg->compressed = 0;
  cg->fmt = '5';
  return cg;
}

/*
  8 bits = 0 (1bit) | run length of NA (7bits)
  8 bits = 1 (1bit)|value (1bit) + 1 (1bit)|value (1bit) + ...
 */
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

void fmt5_decompress(cgdata_t *cg, cgdata_t *expanded) {
  uint64_t i = 0, m = 1<<20,n = 0, j=0;
  uint8_t *s = realloc(expanded->s, m*sizeof(uint8_t));

  for (i=0; i<cg->n; ++i) {
    if (cg->s[i] & (1<<7)) {
      int offset = 6;
      if (n+2>m) {m<<=1; s = realloc(s, m*sizeof(uint8_t));}
      for (offset = 6; offset >= 0; offset -= 2) {
        if ((cg->s[i]>>offset) & 0x2) {
          s[n++] = ((cg->s[i]>>offset) & 0x1);
        } else {
          break;
        }
      }
    } else {
      if (n+cg->s[i]+10>m) {m=n+cg->s[i]+10; m<<=1; s = realloc(s, m*sizeof(uint8_t));}
      for (j=0; j < cg->s[i]; ++j) s[n++] = 2;
    }
  }

  expanded->s = (uint8_t*) s;
  expanded->n = n;
  expanded->compressed = 0;
  expanded->fmt = '5';
}

/* uncompressed: [ M (uint32_t) | U (uint32_t) ] */
cgdata_t* fmt6_read_uncompressed(char *fname, int verbose) {
  gzFile fh = wzopen(fname, 1);
  char *line = NULL;
  uint64_t n = 0, m=1<<10;
  uint64_t *s = calloc(m, sizeof(uint64_t));
  while (gzFile_read_line(fh, &line) > 0) {
    if (n+2>m) { m<<=1; s=realloc(s,m*sizeof(uint64_t)); }
    s[n++] = atol(line);
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
  cg->fmt = '6';
  return cg;
}

/* compressed: assume little endian, TODO: use endian swap
   1byte | L in [0,64] ------ = L (6bit) + 1 (2bit)
   2byte | L in [64,2^14]---- = L (14bit) + 2 (2bit)
   8byte | L in [2^14,2^62]-- =  L (62bit) + 3 (2bit)
*/
void fmt6_compress(cgdata_t *cg) {
  uint8_t *s = NULL;
  uint64_t n = 0;
  uint64_t i = 0;
  uint64_t *s0 = (uint64_t*) cg->s;
  for (i=0; i<cg->n; ++i) {
    uint64_t L = s0[i];
    if (L < (1<<6)-2) {
      s = realloc(s, n+1);
      s[n] = (L<<2) | 0x1;
      n++;
    } else if (L < (1<<14)-2) {
      s = realloc(s, n+2);
      *((uint16_t*) (s+n)) = (L<<2) | 0x2;
      n += 2;
    } else {
      /* cap the counts, not the best solution */
      if (L > (1ul<<62)) L = (1ul<<62) - 1;
      s = realloc(s, n+8);
      *((uint64_t*) (s+n)) = (L<<62) | 3ul;
      n += 8;
    }
  }
  free(cg->s);
  cg->s = s;
  cg->n = n;
  cg->compressed = 1;
}

void fmt6_decompress(cgdata_t *cg, cgdata_t *expanded) {
  uint64_t i = 0, m = 1<<20,n = 0;
  uint64_t *s = realloc(expanded->s, m*sizeof(uint64_t));
  while (i < cg->n) {
    /* if ((cg->s[i] & 0x3) == 0) { */
    /*   l = (((uint16_t*) (cg->s+i))[0])>>2; */
    /*   if (n+l+10>m) {m=n+l+10; m<<=1; s = realloc(s, m*sizeof(uint64_t));} */
    /*   for (j=0; j<l; ++j) s[n++] = 0; */
    /*   i += 2; */
    if ((cg->s[i] & 0x3) == 1) {
      if (n+2>m) {m<<=1; s = realloc(s, m*sizeof(uint64_t));}
      s[n++] = ((cg->s[i])>>2);
      i++;
    } else if ((cg->s[i] & 0x3) == 2) {
      if (n+2>m) {m<<=1; s = realloc(s, m*sizeof(uint64_t));}
      s[n++] = (((uint16_t*) (cg->s+i))[0])>>2;
      i += 2;
    } else {
      if (n+2>m) {m<<=1; s = realloc(s, m*sizeof(uint64_t));}
      s[n++] = (((uint64_t*) (cg->s+i))[0])>>2;
      i += 8;
    }
  }
  expanded->s = (uint8_t*) s;
  expanded->n = n;
  expanded->compressed = 0;
  expanded->fmt = '6';
}

void fmta_tryBinary2byteRLE_ifsmaller(cgdata_t *cg) {
  uint64_t n = 0;
  uint8_t *s = NULL;
  uint64_t i=0; uint16_t l=0; uint8_t u0=0;
  for (i=0, l=0; i<cg->n; ++i) {
    /* unsigned char u = s[i>>3] & (1<<(n&0x7)); */
    uint8_t u = (cg->s[i>>3]>>(i&0x7))&0x1;
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
  if (cg->n>>3 > n) {
    free(cg->s);
    cg->s = s;
    cg->n = n;
    cg->fmt = '1';
    cg->compressed = 1;
  }
}


