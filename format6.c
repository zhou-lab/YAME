#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "cdata.h"

static int is_int(char *s) {
  size_t i;
  for (i=0; i<strlen(s); ++i) {
    if (!isdigit(s[i])) return 0;
  }
  return 1;
}

/* every byte is SUSUSUSU */
cdata_t* fmt6_read_raw(char *fname, int verbose) {

  gzFile fh = wzopen(fname, 1);
  char *line = NULL;
  uint8_t *s = NULL; uint64_t n = 0;
  char **fields; int nfields;
  while (gzFile_read_line(fh, &line) > 0) {
    line_get_fields(line, "\t", &fields, &nfields);
    if (nfields < 2) wzfatal("Number of fields < 2. Abort.");
    if (!is_int(fields[1])) wzfatal("The 2nd column must be integers.");
    s = realloc(s, n/4 + 1);
    // Binarize fields[0] to the first bit and fields[1] to the second bit
    if (n%4 == 0) s[n/4] = 0; // initialize
    if (atoi(fields[1])) s[n/4] |= 1<<((n%4)*2 + 1); // universe
    if (strcmp(fields[0], "1") == 0) s[n/4] |= 1<<((n%4)*2); // set
    n++;
    free_fields(fields, nfields);
  }
  free(line);
  wzclose(fh);
  if (verbose) {
    fprintf(stderr, "[%s:%d] Data of length %lu loaded\n", __func__, __LINE__, n);
    fflush(stderr);
  }
  cdata_t *c = calloc(sizeof(cdata_t),1);
  c->s = s;
  c->n = n;
  c->compressed = 0;
  c->fmt = '6';
  c->unit = 2;
  return c;
}

// do nothing
void fmt6_compress(cdata_t *c) {
  c->compressed = 1;
}

// just copy, nothing else, one can just flip compressed bit if possible
void fmt6_decompress(cdata_t *c, cdata_t *expanded) {
  expanded->s = realloc(expanded->s, cdata_nbytes(c));
  memcpy(expanded->s, c->s, cdata_nbytes(c));
  expanded->n = c->n;
  expanded->compressed = 0;
  expanded->fmt = '6';
  expanded->unit = 2;
}

/* uint64_t fmt6_count_universe(cdata_t c) { */

/*   /\* create a look-up table *\/ */
/*   int byte2cnt[256]; int p; */
/*   for (p=0; p<256; ++p) { */
/*     unsigned char q = p; int ii, cnt = 0; */
/*     for (ii=0; ii<4; ++ii) { */
/*       if (q & 0x2) cnt++; q>>=2; // universe is the higher bit */
/*     } */
/*     byte2cnt[p] = cnt; */
/*   } */
  
/*   size_t i,k,m = 0; */
/*   for (i=0; i<(c.n>>2); ++i) m += byte2cnt[c.s[i]];   // full bytes */
/*   for (k=0; k<(c.n&0x3); ++k) m += (c.s[i]>>((k*2)+1)) & 0x1; // last byte */
/*   return m; */
/* } */

/* uint64_t fmt6_count_query(cdata_t c) { */

/*   /\* create a look-up table *\/ */
/*   int byte2cnt[256]; int p; */
/*   for (p=0; p<256; ++p) { */
/*     unsigned char q = p; int ii, cnt = 0; */
/*     for (ii=0; ii<4; ++ii) { */
/*       if ((q & 0x2) && (q & 0x1)) cnt++; q>>=2; // both in query and universe */
/*     } */
/*     byte2cnt[p] = cnt; */
/*   } */
  
/*   size_t i,k,m = 0; */
/*   for (i=0; i<(c.n>>2); ++i) m += byte2cnt[c.s[i]]; // full bytes */
/*   for (k=0; k<(c.n&0x3); ++k) {                     // last byte */
/*     if (((c.s[i]>>(k*2)) & 0x3) == 0x3) m++; */
/*   } */
/*   return m; */
/* } */
