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
    if (!is_int(fields[0]) || !is_int(fields[1])) wzfatal("First two columns must be integers.");
    s = realloc(s, n/4 + 1);
    // Binarize fields[0] to the first bit and fields[1] to the second bit
    s[n/4] |= (atoi(fields[0]) & 1) << (2*(n%4));
    s[n/4] |= (atoi(fields[1]) & 1) << (2*(n%4) + 1);
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
  expanded->unit = 2;
  expanded->n = c->n;
  expanded->compressed = 0;
  expanded->fmt = '0';
}
