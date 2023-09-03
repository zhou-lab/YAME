#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "cdata.h"

static int is_nonnegative_int(char *s) {
  size_t i;
  for (i=0; i<strlen(s); ++i) {
    if (!isdigit(s[i])) return 0;
  }
  return 1;
}

static void append_chrm(char *chrm, uint8_t **s, uint64_t *n) {
  int dn = strlen(chrm)+1;
  *s = realloc(*s, (*n) + dn);
  strcpy((char*) ((*s) + *n), chrm);
  *n += dn;
}

/* first byte:
  00 00 0000 - 1 byte, max 0x7f
  10 00 0000 - 2 bytes, max 0x3fff
  11 00 0000 - 8 bytes max (1<<62) - 1
*/
static void append_loc(uint64_t loc, uint8_t **s, uint64_t *n) {
  if (loc <= 0x7f) {
    *s = realloc(*s, (*n)+1);
    (*s)[*n] = (uint8_t) loc;
    (*n)++;
  } else if (loc <= 0x3fff) {
    *s = realloc(*s, (*n)+2);
    (*s)[*n] = (0x80 | (loc>>8));
    (*s)[(*n)+1] = (loc & 0xff);
    (*n) += 2;
  } else if (loc <= ((1ul<<62) - 1)) {
    *s = realloc(*s, (*n)+8);
    for (int i=0; i<8; ++i)
      (*s)[(*n)+i] = ((loc>>(8*(7-i)))&0xff);
    (*s)[*n] |= 0xc0;
    (*n) += 8;
  } else {
    fprintf(stderr, "[%s:%d] Inter-loci distance exceeds maximum: %"PRIu64"\n", __func__, __LINE__, loc);
    fflush(stderr);
    exit(1);
  }
}

static void append_end(uint8_t **s, uint64_t *n) {
  *s = realloc(*s, (*n)+1);
  (*s)[*n] = 0xff;
  (*n)++;
}

cdata_t *fmt7_read_raw(char *fname, int verbose) {
  gzFile fh = wzopen(fname, 1);
  char *line = NULL;
  uint8_t *s = NULL; uint64_t n = 0;
  char **fields; int nfields;
  char *chrm = NULL;
  uint64_t last = 0;
  while (gzFile_read_line(fh, &line)>0) {
    line_get_fields(line, "\t", &fields, &nfields);
    if (nfields < 2) wzfatal("Number of fields <2. Abort.");
    if (!is_nonnegative_int(fields[1]))
      wzfatal("Field 1 or 2 is not a nonnegative integer.");

    uint64_t loc = atol(fields[1])+1;
    if (!chrm || strcmp(chrm, fields[0]) != 0) {
      if (chrm) {
        append_end(&s, &n);
        free(chrm);
      }
      chrm = strdup(fields[0]);
      append_chrm(chrm, &s, &n);
      last = 0;
    }
    append_loc(loc-last, &s, &n);
    last = loc;
    free_fields(fields, nfields);
  }
  if (verbose) {
    fprintf(stderr, "[%s:%d] Vector of length %lu loaded\n", __func__, __LINE__, n);
    fflush(stderr);
  }
  
  if (chrm) free(chrm);
  cdata_t *c = calloc(sizeof(cdata_t),1);
  c->s = s;
  c->n = n;
  c->compressed = 1;
  c->fmt = '7';
  c->unit = 1;
  return c;
}

int row_reader_next_loc(row_reader_t *rdr, cdata_t *c) {
  if (rdr->loc >= c->n) return 0;
  if (c->s[rdr->loc] == 0xff || !rdr->index) {
    if (c->s[rdr->loc] == 0xff) rdr->loc++;
    rdr->chrm = (char*) c->s + rdr->loc;
    rdr->loc += strlen(rdr->chrm)+1;
    rdr->value = 0;
  }

  if ((c->s[rdr->loc]>>6) == 3) { // 8 bytes
    uint64_t dn = (((uint64_t) c->s[rdr->loc] & 0x3f)<<(8*7));
    for (int i=1; i<8; ++i)
      dn |= (((uint64_t) c->s[rdr->loc+i])<<(8*(7-i)));
    rdr->value += dn;
    rdr->loc += 8;
  } else if ((c->s[rdr->loc]>>6) == 2) { // 2 bytes
    uint64_t dn = (((uint64_t) c->s[rdr->loc] & 0x3f)<<8);
    dn |= (uint64_t) c->s[rdr->loc+1];
    rdr->value += dn;
    rdr->loc += 2;
  } else {                      // 1 byte
    rdr->value += c->s[rdr->loc]<<1>>1;
    rdr->loc++;
  }
  rdr->index++;
  return 1;
}

void fmt7_decompress(cdata_t *c, cdata_t *inflated) {
  row_reader_t rdr = {0};
  uint64_t n=0;
  while (row_reader_next_loc(&rdr, c)) n++;
  inflated->s = realloc(inflated->s, n*sizeof(uint64_t));
  uint64_t *s = (uint64_t*) inflated->s;
  memset(&rdr, 0, sizeof(row_reader_t));
  while (row_reader_next_loc(&rdr, c)) s[rdr.index-1] = rdr.value;
  
}
