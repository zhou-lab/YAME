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

// beg and end are both 0-based
void fmt7_sliceToBlock(cdata_t *cr, uint64_t beg, uint64_t end, cdata_t *cr2) {
  if (cr->fmt != '7') {
    fprintf(stderr, "[%s:%d] Expect format 7 but got %c.\n", __func__, __LINE__, cr->fmt);
    fflush(stderr);
    exit(1);
  }
  
  uint64_t n0 = fmt7_data_length(cr);
  if (end > n0-1) end = n0-1; // 0-base
  if (beg > n0-1) {
    fprintf(stderr, "[%s:%d] Begin (%"PRIu64") is bigger than the data vector size (%"PRIu64").\n", __func__, __LINE__, beg, n0);
    fflush(stderr);
    exit(1);
  }

  row_reader_t rdr = {0};
  uint64_t n = 0;
  uint64_t i = 0, n_rec = 0;
  char *chrm = NULL; uint64_t last = 0;
  memset(cr2, 0, sizeof(cdata_t));
  while (row_reader_next_loc(&rdr, cr)) {
    if (i >= beg && i <= end) {
      if (chrm != rdr.chrm) {
        if (chrm) append_end(&(cr2->s), &n);
        chrm = rdr.chrm;
        append_chrm(chrm, &(cr2->s), &n);
      }
      append_loc(rdr.value - last, &(cr2->s), &n);
      n_rec++;
      last = rdr.value;
    }
    i++;
  }
  if (n_rec != end - beg + 1) {
    fprintf(stderr, "[%s:%d] row slicing has inconsistent dimension (n: %"PRIu64", expected: %"PRIu64")\n", __func__, __LINE__, n_rec, end - beg + 1);
    fflush(stderr);
    exit(1);
  }
  cr2->unit = cr->unit;
  cr2->fmt = cr->fmt;
  cr2->n = n;
  cr2->compressed = 1;
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

int fmt7_next_bed(cdata_t *c) {
  row_reader_t *rdr;
  if (!c->aux) c->aux = calloc(1, sizeof(row_reader_t));
  rdr = (row_reader_t*) c->aux;
  return row_reader_next_loc(rdr, c);
}

// no decompression, just calculate N and prepare row_reader
void fmt7_decompress(cdata_t *c, cdata_t *inflated) {
  inflated->s = malloc(c->n);
  memcpy(inflated->s, c->s, c->n);
  inflated->unit = 1;
  inflated->aux = calloc(1,sizeof(row_reader_t));
  inflated->compressed = 0;
  inflated->fmt = '7';
  inflated->n = c->n;
}

uint64_t fmt7_data_length(cdata_t *c) {
  row_reader_t rdr = {0};
  uint64_t n = 0;
  while (row_reader_next_loc(&rdr, c)) n++;
  return n;
}
