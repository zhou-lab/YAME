// SPDX-License-Identifier: AGPL-3.0-or-later
/**
 * This file is part of YAME.
 *
 * Copyright (C) 2021-present Wanding Zhou
 *
 * YAME is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * YAME is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with YAME.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "cdata.h"
#include "summary.h"

/** ---- format 7 (genomic coordinates: chromosome + delta-encoded loci) -----
 *
 * Overview
 * --------
 * Format 7 stores a sparse list of genomic coordinates (typically CpG loci).
 * Coordinates are grouped by chromosome, and on disk each chromosome is
 * represented by its name followed by a stream of *delta-encoded* positions:
 *
 *     [ chrm_name '\0' ][ delta1 ][ delta2 ] ... [ 0xff ]
 *
 * where:
 *   - chrm_name is a null-terminated string (e.g. "chr1")
 *   - each delta encodes the *difference* (loc_i - loc_{i-1}) with variable
 *     byte width (1/2/8 bytes)
 *   - 0xff marks the end of a chromosome block
 *
 * The original (on-disk) byte stream is always in this compressed form and is
 * decoded sequentially by row_reader_next_loc().
 *
 *
 * Position and delta encoding (compressed representation)
 * -------------------------------------------------------
 * Let loc_i be the 1-based genomic coordinate of entry i.  Instead of storing
 * loc_i outright, we store dn_i = loc_i - loc_{i-1}.  dn_i is encoded as:
 *
 *   1-byte delta (0 ≤ dn ≤ 0x7f):
 *       [0bbb_bbbb]          (MSB two bits = 00)
 *
 *   2-byte delta (0x80 ≤ dn ≤ 0x3fff):
 *       [10xx_xxxx][yyyy_yyyy]
 *       top bits = 10, lower 14 bits carry dn
 *
 *   8-byte delta (dn ≤ (1<<62)-1):
 *       [11......][........] (total 8 bytes, MSB two bits = 11)
 *       remaining 62 bits store dn (big-endian in the byte stream)
 *
 * These rules match fmt7c_append_loc() and the decoding logic in
 * row_reader_next_loc().
 *
 *
 * Chromosome boundaries (compressed representation)
 * -------------------------------------------------
 * After all deltas for a chromosome:
 *
 *     0xff
 *
 * is appended as an end-of-chromosome marker.  Next comes the next chromosome
 * name + '\0', followed by its own delta stream.
 *
 *
 * Original representation vs. cdata_t (compressed)
 * ------------------------------------------------
 * A format-7 cdata_t loaded from disk (e.g. via fmt7_read_raw()) looks like:
 *
 *     c->fmt        = '7'
 *     c->compressed = 1
 *     c->unit       = 1
 *     c->s          = [chrom\0][delta...][0xff][chrom\0][delta...]...
 *
 * There is no fixed-width inflated layout on disk.  Iteration always proceeds
 * by delta decoding with row_reader_next_loc().
 *
 * fmt7_data_length(c):
 *     counts how many coordinates exist by scanning the compressed stream.
 *
 * fmt7_next_bed(c):
 *     like row_reader_next_loc(), but updates c->aux state for BED-style use.
 *
 *
 * Decompressed/indexed representation
 * ---------------------------------------
 * For downstream tasks that need random access or per-row indexing, we provide
 * a *decompressed/indexed* in-memory representation via fmt7_decompress().
 * This is an alternative layout stored directly in cdata_t::s:
 *
 *     [chr1\0chr2\0...chrN\0\0][entry_0][entry_1]...[entry_{n-1}]
 *
 * where:
 *
 *   1) Chromosome name block:
 *        A sequence of NUL-terminated chromosome name strings, ending with an
 *        extra '\0'.  The "\0\0" double NUL sentinel marks the end of the
 *        name list:
 *
 *           "chr1\0chr2\0...\0chrN\0\0"
 *
 *   2) Packed coordinate entries:
 *        Each CpG/row is stored as a fixed 8-byte record:
 *
 *           byte 0..1 : uint16_t chr_id        (index into the name list)
 *           byte 2..7 : 48-bit coordinate      (beg1, little-endian)
 *
 *        Conceptually:
 *
 *           [ chr_id(2 bytes) | loc(6 bytes) ]
 *
 * fmt7_decompress():
 *     - Iterates over the original compressed stream using
 *       row_reader_next_loc().
 *     - Assigns each chromosome a sequential ID (0, 1, 2, ...).
 *     - Packs (chr_id, full_coordinate) for each CpG into the 8-byte layout.
 *     - Returns a cdata_t with:
 *
 *           c->fmt        = '7'
 *           c->compressed = 0
 *           c->unit       = 8
 *           c->s          = [names...][packed entries...]
 *           c->n          = cpg number (not the byte number)
 *
 * IMPORTANT:
 *     This “decompressed” format is not a traditional per-base inflation of
 *     the genome.  It is an alternate, index-friendly representation of the
 *     *same sparse coordinate list*, with fixed-width records and an explicit
 *     chromosome-name table.  It does not contain deltas or 0xff markers.
 *
 *
 * Accessing the decompressed representation
 * -----------------------------------------
 * The packed entries are 8-byte records that must be decoded via helper
 * macros/functions, not by casting to uint64_t* (the 48-bit coordinate is
 * explicitly stored in little-endian byte order).
 *
 * Example (conceptual):
 *
 *     uint16_t chr_id;
 *     uint64_t loc;
 *     FMT7_GET_LOC(c->s, i, &chr_id, &loc);   // decode entry i
 *
 * The chromosome name corresponding to chr_id is looked up from the leading
 * name block.
 *
 *
 * Relationship between compressed and decompressed forms
 * ------------------------------------------------------
 * • fmt7_read_raw()       : produces the compressed/delta representation  
 * • row_reader_next_loc() : decodes compressed format 7 sequentially  
 * • fmt7_data_length()    : counts CpGs by scanning (works on compressed)  
 * • fmt7_decompress()     : builds an indexed, fixed-width in-memory layout  
 *
 * Decompressed format-7 is intended for in-memory querying, slicing, and
 * summarization. The on-disk format remains the delta-encoded stream.
 *
 *
 * Slicing utilities (compressed side)
 * -----------------------------------
 * fmt7_sliceToBlock(), fmt7_sliceToIndices(), fmt7_sliceToMask():
 *     - walk the decoded order of coordinates (via row_reader_next_loc())
 *     - rebuild a new format-7 byte stream containing only selected entries
 *     - always output cdata_t with fmt='7' and compressed=1 (delta-encoded)
 *
 */

static void fmt7c_append_chrm(char *chrm, uint8_t **s, uint64_t *n) {
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
static void fmt7c_append_loc(uint64_t loc, uint8_t **s, uint64_t *n) {
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

static void fmt7c_append_end(uint8_t **s, uint64_t *n) {
  *s = realloc(*s, (*n)+1);
  (*s)[*n] = 0xff;
  (*n)++;
}

static int is_nonnegative_int(char *s) {
  size_t i;
  for (i=0; i<strlen(s); ++i) {
    if (!isdigit(s[i])) return 0;
  }
  return 1;
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
    if (!chrm || strcmp(chrm, fields[0]) != 0 ||
        loc < last) { // if unsorted, this will treat as a new chromosome.
      if (chrm) {
        fmt7c_append_end(&s, &n);
        free(chrm);
      }
      chrm = strdup(fields[0]);
      fmt7c_append_chrm(chrm, &s, &n);
      last = 0;
    }
    fmt7c_append_loc(loc-last, &s, &n);
    last = loc;
    free_fields(fields, nfields);
  }
  if (verbose) {
    fprintf(stderr, "[%s:%d] Vector of length %"PRIu64" loaded\n", __func__, __LINE__, n);
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

/**
 * row_reader_next_loc()
 * ----------------------
 * Sequentially decode the next coordinate from a compressed fmt7 cdata_t.
 *
 * This function is the core decoder for the delta-encoded (compressed)
 * format-7 representation. It advances a row_reader_t cursor (`rdr`)
 * across the compressed byte stream `c->s`, reconstructing one genomic
 * coordinate at a time.
 *
 * Format-7 compressed layout recap
 * --------------------------------
 * The compressed representation stored in `cdata_t::s` has this structure:
 *
 *   [ "chr1\0" ][ delta1 ][ delta2 ] ... [ 0xff ]
 *   [ "chr2\0" ][ deltas ... ] ...
 *
 * • chromosome name = NUL-terminated string
 * • each delta is encoded as 1, 2, or 8 bytes:
 *
 *     1 byte: 00xxxxxx                (dn  = 0 .. 0x7f)
 *     2 byte: 10xxxxxx yyyy yyyy      (dn  = 0x80 .. 0x3fff)
 *     8 byte: 11xxxxxx [7 more bytes] (dn  = up to 2^62−1)
 *
 * • dn = delta = loc_i – loc_{i−1}
 * • A single byte 0xff marks end-of-chromosome.
 *
 * Behavior of this decoder
 * ------------------------
 * Given a row_reader_t with fields:
 *
 *     rdr->loc   : byte offset into c->s for the next object to decode
 *     rdr->index : 1-based record index (starts at 0)
 *     rdr->value : current genomic coordinate
 *     rdr->chrm  : current chromosome name
 *
 * this function:
 *
 *  1. Detects start-of-chromosome or end-of-chromosome via:
 *        - rdr->index == 0     (first call)
 *        - c->s[rdr->loc] == 0xff (end-of-chr marker)
 *
 *  2. Loads `rdr->chrm` from c->s (a pointer into the data buffer).
 *
 *  3. Reads one delta (1/2/8 bytes), adds it to `rdr->value`, and advances
 *     `rdr->loc` accordingly.
 *
 *  4. Increments `rdr->index` and returns 1.
 *
 *  5. Returns 0 only when rdr->loc has passed `c->n`, i.e. end-of-stream.
 *
 *
 * Return value
 * ------------
 *   1 : successfully decoded the next coordinate
 *   0 : end of full fmt7 stream (no more chromosomes)
 *
 *
 * Usage example
 * -------------
 * Iterate through all coordinates in a compressed fmt7 vector:
 *
 *     cdata_t *c = fmt7_read_raw("coords.bed.gz", 0);  // compressed format-7
 *     row_reader_t rdr = {0};                          // initialize cursor
 *
 *     while (row_reader_next_loc(&rdr, c)) {
 *         printf("%s\t%" PRIu64 "\t(index=%" PRIu64 ")\n",
 *                rdr.chrm, rdr.value, rdr.index);
 *     }
 *
 * Typical output:
 *
 *     chr1    100    (index=1)
 *     chr1    105    (index=2)
 *     chr1    220    (index=3)
 *     chr2     50    (index=4)
 *     chr2     55    (index=5)
 *
 * Notes
 * -----
 * • `rdr->chrm` is always a pointer *into c->s* → do not free it.
 * • Decoding is strictly sequential; random access requires fmt7_decompress().
 * • `rdr->value` resets to 0 whenever a new chromosome begins.
 * • The caller is expected to initialize rdr = {0} before the first call.
 */
int row_reader_next_loc(row_reader_t *rdr, const cdata_t *c) {
  if (rdr->loc >= c->n) return 0; // past chromosome length
  if (c->s[rdr->loc] == 0xff || !rdr->index) { // hit end of chromosome
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

/**
 * Return the CpG number from iterating the whole cdata_t
 */
uint64_t fmt7_data_length(const cdata_t *c) {
  row_reader_t rdr = {0};
  uint64_t n = 0;
  while (row_reader_next_loc(&rdr, c)) n++;
  return n;
}

/**
 * ---------------------------------------------------------------------------
 *  FMT7_SET_LOC / FMT7_GET_LOC
 * ---------------------------------------------------------------------------
 * These helper macros encode and decode a single fixed-width coordinate
 * record in the *decompressed/indexed* fmt7 representation.
 *
 * Background
 * ----------
 * After fmt7_decompress(), the in-memory layout of an uncompressed format-7
 * vector is:
 *
 *     [chr1\0chr2\0 ... chrN\0\0][entry_0][entry_1]...[entry_{n-1}]
 *
 * where each entry_i is an 8-byte packed structure:
 *
 *     byte 0..1 : uint16_t chr_id       (index into the chromosome name list)
 *     byte 2..7 : uint48  loc           (full genomic coordinate, beg1)
 *
 * Both fields are encoded in **little-endian** byte order.  The coordinate
 * uses only 48 bits and is masked to ensure safety:
 *
 *     loc & 0xFFFFFFFFFFFFULL
 *
 * These macros perform byte-level encoding/decoding.  They do *not* rely on
 * casting to uint64_t*, which avoids alignment requirements and is portable
 * across architectures.
 *
 *
 * FMT7_SET_LOC(buf, i, chr_id, loc)
 * ---------------------------------
 * Write the (chr_id, loc) pair into the i-th 8-byte slot of `buf`.
 *
 * Parameters:
 *     buf     : pointer to start of packed entry array (uint8_t*)
 *     i       : 0-based index of the entry to write
 *     chr_id  : 16-bit chromosome index
 *     loc     : 48-bit genomic coordinate (beg1)
 *
 * Layout written:
 *
 *     buf[i*8 + 0] = chr_id low byte
 *     buf[i*8 + 1] = chr_id high byte
 *     buf[i*8 + 2..7] = loc in little-endian order
 *
 * FMT7_SET_LOC never allocates; it assumes buf has at least (i+1)*8 bytes.
 *
 *
 * FMT7_GET_LOC(buf, i, chr_id_out, loc_out)
 * -----------------------------------------
 * Read the (chr_id, loc) pair stored at the i-th slot of `buf`.
 *
 * Parameters:
 *     buf        : pointer to packed entry array (uint8_t*)
 *     i          : 0-based entry index
 *     chr_id_out : pointer to uint16_t to receive chromosome ID
 *     loc_out    : pointer to uint64_t to receive full coordinate
 *
 * Reads:
 *
 *     chr_id = buf[i*8 + 0] + (buf[i*8 + 1] << 8)
 *     loc    = buf[i*8 + 2] + (buf[i*8 + 3] << 8) + ... + (buf[i*8 + 7] << 40)
 *
 * loc_out receives a full uint64_t but only 48 bits are meaningful.
 *
 *
 * Usage example
 * -------------
 *     // Suppose 'p' is the pointer to the packed entry array (after the
 *     // chromosome-name block inside the decompressed fmt7 buffer).
 *
 *     uint16_t chr_id;
 *     uint64_t loc;
 *
 *     // decode entry #123
 *     FMT7_GET_LOC(p, 123, &chr_id, &loc);
 *     printf("chr_id=%u, loc=%" PRIu64 "\n", chr_id, loc);
 *
 *     // modify and re-encode it
 *     chr_id = 2;
 *     loc    = loc + 100;
 *     FMT7_SET_LOC(p, 123, chr_id, loc);
 *
 *
 * Notes
 * -----
 * • These macros operate safely on unaligned memory.
 * • They define the canonical binary layout for the decompressed fmt7 format.
 * • They must be used whenever reading or writing the uncompressed entries.
 * • Do not treat the entries as uint64_t*; only these macros define the
 *   correct byte order and bit layout.
 *
 */
#define FMT7_SET_LOC(buf, i, chr_id, loc)                   \
  do {                                                      \
    uint8_t *p = (buf) + ((i) * 8);                         \
    /* chr_id uses p[0], p[1] (16 bits little-endian) */    \
    uint16_t cid = (uint16_t)(chr_id);                      \
    p[0] = (uint8_t)(cid);                                  \
    p[1] = (uint8_t)(cid >> 8);                             \
                                                            \
    /* loc occupies p[2]..p[7] (48 bits little-endian) */   \
    uint64_t pos48 = (uint64_t)(loc) & 0xFFFFFFFFFFFFULL;   \
    for (int k = 0; k < 6; ++k) {                           \
      p[2 + k] = (uint8_t)(pos48 >> (8 * k));               \
    }                                                       \
  } while (0)

#define FMT7_GET_LOC(buf, i, chr_id_out, loc_out)           \
  do {                                                      \
    const uint8_t *p = (buf) + ((i) * 8);                   \
                                                            \
    /* read chr_id from p[0], p[1] */                       \
    uint16_t cid =  (uint16_t)p[0] | ((uint16_t)p[1] << 8); \
    (chr_id_out) = cid;                                    \
                                                            \
    /* read 48-bit coordinate from p[2..7] */               \
    uint64_t pos48 = 0;                                     \
    for (int k = 0; k < 6; ++k)                             \
      pos48 |= ((uint64_t)p[2 + k] << (8 * k));             \
                                                            \
    (loc_out) = pos48;                                     \
  } while (0)


/**
 * @brief Prepare access to an uncompressed/indexed fmt7 cdata_t.
 *
 * The uncompressed representation produced by fmt7_decompress() has layout:
 *
 *   [chr1\0chr2\0...chrN\0\0][8-byte entries...]
 *
 * This helper:
 *   - parses the leading chromosome-name block
 *   - allocates an array of pointers to those names
 *   - returns a pointer to the start of the packed entries

 * usage:
 *   uint8_t *locs_beg = NULL;
 *   char **chrms = NULL; int nchrms = 0;
 *   fmt7_prep(c, &locs_beg, &chrms, &nchrms);
 *
 * All chromosome name pointers in *chrms point directly into c->s and must
 * not be freed individually. Free only the *chrms array itself when done.
 *
 * @param c        Uncompressed fmt7 cdata_t (fmt=='7', compressed==0).
 * @param locs_beg Output: pointer to first packed 8-byte entry.
 * @param chrms    Output: allocated array of chromosome-name pointers.
 * @param nchrms   Output: number of chromosomes.
 */
static void fmt7_prep(const cdata_t c, uint8_t **locs_beg, char ***chrms, int *nchrms) {

  if (c.fmt != '7') wzfatal("[fmt7_prep] Expect format 7 but got %c.\n", c.fmt);
  if (c.compressed) wzfatal("[fmt7_prep] Expect uncompressed/indexed fmt7 (compressed=0).\n");

  /* ------------------------------------------------------------
   * First pass: count chromosome names.
   * Names are: "chr1\0chr2\0...\0chrN\0\0"
   * We count every '\0' until we hit the double-NUL sentinel.
   * ---------------------------------------------------------- */
  int count = 0;
  uint8_t *p = c.s;
  while (1) {
    if (p[0] == '\0') {         /* one chromosome name ended */
      count++;
      if (p[1] == '\0') break;  /* reached "\0\0" sentinel */
    }
    p++;
  }
  *nchrms = count;

  /* If there were no chromosomes, locs starts right after the sentinel,
     or there may be no entries at all. */
  if (count == 0) {
    *locs_beg = p + 2;  /* skip the double NUL, though there may be nothing after */
    *chrms = NULL;
    return;
  }

  /* ------------------------------------------------------------
   * Second pass: fill chrms[] with pointers to each name and
   * set *locs_beg to the start of the packed entries.
   * ---------------------------------------------------------- */
  char **names = calloc(count, sizeof(char *));
  if (!names) wzfatal("[fmt7_prep] Out of memory allocating chrms (%d names).\n", count);

  p = c.s;
  int idx = 0;
  while (p[0] != '\0') {
    names[idx++] = (char *)p;              /* name starts here */
    p += strlen((char *)p) + 1;            /* skip to next name */
  }

  /* p currently points at the second '\0' of the "\0\0" sentinel. */
  *locs_beg = p+1;
  *chrms = names;
}

/**
 * fmt7_decompress
 *
 * Re-encodes a fmt7 cdata_t (row coordinates) into a more index-friendly
 * layout stored directly in cdata_t::s:
 *
 *   [chr1\0chr2\0...chrN\0\0][packed entries...]
 *
 * where:
 *   - The first block is a sequence of chromosome names as C-strings,
 *     terminated by an extra '\0' (so "\0\0" marks the end of the list).
 *   - The second block is an array of fixed 8-byte records, one per CpG:
 *
 *       byte 0..1 : uint16_t chr_id   (index into the name list)
 *       byte 2..7 : 48-bit coordinate (beg1, little-endian)
 *
 * This replaces the old "just memcpy the compressed blob" behavior with
 * a more structured format that can be decoded with FMT7_GET_LOC().
 */
cdata_t fmt7_decompress(const cdata_t c) {

  cdata_t out = (cdata_t){0};

  /* number of CpGs in the compressed encoding */
  uint64_t n = fmt7_data_length((cdata_t *)&c);

  /* allocate coordinate buffer: n records * 8 bytes each */
  uint8_t *chrmlocs = calloc(n, 8 * sizeof(uint8_t));
  if (!chrmlocs) {
    fprintf(stderr, "[fmt7_decompress] Out of memory allocating %"PRIu64" entries\n", n);
    exit(1);
  }

  row_reader_t rdr = {0};
  char *chrm = NULL;
  uint64_t i = 0;

  char **chrms = NULL;
  int nchrms = 0;

  /* ------------------------------------------------------------------
   * Pass 1: iterate original fmt7 data, collect chromosomes and pack
   *         (chr_id, coordinate) into chrmlocs[].
   * ------------------------------------------------------------------ */
  while (row_reader_next_loc(&rdr, (cdata_t *)&c)) {
    if (chrm != rdr.chrm) {   /* new chromosome encountered */
      chrm = rdr.chrm;

      char **tmp = realloc(chrms, (nchrms + 1) * sizeof(char *));
      if (!tmp) {
        fprintf(stderr, "[fmt7_decompress] Out of memory realloc chrms\n");
        exit(1);
      }
      chrms = tmp;

      chrms[nchrms] = strdup(chrm);  /* own a copy of the name */
      if (!chrms[nchrms]) {
        fprintf(stderr, "[fmt7_decompress] Out of memory strdup(chrm)\n");
        exit(1);
      }
      nchrms++;
    }

    if (i >= n) {
      wzfatal("[fmt7_decompress] Internal error: more rows than fmt7_data_length (i=%"PRIu64", n=%"PRIu64")\n", i, n);
    }

    uint64_t chr_id = (uint64_t)(nchrms - 1);  /* index of current chromosome */
    FMT7_SET_LOC(chrmlocs, i, chr_id, rdr.value);
    i++;
  }

  if (i != n) {
    wzfatal("[fmt7_decompress] Warning: collected %"PRIu64" rows, fmt7_data_length reported %"PRIu64"\n", i, n);
  }

  /* ------------------------------------------------------------------
   * Build the output buffer layout:
   *   [chr1\0chr2\0...chrN\0\0][n * 8 bytes of packed entries]
   * ------------------------------------------------------------------ */

  /* 1) compute size of chromosome name section */
  uint64_t names_bytes = 0;
  for (int j = 0; j < nchrms; ++j) {
    names_bytes += (uint64_t)strlen(chrms[j]) + 1;  /* +1 for '\0' */
  }
  names_bytes += 1;  /* extra '\0' as sentinel: double NUL terminator */

  uint64_t total_bytes = names_bytes + n * 8;

  out.s = malloc(total_bytes);
  if (!out.s) {
    fprintf(stderr, "[fmt7_decompress] Out of memory allocating %"PRIu64" bytes\n", total_bytes);
    exit(1);
  }

  uint8_t *p = out.s;

  /* 2) write chromosome names as consecutive C-strings */
  for (int j = 0; j < nchrms; ++j) {
    size_t len = strlen(chrms[j]);
    memcpy(p, chrms[j], len);
    p[len] = '\0';
    p += len + 1;
  }

  /* extra NUL sentry: two consecutive '\0' means “end of name list” */
  *p++ = '\0';

  /* 3) copy packed (chr_id, coord) entries */
  memcpy(p, chrmlocs, n * 8);

  /* ------------------------------------------------------------------
   * Fill cdata_t metadata and clean temp buffers
   * ------------------------------------------------------------------ */

  out.fmt        = '7';
  out.n          = n;           /* total bytes in out.s */
  out.compressed = 0;           /* now in the “decompressed/indexed” form */
  out.unit       = 8;           /* logical data unit: 8-byte packed entries */
  out.aux        = NULL;

  /* free temporary buffers */
  free(chrmlocs);
  for (int j = 0; j < nchrms; ++j) free(chrms[j]);
  free(chrms);

  return out;
}


// beg and end are both 0-based
cdata_t fmt7_sliceToBlock(cdata_t *cr, uint64_t beg, uint64_t end) {
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
  cdata_t cr2 = {0};
  while (row_reader_next_loc(&rdr, cr)) {
    if (i >= beg && i <= end) {
      if (chrm != rdr.chrm) {
        if (chrm) fmt7c_append_end(&(cr2.s), &n);
        chrm = rdr.chrm;
        fmt7c_append_chrm(chrm, &(cr2.s), &n);
        last = 0;
      }
      fmt7c_append_loc(rdr.value - last, &(cr2.s), &n);
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
  cr2.unit = cr->unit;
  cr2.fmt = cr->fmt;
  cr2.n = n;
  cr2.compressed = 1;
  return cr2;
}

cdata_t fmt7_sliceToIndices(cdata_t *c, int64_t *row_indices, int64_t n_indices) {

  cdata_t inflated = fmt7_decompress(*c);
  uint8_t *locs_beg = NULL;      // point to inflated.s
  char **chrms = NULL; int nchrms = 0;
  fmt7_prep(inflated, &locs_beg, &chrms, &nchrms);

  char *chrm = NULL; uint64_t last = 0;
  uint8_t *s_out = NULL; uint64_t n_out = 0;
  for (uint64_t i=0; i<(uint64_t) n_indices; ++i) {
    uint16_t ichrm = 0; uint64_t loc = 0;
    FMT7_GET_LOC(locs_beg, row_indices[i]-1, ichrm, loc);
    if (!chrm || chrm != chrms[ichrm] || loc < last) {
      if (chrm) fmt7c_append_end(&s_out, &n_out);
      fmt7c_append_chrm(chrms[ichrm], &s_out, &n_out);
      chrm = chrms[ichrm];
      last = 0;
    }
    fmt7c_append_loc(loc - last, &s_out, &n_out);
    last = loc;
  }
  free(chrms);

  cdata_t c_out = {0};
  c_out.s = s_out;
  c_out.n = n_out;
  c_out.fmt = '7';
  c_out.compressed = 1;
  c_out.unit = 1;
  return c_out;
}

cdata_t fmt7_sliceToMask(cdata_t *cr, cdata_t *c_mask) {

  row_reader_t rdr = {0};
  uint64_t n = 0;
  uint64_t i = 0;
  char *chrm = NULL; uint64_t last = 0;
  cdata_t cr2 = {0};
  while (row_reader_next_loc(&rdr, cr)) {
    if (c_mask->s[i>>3]&(1<<(i&0x7))) {
      if (chrm != rdr.chrm) {
        if (chrm) fmt7c_append_end(&(cr2.s), &n);
        chrm = rdr.chrm;
        fmt7c_append_chrm(chrm, &(cr2.s), &n);
        last = 0;
      }
      fmt7c_append_loc(rdr.value - last, &(cr2.s), &n);
      last = rdr.value;
    }
    i++;
  }
  cr2.unit = cr->unit;
  cr2.fmt = cr->fmt;
  cr2.n = n;
  cr2.compressed = 1;
  return cr2;
}

stats_t* summarize1_queryfmt7(
  cdata_t *c, cdata_t *c_mask, uint64_t *n_st, char *sm, char *sq, config_t *config) {

  stats_t *st = NULL;
  uint64_t n = c->n;
  
  /* -------------------------------------------------------------
   * No mask: just per-chromosome counts.
   * ----------------------------------------------------------- */
  if (c_mask->n == 0) {

    uint8_t *locs_beg = NULL;
    char **chrms = NULL; int nchrms = 0;
    fmt7_prep(*c, &locs_beg, &chrms, &nchrms);

    *n_st = (uint64_t) nchrms;
    st = calloc(*n_st, sizeof(stats_t));

    uint64_t* chrm_cnts = calloc(nchrms, sizeof(uint64_t));
    for (uint64_t i=0; i<c->n; ++i) {
      uint16_t ichr = 0; uint64_t loc = 0;
      FMT7_GET_LOC(locs_beg, i, ichr, loc);
      chrm_cnts[ichr]++;
    }

    for (int ichr = 0; ichr < nchrms; ++ichr) {
      stats_t *s = &st[ichr];
      s->n_u   = n;                    /* universe = all CpGs */
      s->n_q   = chrm_cnts[ichr];      /* CpGs on this chromosome */
      s->n_m   = 0;
      s->n_o   = 0;
      s->beta  = -1.0;                 /* force Beta to NA; sum_depth stays 0 */

      s->sm = strdup(sm);
      if (config->section_name) {
        kstring_t tmp = {0};
        ksprintf(&tmp, "%s-%s", sq, chrms[ichr] ? chrms[ichr] : "");
        s->sq = tmp.s;
      } else {
        s->sq = strdup(chrms[ichr] ? chrms[ichr] : "");
      }
    }
    free(chrms);

  } else {  /* other mask formats unsupported for fmt7 */

    fprintf(stderr, "[%s:%d] Mask format %c unsupported for query format 7.\n",
            __func__, __LINE__, c_mask->fmt);
    fflush(stderr);
    exit(1);
  }

  return st;
}
