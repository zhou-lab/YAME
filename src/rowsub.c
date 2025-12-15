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

#include <sys/stat.h>
#include <sys/types.h>
#include <string.h>
#include "cfile.h"
#include "cdata.h"

/**
 * yame rowsub
 * ===========
 *
 * Goal
 * ----
 * Extract a subset of rows from each dataset (cdata record) in an input .cx stream,
 * writing the resulting .cx stream to stdout.
 *
 * "Rows" correspond to vector positions in each cdata record:
 *   - For fmt 0/1/6: rows are bit/2-bit positions
 *   - For fmt 3/4:   rows are fixed-width units (unit bytes per row)
 *   - For fmt 2:     rows are state entries, but fmt2 has a key section that must be preserved
 *   - For fmt 7:     rows are coordinate entries; slicing uses fmt7_sliceTo* helpers
 *
 * Selection modes
 * ---------------
 * rowsub supports several ways to specify which rows to keep:
 *
 * (A) Explicit row indices list (-l)
 *   - Input file: one integer [index1] per line (1-based).
 *   - No sorting required; order is preserved.
 *   - Applied via sliceToIndices() (or fmt7_sliceToIndices for fmt7).
 *
 * (B) Coordinate list resolved through a row coordinate table (-R + -L)
 *   - -R supplies a format-7 "row coordinate dataset" that maps rows to coordinates.
 *   - -L supplies a list of coordinates, one per line: "chrm_beg1".
 *   - We resolve each coordinate to a row index using row_finder_search(), then slice.
 *   - -1 optionally prepends the subsetted coordinate dataset as the first output record.
 *
 * (C) Binary mask (-m)
 *   - Mask must be fmt 0/1; it is converted to fmt0 bitset.
 *   - Keep rows where mask bit is 1.
 *   - For fmt2, output preserves the key section and filters only the data section
 *     (see sliceToMask()).
 *
 * (D/E) Contiguous block (-B or -I)
 *   - -B: keep a range by absolute indices: beg0 and optional end1 (exclusive, 1-based).
 *   - -I: keep a block by blockIndex0 and optional blockSize; computes beg0/end0.
 *   - If neither is provided, defaults are config.beg/config.end.
 *
 * Precedence
 * ----------
 * If multiple selection mechanisms are provided, rowsub behaves as:
 *   explicit indices (-l / -L) > mask (-m) > block (-B / -I / default).
 * This is implemented by checking row_indices first, then c_mask, then block slicing.
 *
 * Format-specific slicing
 * -----------------------
 * - Non-fmt7 data are decompressed first (cdata_t c2 = decompress(c)), sliced in memory,
 *   then re-compressed before writing.
 * - fmt2 requires preserving the key section when slicing blocks or masks:
 *     [keys...][\\0][filtered data rows...]
 * - fmt7 records are handled via fmt7_sliceToIndices / fmt7_sliceToMask / fmt7_sliceToBlock.
 *
 * Output
 * ------
 * Writes a valid .cx stream to stdout. Caller can redirect to a file.
 * No index is written by this subcommand.
 */

typedef struct config_t {
  char *fname_rindex;
  uint64_t beg;
  uint64_t end;
  int64_t index;
  int64_t isize;
} config_t;

// rowsub supports multiple mutually-exclusive row selection modes.
// Precedence in code is: -l/-L (explicit indices) > -m (mask) > -B/-I (block) > default block.
static int usage(config_t *config) {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "  yame rowsub [options] <in.cx> > out.cx\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Purpose:\n");
  fprintf(stderr, "  Subset (slice) rows from each dataset (record) in a CX stream.\n");
  fprintf(stderr, "  Output is always written to stdout.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Row selection modes (choose one):\n");
  fprintf(stderr, "  (A) Explicit row indices (1-based list):\n");
  fprintf(stderr, "      -l <idx.txt>     One [index1] per line (1-based). Order preserved; no sorting required.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  (B) Explicit genomic coordinates via row coordinate table (format 7):\n");
  fprintf(stderr, "      -R <rows.cx>     Row coordinate dataset (format 7; e.g. BED-like coordinates).\n");
  fprintf(stderr, "      -L <coord.txt>   One [chrm]_[beg1] per line (1-based beg). Requires -R.\n");
  fprintf(stderr, "                       Order preserved; no sorting required.\n");
  fprintf(stderr, "      -1               If -R is provided, emit the subsetted row coordinates as the FIRST dataset.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  (C) Mask-based filtering (binary mask):\n");
  fprintf(stderr, "      -m <mask.cx>     Mask file (format 0/1 only). Rows with bit=1 are kept.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  (D) Contiguous block by absolute row range (0-based):\n");
  fprintf(stderr, "      -B <beg0>[_<end1>]\n");
  fprintf(stderr, "         Keep rows in [beg0, end0] where end0 = end1-1.\n");
  fprintf(stderr, "         If <end1> is omitted, keep a single row at beg0.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  (E) Contiguous block by block index and size (0-based):\n");
  fprintf(stderr, "      -I <blockIndex0>[_<blockSize>]\n");
  fprintf(stderr, "         Keep rows:\n");
  fprintf(stderr, "           beg0 = blockIndex0 * blockSize\n");
  fprintf(stderr, "           end0 = (blockIndex0+1)*blockSize - 1\n");
  fprintf(stderr, "         If <blockSize> is omitted, default blockSize=%"PRIu64".\n", (uint64_t)config->isize);
  fprintf(stderr, "\n");
  fprintf(stderr, "Other options:\n");
  fprintf(stderr, "  -h               Show this help message.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Index conventions:\n");
  fprintf(stderr, "  - '0' suffix means 0-based (beg0, blockIndex0).\n");
  fprintf(stderr, "  - '1' suffix means 1-based (index1, beg1, end1).\n");
  fprintf(stderr, "  - For -B, end is provided as end1 (exclusive, 1-based), internally converted to end0.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Notes:\n");
  fprintf(stderr, "  * For format 2 (state data), the key section is preserved when slicing.\n");
  fprintf(stderr, "  * Format 7 (row coordinates) is sliced with fmt7_* helpers.\n");
  fprintf(stderr, "  * If multiple selection options are given, the effective precedence is:\n");
  fprintf(stderr, "      -l/-L  >  -m  >  -B/-I  >  default.\n");
  fprintf(stderr, "\n");
  return 1;
}

static int64_t *load_row_indices(char *fname, int64_t *n) {

  int64_t *indices = NULL;
  /* snames_t snames = {0}; */
  if (fname == NULL) { *n = 0; return indices; }
  gzFile fp;
  if (strcmp(fname, "-") == 0) {
    fp = gzdopen(fileno(stdin), "r");
  } else {
    fp = gzopen(fname, "r");
    if (!fp) {
      fprintf(stderr, "[%s:%d] Fatal, cannot open file: %s\n",
              __func__, __LINE__, fname);
      fflush(stderr);
      exit(1);
    }
  }
  
  if (fp == NULL) return indices;

  char *line = NULL; *n=0;
  while (gzFile_read_line(fp, &line) > 0) {
    char *field = NULL;
    if (line_get_field(line, 0, "\t", &field)) {
      indices = realloc(indices, sizeof(int64_t)*((*n)+1));
      if (indices == NULL) {
        fprintf(stderr, "Failed to allocate memory\n");
        fflush(stderr);
        exit(1);
      }
      indices[(*n)++] = strtoll(field, NULL, 10);
      free(field);
    }
  }
  free(line);
  gzclose(fp);
  return indices;
}

static int split_string_and_number(const char* input, char** out_string, uint64_t* out_number) {
  // Find the position of the underscore ('_')
  const char *underscore = strchr(input, '_');
  if (!underscore) return -1; // underscore not found

  // Allocate memory for the string part (including the null-terminator)
  *out_string = malloc(underscore - input + 1);
  if (!*out_string) return -2; // memory allocation failed

  // Copy the string part to the output
  strncpy(*out_string, input, underscore - input);
  (*out_string)[underscore - input] = '\0'; // Null-terminate

  // Convert the number part to uint64_t
  char *endptr;
  *out_number = strtoull(underscore + 1, &endptr, 10);

  // Check if conversion was successful
  if (*endptr != '\0') { free(*out_string); return -3; }
  return 0; // success
}

static int64_t *load_row_indices_by_names(char *fname_rnindex, cdata_t *cr, int64_t *n_indices) {
  int64_t *indices = NULL;
  /* snames_t snames = {0}; */
  if (fname_rnindex == NULL) { *n_indices = 0; return indices; }
  gzFile fp;
  if (strcmp(fname_rnindex, "-") == 0) {
    fp = gzdopen(fileno(stdin), "r");
  } else {
    fp = gzopen(fname_rnindex, "r");
    if (!fp) {
      fprintf(stderr, "[%s:%d] Fatal, cannot open file: %s\n",
              __func__, __LINE__, fname_rnindex);
      fflush(stderr);
      exit(1);
    }
  }
  
  if (fp == NULL) return indices;
  
  row_finder_t fdr = init_finder(cr);
  *n_indices = 0;
  char *line = NULL;
  while (gzFile_read_line(fp, &line) > 0) {
    char *chrm = NULL;
    uint64_t beg1 = 0;
    if (split_string_and_number(line, &chrm, &beg1) < 0) {
      fprintf(stderr, "Failed to extract coordinate: %s\n", line);
      fflush(stderr);
      exit(1);
    }
    indices = realloc(indices, ((*n_indices)+1)*sizeof(int64_t));
    indices[(*n_indices)] = row_finder_search(chrm, beg1, &fdr, cr);
    if (!indices[(*n_indices)]) {
      fprintf(stderr, "[%s:%d] Cannot find coordinate: %s\n", __func__, __LINE__, line);
      fflush(stderr);
      exit(1);
    }
    (*n_indices)++;
    free(chrm);
  }
  free_row_finder(&fdr);
  free(line);
  gzclose(fp);
  return indices;
}

static cdata_t sliceToIndices(cdata_t *c, int64_t *row_indices, int64_t n) {
  if (c->compressed) {
    fprintf(stderr, "[%s:%d] Slicing compressed data.\n", __func__, __LINE__);
    fflush(stderr);
    exit(1);
  }
  cdata_t c2 = {0};
  c2.unit = c->unit;
  c2.s = realloc(c2.s, n*c2.unit);
  c2.fmt = c->fmt;
  int64_t i;
  for (i=0; i<n; ++i) {
    memcpy(c2.s+c2.unit*i, c->s+c->unit*(row_indices[i]-1), c->unit);
  }
  c2.n = n;
  c2.compressed = 0;
  return c2;
}

/**
 * sliceToBlock()
 * ----------------
 * Extract a contiguous block of uncompressed rows from a cdata_t.
 *
 * Given an uncompressed vector `c` (compressed == 0), return a new cdata_t
 * containing rows [beg, end] (0-based, inclusive).  The output preserves:
 *     • fmt     (same format as input)
 *     • unit    (same per-row byte width)
 *     • compressed = 0
 *
 * Behavior:
 *   • For most formats (unit > 0), this is a simple memcpy of the slice:
 *
 *         out.s = c->s + unit*beg  →  unit*(end-beg+1) bytes
 *
 *   • For format 2 (state data), the key section must be preserved.
 *     The output buffer is:
 *
 *         [key strings ... '\0'][data slice...]
 *
 *     where the key region length is obtained via fmt2_get_keys_nbytes().
 *
 * Arguments:
 *   c    : input cdata_t (must be uncompressed)
 *   beg  : first row index to keep (0-based)
 *   end  : last row index to keep (0-based; clipped to n-1)
 *
 * Returns:
 *   A new cdata_t containing just the requested row block.
 *
 * Notes:
 *   • Caller must free the returned cdata_t.s.
 *   • beg and end refer to logical row indices, not byte offsets.
 */
static cdata_t sliceToBlock(cdata_t *c, uint64_t beg, uint64_t end) {
  assert(!c->compressed);
  if (end > c->n-1) end = c->n-1; // 0-base
  if (beg > c->n-1)
    wzfatal("[%s:%d] Begin (%"PRIu64") is bigger than the data vector size (%"PRIu64").\n", __func__, __LINE__, beg, c->n);

  cdata_t c_out = {0};
  c_out.unit = c->unit;
  c_out.fmt = c->fmt;
  if (c_out.fmt == '2') {
    uint64_t keys_nb = fmt2_get_keys_nbytes(c);
    c_out.s = calloc(1, (end-beg+1)*c_out.unit + keys_nb + 1);
    memcpy(c_out.s, c->s, keys_nb + 1);
    memcpy(c_out.s+keys_nb+1, c->s+keys_nb+1+c->unit*beg, c->unit*(end-beg+1));
    c_out.n = end-beg+1;
    // TODO: format 7 should be merged here, I have a separate fmt7_sliceToBlock
  } else {
    c_out.s = realloc(c_out.s, (end-beg+1)*c_out.unit);
    memcpy(c_out.s, c->s+c->unit*beg, c->unit*(end-beg+1));
    c_out.n = end-beg+1;
  }
  c_out.compressed = 0;
  
  return c_out;
}

/* static cdata_t sliceToMask(cdata_t *c, cdata_t *c_mask) { */
/*   assert(!c->compressed); */
/*   if (c->n != c_mask->n) */
/*     wzfatal("[%s:%d] Mask (N=%"PRIu64") and data (N=%"PRIu64") are of different lengths.\n", __func__, __LINE__, c_mask->n, c->n); */

/*   uint64_t n = 0; */
/*   for (uint64_t i=0; i<c->n; ++i) */
/*     if (FMT0_IN_SET(*c_mask, i)) n++; */

/*   cdata_t c_out = {0}; */
/*   c_out.unit = c->unit; */
/*   c_out.s = realloc(c_out.s, n*c_out.unit); */
/*   c_out.fmt = c->fmt; */
/*   for (uint64_t i=0, k=0; i<c->n; ++i) { */
/*     if (FMT0_IN_SET(*c_mask, i)) */
/*       memcpy(c_out.s+(k++)*c->unit, c->s+i*c->unit, c->unit); */
/*   } */
/*   c_out.n = n; */
/*   c_out.compressed = 0; */
/*   return c_out; */

/* } */

static cdata_t sliceToMask(cdata_t *c, cdata_t *c_mask) {
  assert(!c->compressed);
  if (c->n != c_mask->n)
    wzfatal("[%s:%d] Mask (N=%"PRIu64") and data (N=%"PRIu64") are of different lengths.\n", __func__, __LINE__, c_mask->n, c->n);

  /* Count how many rows will be kept */
  uint64_t n = 0;
  for (uint64_t i = 0; i < c->n; ++i)
    if (FMT0_IN_SET(*c_mask, i)) n++;

  cdata_t c_out = (cdata_t){0};
  c_out.unit = c->unit;
  c_out.fmt  = c->fmt;
  c_out.compressed = 0;
  c_out.n = n;  /* number of rows (states), not including key bytes */
  if (c_out.fmt == '2') {
    /* layout: [keys...][\0][filtered data rows...] */
    uint64_t keys_nb = fmt2_get_keys_nbytes(c); // no trailing '\0'
    c_out.s = calloc(1, keys_nb + 1 + n * c_out.unit);
    memcpy(c_out.s, c->s, keys_nb + 1); // copy key section + '\0'
    uint8_t *dst      = c_out.s + keys_nb + 1;
    uint8_t *src_data = fmt2_get_data(c);  /* start of original data section */
    for (uint64_t i = 0; i < c->n; ++i) {
      if (FMT0_IN_SET(*c_mask, i)) {
        memcpy(dst, src_data + i * c->unit, c->unit);
        dst += c->unit;
      }
    }
  } else { // all other formats
    if (n > 0) {
      c_out.s = malloc(n * c_out.unit);
      for (uint64_t i = 0, k = 0; i < c->n; ++i)
        if (FMT0_IN_SET(*c_mask, i))
          memcpy(c_out.s + (k++) * c->unit, c->s + i * c->unit, c->unit);
    }
  }
  
  return c_out;
}


int main_rowsub(int argc, char *argv[]) {

  config_t config = {
    .fname_rindex = NULL,
    .index = -1, .isize = 1000000, .beg = 0, .end = 1};
  int c; char *fname_row = NULL; char *fname_mask = NULL;
  char *fname_rnindex = NULL; int add_row_coordinates = 0;
  char *B_option = NULL, *I_option = NULL;
  while ((c = getopt(argc, argv, "1R:m:l:L:B:I:h"))>=0) {
    switch (c) {
    case '1': add_row_coordinates = 1; break;
    case 'R': fname_row = strdup(optarg); break;
    case 'm': fname_mask = strdup(optarg); break;
    case 'l': config.fname_rindex = strdup(optarg); break;
    case 'L': fname_rnindex = strdup(optarg); break;
    case 'B': B_option = strdup(optarg); break;
    case 'I': I_option = strdup(optarg); break;
    case 'h': return usage(&config); break;
    default: usage(&config); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (B_option) {               // [rowIndexBeg0]_(rowIndexEnd1)
    char *token = strtok(B_option, "_");
    if (token != NULL) {
      config.beg = atoi(token);
      token = strtok(NULL, "_");
      if (token != NULL) {
        config.end = atoi(token) - 1; // convert to 0-base
      } else {
        config.end = config.beg;
      }
    }
  }

  if (I_option) {               // [blockIndex0]_(blockSize)
    char *token = strtok(I_option, "_");
    if (token != NULL) {
      int64_t blockIndex0 = atoi(token);
      int64_t blockSize = config.isize;
      token = strtok(NULL, "_");
      if (token != NULL) {
        blockSize = atoi(token);
      }
      config.beg = blockIndex0 * blockSize;
      config.end = (blockIndex0+1) * blockSize - 1;
    }
  }
  
  if (optind + 1 > argc) {
    usage(&config); 
    wzfatal("Please supply input files.\n");
  }
  char *fname = argv[optind];

  int64_t *row_indices = NULL;
  int64_t n_indices=0;
  if (config.fname_rindex) {
    row_indices = load_row_indices(config.fname_rindex, &n_indices);
  }

  cdata_t c_mask = {0};
  if (fname_mask) {
    cfile_t cf_mask = open_cfile(fname_mask);
    c_mask = read_cdata1(&cf_mask);
    if (c_mask.fmt >= '2') {
      fprintf(stderr, "[%s:%d] Mask is not binary.\n", __func__, __LINE__);
      fflush(stderr);
      exit(1);
    }
    convertToFmt0(&c_mask);
    bgzf_close(cf_mask.fh);
  }

  cfile_t cf = open_cfile(fname);
  BGZF *fp_out = bgzf_dopen(fileno(stdout), "w");
  if (fp_out == NULL) {
    fprintf(stderr, "[%s:%d] Cannot open output stream.\n", __func__, __LINE__);
    fflush(stderr);
    exit(1);
  }

  if (fname_rnindex && !fname_row) {
    fprintf(stderr, "[%s:%d] Missing -R for BED coordinates.\n", __func__, __LINE__);
    fflush(stderr);
    exit(1);
  }

  if (fname_row) {
    cfile_t cf_row = open_cfile(fname_row);
    cdata_t cr = read_cdata1(&cf_row);
    if (!row_indices && fname_rnindex) {
      row_indices = load_row_indices_by_names(fname_rnindex, &cr, &n_indices);
    }
    if (add_row_coordinates) {
      cdata_t cr2;
      if (row_indices) cr2 = fmt7_sliceToIndices(&cr, row_indices, n_indices);
      else if (c_mask.n) cr2 = fmt7_sliceToMask(&cr, &c_mask);
      else cr2 = fmt7_sliceToBlock(&cr, config.beg, config.end);
      cdata_write1(fp_out, &cr2);
      free_cdata(&cr2);
    }
    free_cdata(&cr);
    bgzf_close(cf_row.fh);
  }
  
  while (1) {
    cdata_t c = read_cdata1(&cf);
    if (c.n == 0) break;
    if (c.fmt == '7') {
      cdata_t c2;
      if (row_indices) c2 = fmt7_sliceToIndices(&c, row_indices, n_indices);
      else if (c_mask.n) c2 = fmt7_sliceToMask(&c, &c_mask);
      else c2 = fmt7_sliceToBlock(&c, config.beg, config.end);
      cdata_write1(fp_out, &c2);
      free_cdata(&c2);
    } else {
      cdata_t c2 = decompress(c);
      cdata_t c3 = {0};
      if (row_indices) c3 = sliceToIndices(&c2, row_indices, n_indices);
      else if (c_mask.n) c3 = sliceToMask(&c2, &c_mask);
      else c3 = sliceToBlock(&c2, config.beg, config.end);
      cdata_compress(&c3);
      cdata_write1(fp_out, &c3);
      free(c3.s); free(c2.s);
    }
    free_cdata(&c);
  }
  bgzf_close(cf.fh);
  bgzf_close(fp_out);

  if (n_indices) free(row_indices);
  if (fname_row) free(fname_row);
  if (fname_rnindex) free(fname_rnindex);
  
  return 0;
}
