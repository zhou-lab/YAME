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

#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdio.h>
#include "cfile.h"
#include "vector.h"
#include "snames.h"

/**
 * yame unpack
 * ===========
 *
 * Goal
 * ----
 * Convert selected CX records into a tab-delimited table on stdout:
 *   - each output line corresponds to one row position i
 *   - each output column corresponds to one selected cdata record (sample)
 * Optionally prepend a coordinate column derived from a separate row-coordinate dataset (-R),
 * or print coordinates directly when the selected dataset is format 7.
 *
 * Sample selection
 * ----------------
 * The tool loads a vector of records (cdata_v *cs) using the following priority:
 *
 *   1) If an index exists and snames.n > 0:
 *        cs = read_cdata_with_snames(&cf, idx, &snames)
 *   2) Else if -a:
 *        cs = read_cdata_all(&cf)
 *   3) Else if -H N:
 *        cs = read_cdata_from_head(&cf, N)
 *   4) Else if -T N (requires index):
 *        cs = read_cdata_from_tail(&cf, idx, N)
 *   5) Else:
 *        cs = read_cdata_from_head(&cf, 1)
 *
 * Random access requirement
 * -------------------------
 * Selecting by explicit sample names or selecting from tail requires random access,
 * which in turn requires an input index (.cxi), unless the input is stdin.
 *
 * Coordinates / left column (-R / format 7)
 * -----------------------------------------
 * - If -R is provided, a separate CX record is read as cr and printed as the first
 *   column(s) for every row i (print_cdata1(&cr, i, pfmt)).
 * - If the first selected dataset is format 7, unpack treats column 1 as coordinates
 *   (col1_is_row_index) and prints coordinates from that dataset.
 * - Coordinate formatting is controlled by pfmt.ref (-r):
 *     0: chrm beg0 end1
 *     1: chrm beg0 end0
 *     else: chrm_beg1
 *
 * Value printing (print_cdata1)
 * -----------------------------
 * Printing is format-specific:
 * - fmt0: bit (0/1)
 * - fmt1: raw byte/ASCII value
 * - fmt2: state label string (f2_get_string)
 * - fmt3: controlled by pfmt.data (-f):
 *     0 : print packed MU uint64
 *    <0 : print "M<TAB>U"
 *    >0 : print beta; "NA" if cov < pfmt.data or cov==0
 * - fmt4: float; negative values printed as NA
 * - fmt5: ternary; value 2 printed as NA
 * - fmt6: ALSO controlled by pfmt.data (-f):
 *    <0 : print "value<TAB>universe" (e.g., 1<tab>1 / 0<tab>1 / NA<tab>0)
 *     0 : print 0/1, NA coded as '2'
 *    >0 : print raw 2-bit code (FMT6_2BIT)
 * - fmt7: prints coordinates via row_reader_t (fmt7_next_bed), not a scalar value.
 *
 * Chunked printing (-c / -s)
 * --------------------------
 * In chunk mode, unpack repeatedly decompresses each selected record, slices a contiguous
 * row block of size s, and prints it. This reduces peak memory compared to inflating the
 * full vector. Format 7 chunking is not supported (explicitly rejected).
 *
 * Header printing (-C)
 * --------------------
 * If -C is requested, unpack prints a header line. When sample names were not explicitly
 * provided, it derives names from the input index (requires .cxi). If a coordinate column
 * is printed (either via -R or fmt7), the header includes the coordinate field names.
 *
 * Unit override (-u)
 * ------------------
 * Sets c->unit for each selected record before decompression (0=auto-infer).
 * Allowed values are {0,1,2,4,6,8}; smaller units reduce memory but may be lossy.
 */

static int usage(void) {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "  yame unpack [options] <in.cx> [sample1 sample2 ...]\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Purpose:\n");
  fprintf(stderr, "  Print selected records from a .cx file as a tab-delimited table.\n");
  fprintf(stderr, "  Each output row is a genomic row index; each output column is a selected sample/record.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Sample selection (default: first record):\n");
  fprintf(stderr, "  -a            Output all records in the file.\n");
  fprintf(stderr, "  -l <list>     Sample list file (one name per line).\n");
  fprintf(stderr, "                Ignored if sample names are provided as trailing arguments.\n");
  fprintf(stderr, "  -H <N>        Output the first N samples.\n");
  fprintf(stderr, "  -T <N>        Output the last  N samples (requires index).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Row coordinates (optional first column):\n");
  fprintf(stderr, "  -R <rows.cx>  Row coordinate dataset (CX; typically format 7).\n");
  fprintf(stderr, "  -r <mode>     Coordinate print mode (default: 0):\n");
  fprintf(stderr, "                0: chrm<tab>beg0<tab>end1   (cg-style)\n");
  fprintf(stderr, "                1: chrm<tab>beg0<tab>end0   (allc-style)\n");
  fprintf(stderr, "                else: chrm_beg1\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Output formatting:\n");
  fprintf(stderr, "  -C            Print a header line (column names).\n");
  fprintf(stderr, "  -u <bytes>    Inflated unit-size override (0=auto; allowed: 1,2,4,6,8).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Value printing (-f):\n");
  fprintf(stderr, "  -f <N>        Print mode for certain formats (default: 0):\n");
  fprintf(stderr, "                For format 3 (MU):\n");
  fprintf(stderr, "                  N == 0 : print packed MU (uint64)\n");
  fprintf(stderr, "                  N  < 0 : print M<tab>U (two columns)\n");
  fprintf(stderr, "                  N  > 0 : print beta; print NA if cov < N or cov==0\n");
  fprintf(stderr, "                For format 6 (set+universe):\n");
  fprintf(stderr, "                  N == 0 : print 0/1, NA coded as '2'\n");
  fprintf(stderr, "                  N  < 0 : print value<tab>universe  (e.g., 1<tab>1, 0<tab>1, NA<tab>0)\n");
  fprintf(stderr, "                  N  > 0 : print raw 2-bit code (FMT6_2BIT)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Chunked printing:\n");
  fprintf(stderr, "  -c            Enable chunked printing (reduces peak memory).\n");
  fprintf(stderr, "  -s <rows>     Chunk size in rows (default: 1000000).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Other:\n");
  fprintf(stderr, "  -h            Show this help message.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Notes:\n");
  fprintf(stderr, "  * Selecting by sample name or using -T requires an index (.cxi) unless reading from stdin.\n");
  fprintf(stderr, "  * Chunking does not support format 7 datasets.\n");
  fprintf(stderr, "\n");
  return 1;
}

typedef struct cdata_pfmt_t {
  int data;
  int ref;
} cdata_pfmt_t;

int fmt7_next_bed(cdata_t *c);
static void print_cdata1(cdata_t *c, uint64_t i, cdata_pfmt_t pfmt) {
  switch (c->fmt) {
  case '0': {
      fputc(((c->s[i>>3]>>(i&0x7))&0x1)+'0', stdout);
      break;
  }
  case '1': {
    fputc(c->s[i], stdout);
    break;
  }
  case '2': {
    fprintf(stdout, "%s", f2_get_string(c, i));
    break;
  }
  case '3': {
    uint64_t mu = f3_get_mu(c, i);
    if (pfmt.data == 0)
      fprintf(stdout, "%"PRIu64"", mu);
    else if (pfmt.data < 0)
      fprintf(stdout, "%"PRIu64"\t%"PRIu64"",mu>>32, mu<<32>>32);
    else {
      uint64_t M = mu>>32;
      uint64_t U = mu<<32>>32;
      if ((M==0 && U==0) || (M+U) < (uint64_t) pfmt.data) fputs("NA", stdout);
      else fprintf(stdout, "%1.3f", (double) M/(M+U));
    }
    break;
  }
  case '4': {
    if (((float_t*) (c->s))[i]<0) {
      fputs("NA", stdout);
    } else {
      fprintf(stdout, "%1.3f", ((float_t*) (c->s))[i]);
    }
    break;
  }
  case '5': {
    if (c->s[i] == 2) {
      fputs("NA", stdout);
    } else {
      fputc(c->s[i]+'0', stdout);
    }
    break;
  }
  case '6': {
    if (pfmt.data < 0) {
      if (FMT6_IN_UNI(*c, i)) {
        if (FMT6_IN_SET(*c, i)) {
          fputs("1\t1", stdout);
        } else {
          fputs("0\t1", stdout);
        }
      } else {
        fputs("NA\t0", stdout);
      }
    } else if (pfmt.data == 0) {
      if (FMT6_IN_UNI(*c, i)) {
        if (FMT6_IN_SET(*c, i)) {
          fputc('1', stdout);
        } else {
          fputc('0', stdout);
        }
      } else {
        fputc('2', stdout);
      }
    } else {
      fputc('0'+FMT6_2BIT(*c, i), stdout);
    }
    break;
  }
  case '7': {
    if (!fmt7_next_bed(c)) {
      fprintf(stderr, "[%s:%d] next BED record unfound.\n", __func__, __LINE__);
      fflush(stderr);
      exit(1);
    }
    row_reader_t *rdr = (row_reader_t*) c->aux;
    if (rdr->index != i+1) {
      fprintf(stderr, "[%s:%d] row reader index mismatch (i=%"PRIu64", rdr.index=%"PRIu64").\n",
              __func__, __LINE__, i, rdr->index);
      fflush(stderr);
      exit(1);
    }
    if (pfmt.ref == 0) {
      fprintf(stdout, "%s\t%"PRIu64"\t%"PRIu64"", rdr->chrm, rdr->value-1, rdr->value+1);
    } else if (pfmt.ref == 1) {
      fprintf(stdout, "%s\t%"PRIu64"\t%"PRIu64"", rdr->chrm, rdr->value-1, rdr->value);
    } else {
      fprintf(stdout, "%s_%"PRIu64"", rdr->chrm, rdr->value);
    }
    break;
  }
  default: usage(); wzfatal("Unrecognized format: %c.\n", c->fmt);
  }
}

static void print_cdata_chunk(cdata_v *cs, uint64_t s, cdata_pfmt_t pfmt) {

  if (ref_cdata_v(cs, 0)->fmt == '7') {
    fprintf(stderr, "[%s:%d] Unpack does not support format 7 chunking.\n", __func__, __LINE__);
    fflush(stderr);
    exit(1);
  }
  
  uint64_t i,m, k, kn = cs->size;
  cdata_t expanded = decompress(*ref_cdata_v(cs, 0));
  uint64_t n = expanded.n;
  cdata_t *sliced = calloc(kn, sizeof(cdata_t));
  for (m=0; m <= n/s; ++m) {
    for (k=0; k<kn; ++k) {
      expanded = decompress(*ref_cdata_v(cs, k));
      slice(&expanded, m*s, (m+1)*s-1, &sliced[k]);
    }
    for (i=0; i<sliced[0].n; ++i) {
      for (k=0; k<kn; ++k) {
        if(k) fputc('\t', stdout);
        print_cdata1(&sliced[k], i, pfmt);
      }
      fputc('\n', stdout);
    }
  }

  for (k=0; k<kn; ++k) free(sliced[k].s);
  free(expanded.s); free(sliced);
}

static void print_cdata(cdata_v *cs, cdata_pfmt_t pfmt, char *fname_row) {
  uint64_t i, k, kn = cs->size;
  cdata_t *inflated = calloc(kn, sizeof(cdata_t));
  for (k=0; k<kn; ++k) {
    cdata_t *c = ref_cdata_v(cs,k);
    if (c->fmt == '7') { // inflating format 7 is expensive
      memcpy(inflated+k, c, sizeof(cdata_t));
      inflated[k].s = malloc(c->n);
      memcpy(inflated[k].s, c->s, c->n);
    } else {
      inflated[k] = decompress(*c);
    }
  }

  uint64_t n = 0;
  if (inflated[0].fmt == '7') {
    n = fmt7_data_length(&inflated[0]);
  } else {
    n = inflated[0].n;
  }

  cdata_t cr = {0};
  if (fname_row) {
    cfile_t cf_row = open_cfile(fname_row);
    cr = read_cdata1(&cf_row);
  }
  for (i=0; i<n; ++i) {
    if (cr.s) print_cdata1(&cr, i, pfmt);
    for (k=0; k<kn; ++k) {
      if(k || cr.s) fputc('\t', stdout);
      print_cdata1(inflated+k, i, pfmt);
    }
    fputc('\n', stdout);
  }
  if (cr.s) free_cdata(&cr);
  for (k=0; k<kn; ++k) free_cdata(&inflated[k]);
  free(inflated);
}

int main_unpack(int argc, char *argv[]) {

  int c, read_all = 0, chunk = 0;
  cdata_pfmt_t pfmt = {0};
  uint64_t chunk_size = 1000000; char *fname_snames = NULL;
  int head = -1, tail = -1;
  uint8_t unit = 0; // default: auto-inferred
  int print_column_names = 0;
  char *fname_row = NULL;
  while ((c = getopt(argc, argv, "cs:l:H:T:f:u:CR:r:ah"))>=0) {
    switch (c) {
    case 'c': chunk = 1; break;
    case 's': chunk_size = atoi(optarg); break;
    case 'l': fname_snames = strdup(optarg); break;
    case 'H': head = atoi(optarg); break;
    case 'T': tail = atoi(optarg); break;
    case 'u': unit = atoi(optarg); break;
    case 'C': print_column_names = 1; break;
    case 'R': fname_row = strdup(optarg); break;
    case 'r': pfmt.ref = atoi(optarg); break;
    case 'a': read_all = 1; break;
    case 'f': pfmt.data = atoi(optarg); break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 1 > argc) { 
    usage(); 
    wzfatal("Please supply input file.\n"); 
  }

  char *fname_in = strdup(argv[optind]);
  cfile_t cf = open_cfile(fname_in);
  char *fname_index = get_fname_index(fname_in);
  index_t *idx = loadIndex(fname_index);
  if (fname_index) free(fname_index);

  snames_t snames = {0};
  if (optind + 1 < argc) {      // The requested sample names from command line
    for(int i = optind + 1; i < argc; ++i) {
      snames.s = realloc(snames.s, (snames.n+1));
      snames.s[snames.n++] = strdup(argv[i]);
    }
  } else {                      // from a file list
    snames = loadSampleNames(fname_snames, 1);
  }

  // check if we have index
  if ((tail > 0 && !idx) || (snames.n > 0 && strcmp(fname_in, "-") != 0 && !idx)) {
    fprintf(stderr, "Error, the cx file needs indexing for random sample access.\n");
    fflush(stderr);
    exit(1);
  }

  // read in the cdata
  cdata_v *cs = NULL;
  if (idx && snames.n > 0) {
    cs = read_cdata_with_snames(&cf, idx, &snames);
  } else if (read_all) {
    cs = read_cdata_all(&cf);
  } else if (head > 0) {
    cs = read_cdata_from_head(&cf, head);
  } else if (tail > 0) {
    cs = read_cdata_from_tail(&cf, idx, tail);
  } else {
    cs = read_cdata_from_head(&cf, 1);
  }

  // set unit size
  if (unit > 8 || ((unit&0x1) && unit != 1 && unit)) {
    fprintf(stderr, "[%s:%d] Unit size (%u) can only be 1,2,4,6,8.\n", __func__, __LINE__, unit);
    fflush(stderr);
    exit(1);
  }
  for (uint64_t i=0; i<cs->size; ++i) ref_cdata_v(cs, i)->unit = unit;

  int col1_is_row_index=0;
  if (ref_cdata_v(cs, 0)->fmt == '7') col1_is_row_index = 1;
  
  // output headers
  if (print_column_names) {
    if (!snames.n) {
      if (idx) {
        int n0 = 0;
        index_pair_t *idx_pairs = index_pairs(idx, &n0);
        if (read_all) {
          snames.n = n0;
          snames.s = calloc(snames.n, sizeof(char*));
          for (int i=0; i<snames.n; ++i) snames.s[i] = idx_pairs[i].key;
        } else if (head > 0) {
          snames.n = head;
          snames.s = calloc(snames.n, sizeof(char*));
          for (int i=0; i<snames.n; ++i) snames.s[i] = idx_pairs[i].key;
        } else if (tail > 0) {
          snames.n = tail;
          snames.s = calloc(snames.n, sizeof(char*));
          for (int i=0; i<tail; ++i) snames.s[i] = idx_pairs[n0-tail+i].key;
        } else {
          snames.n = 1; snames.s = calloc(1, sizeof(char*));
          snames.s[0] = idx_pairs[0].key;
        }
        free(idx_pairs);          // ownership of keys are transfered to snames.s
      } else {
        fprintf(stderr, "[%s:%d] Error, index file is missing for printing sample names.\n", __func__, __LINE__);
        fflush(stderr);
        exit(1);
      }
    }
    if (fname_row || col1_is_row_index) {
      if (pfmt.ref == 0) fputs("chrm\tbeg0\tend1", stdout);
      else if (pfmt.ref == 1) fputs("chrm\tbeg0\tend0", stdout);
      else fputs("chrm_beg1", stdout);
    }
    for (int i=0; i<snames.n; ++i) {
      if (fname_row || col1_is_row_index || i) fputc('\t', stdout);
      fputs(snames.s[i], stdout);
    }
    fputc('\n', stdout);
  }

  // output the cs
  if (chunk) print_cdata_chunk(cs, chunk_size, pfmt); // TODO: chunking is a little redundant to rowsub
  else print_cdata(cs, pfmt, fname_row);

  // clean up
  for (uint64_t i=0; i<cs->size; ++i) free_cdata(ref_cdata_v(cs,i));
  free_cdata_v(cs);
  bgzf_close(cf.fh);
  free(fname_in);
  if (idx) cleanIndex(idx);
  cleanSampleNames(&snames);
  
  return 0;
}
