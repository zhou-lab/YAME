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

#include "cfile.h"

/*
 * yame subset
 * ===========
 *
 * Overview
 * --------
 * This subcommand has two distinct modes:
 *
 *   (A) Sample subsetting (default): subset a multi-sample .cx by sample name.
 *   (B) Format-2 state subsetting (-s): convert a single fmt2 state vector into
 *       one binary (fmt0) track per requested state term.
 *
 * The requested names can be provided as trailing command-line arguments or
 * via -l <list.txt> (one name per line). If no names are provided, the default
 * is to take head/tail samples from the input index (see below).
 *
 *
 * Mode A: subset_samples()
 * -----------------------
 * Purpose:
 *   Extract specific samples (cdata records) from a .cx stream by name.
 *
 * Requirements:
 *   The input must be indexed (.cxi). The index maps sample name -> BGZF byte offset.
 *   Without an index, random access via bgzf_seek() is not possible, so the code
 *   errors out.
 *
 * Name selection:
 *   - If snames.n > 0 (explicit list), use that list in the given order.
 *   - Else, build a list from the input index:
 *       * if tail > 0: take last  tail samples
 *       * else:        take first head samples (head defaults to 1 if < 1)
 *
 * Extraction:
 *   For each requested sample name:
 *     - lookup offset = getIndex(idx, name)
 *     - bgzf_seek(cf.fh, offset, SEEK_SET)
 *     - read_cdata2(&cf, &c)
 *     - cdata_write1(fp_out, &c)
 *
 * Output indexing:
 *   If -o is provided (writing to a file), we generate an output index by
 *   re-reading the output file sequentially and recording bgzf_tell() offsets
 *   for each emitted record, keyed by the corresponding sample name.
 *   If output is stdout, no index is written.
 *
 *
 * Mode B: subset_fmt2_states()  [enabled by -s]
 * ------------------------------------------------
 * Purpose:
 *   Given a single format-2 (state) dataset, produce a separate binary (fmt0)
 *   mask for each requested term/state name.
 *
 * Input expectations:
 *   - Reads ONE cdata record from the input stream and requires c.fmt == '2'.
 *   - The fmt2 record is decompressed (decompress_in_situ).
 *   - fmt2_set_aux() is called if needed to populate aux->keys[] (the list of
 *     state labels/terms).
 *
 * Term resolution:
 *   For each requested term name:
 *     - scan aux->keys[] to find a unique matching key index (i_term)
 *     - error if missing or if multiple matches occur
 *
 * Mask creation:
 *   Creates a fmt0 output vector c0 with length c.n (same number of rows as input).
 *   For each row ii:
 *     - if f2_get_uint64(&c, ii) == i_term, set bit ii in c0
 *   Writes c0 as one output record per requested term.
 *
 * Output indexing:
 *   If -o is used, an output index is written mapping term name -> record offset,
 *   using the same “re-read output file and bgzf_tell() offsets” approach.
 *
 *
 * Practical notes / gotchas
 * -------------------------
 * - In default sample-subset mode, the tool is fundamentally index-driven.
 * - In -s mode, the term list is taken from the fmt2 key dictionary (aux->keys),
 *   and the output contains *multiple* fmt0 records (one per requested term).
 * - If you pass both explicit names and -H/-T, the explicit list wins because
 *   snames.n is non-zero and head/tail logic is skipped.
 */

void subset_fmt2_states(cfile_t cf, snames_t snames, char *fname_out) {
  cdata_t c = read_cdata1(&cf);
  decompress_in_situ(&c);
  if (c.fmt != '2') {
    wzfatal("To subset states, please provide a format 2 input. Give %c", c.fmt);
  }

  // output
  BGZF *fp;
  if (fname_out) fp = bgzf_open2(fname_out, "w");
  else fp = bgzf_dopen(fileno(stdout), "w");
  if (fp == NULL) {
    fprintf(stderr, "Error opening file for writing: %s\n", fname_out);
    exit(1);
  }

  if (!c.aux) fmt2_set_aux(&c);
  f2_aux_t *aux = (f2_aux_t*) c.aux;
  cdata_t c0 = {.n = c.n, .fmt = '0', .compressed=1}; // output data
  c0.s = calloc(cdata_nbytes(&c0), 1);
  for (int64_t i = 0; i<snames.n; ++i) {
    uint64_t i_term = 0; int found = 0;
    for (uint64_t j = 0; j<aux->nk; ++j) {
      if (strcmp(snames.s[i], aux->keys[j]) == 0) {
        if (found) wzfatal("Multiple match for %s.", snames.s[i]);
        else i_term = j;
        found = 1;
      }
    }
    if (!found) wzfatal("Cannot find term %s.", snames.s[i]);

    memset(c0.s, 0, cdata_nbytes(&c0));
    for (uint64_t ii = 0; ii < c.n; ++ii) {
      if (f2_get_uint64(&c, ii) == i_term) FMT0_SET(c0, ii);
    }
    cdata_write1(fp, &c0);
  }
  free_cdata(&c);
  free_cdata(&c0);
  bgzf_close(fp);
  
  if (fname_out) {              // output index
    cfile_t cf2 = open_cfile(fname_out);
    index_t *idx2 = kh_init(index);
    int64_t addr = bgzf_tell(cf2.fh);
    cdata_t c_tmp = {0};
    for (int i=0; i< snames.n; ++i) {
      if (!read_cdata2(&cf2, &c_tmp)) {
        fprintf(stderr, "[Error] Data is shorter than the sample name list.\n");
        fflush(stderr);
        exit(1);
      }
      insert_index(idx2, snames.s[i], addr);
      addr = bgzf_tell(cf2.fh);
    }
    free_cdata(&c_tmp);

    char *fname_index2 = get_fname_index(fname_out);
    FILE *out = fopen(fname_index2, "w");
    writeIndex(out, idx2);
    fclose(out);
    free(fname_index2);
    free(fname_out);
    bgzf_close(cf2.fh);
    freeIndex(idx2);
  }
}

void subset_samples(cfile_t cf, index_t *idx, snames_t snames, char *fname_out, int head, int tail) {

  // check if we have index
  if (!idx) {
    fprintf(stderr, "Error, the cx file needs indexing for subsetting.\n");
    fflush(stderr);
    exit(1);
  }

  // if sample names are not explicitly given, use head and tails
  if (snames.n == 0) {
    int npairs = 0;
    index_pair_t *pairs = index_pairs(idx, &npairs);
    if (tail > 0) {
      if (tail > npairs) tail = npairs;
      for (int i=0; i<tail; ++i) {
        snames.s = realloc(snames.s, (snames.n+1)*sizeof(const char*));
        snames.s[snames.n++] = strdup(pairs[npairs-tail+i].key);
      }
    } else {
      if (head < 1) head = 1;   // default to head 1
      if (head > npairs) head = npairs;
      for (int i=0; i<head; ++i) {
        snames.s = realloc(snames.s, (snames.n+1)*sizeof(const char*));
        snames.s[snames.n++] = strdup(pairs[i].key);
      }
    }
    clean_index_pairs(pairs, npairs);
  }

  // output
  BGZF *fp;
  if (fname_out) fp = bgzf_open2(fname_out, "w");
  else fp = bgzf_dopen(fileno(stdout), "w");
  if (fp == NULL) {
    fprintf(stderr, "Error opening file for writing: %s\n", fname_out);
    exit(1);
  }
  cdata_t c = {0};              // output data
  for (int i=0; i<snames.n; ++i) {
    int64_t index = getIndex(idx, snames.s[i]);
    assert(index >= 0);
    assert(bgzf_seek(cf.fh, index, SEEK_SET) == 0);
    read_cdata2(&cf, &c);
    if (c.n <= 0) {
      fprintf(stderr, "[%s:%d] Error, cannot find %s.\n", __func__, __LINE__, snames.s[i]);
      fflush(stderr);
      exit(1);
    }
    cdata_write1(fp, &c);
  }
  bgzf_close(fp);

  if (fname_out) {              // output index
    cfile_t cf2 = open_cfile(fname_out);
    index_t *idx2 = kh_init(index);
    int64_t addr = bgzf_tell(cf2.fh);
    for (int i=0; i< snames.n; ++i) {
      if (!read_cdata2(&cf2, &c)) {
        fprintf(stderr, "[Error] Data is shorter than the sample name list.\n");
        fflush(stderr);
        exit(1);
      }
      insert_index(idx2, snames.s[i], addr);
      addr = bgzf_tell(cf2.fh);
    }
    free(c.s);

    char *fname_index2 = get_fname_index(fname_out);
    FILE *out = fopen(fname_index2, "w");
    writeIndex(out, idx2);
    fclose(out);
    free(fname_index2);
    free(fname_out);
    bgzf_close(cf2.fh);
    freeIndex(idx2);
  }
}

static int usage(void) {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "  yame subset [options] <in.cx> [sample1 sample2 ...] > out.cx\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Purpose:\n");
  fprintf(stderr, "  Subset a multi-sample .cx by sample names (requires an index), or\n");
  fprintf(stderr, "  (with -s) convert a format-2 state track into one binary track per state.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Modes:\n");
  fprintf(stderr, "  (A) Sample subsetting (default):\n");
  fprintf(stderr, "      Select named samples from <in.cx> and emit them in the given order.\n");
  fprintf(stderr, "      Requires <in.cx>.cxi index.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  (B) Subset format-2 states (-s):\n");
  fprintf(stderr, "      Interpret <in.cx> as a single format-2 dataset (must be fmt2).\n");
  fprintf(stderr, "      For each requested state name, emit one format-0 bitset where\n");
  fprintf(stderr, "      bit=1 iff row state == that term.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Input sample list:\n");
  fprintf(stderr, "  Provide sample names either:\n");
  fprintf(stderr, "    * as trailing arguments on the command line, OR\n");
  fprintf(stderr, "    * via -l <list.txt> (one name per line).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "  -o <out.cx>  Write output to a file. If provided, an output index (.cxi)\n");
  fprintf(stderr, "              is also generated. If omitted, writes to stdout (no index).\n");
  fprintf(stderr, "  -l <list>    Path to sample/state list. Ignored if names are provided as\n");
  fprintf(stderr, "              trailing command-line arguments.\n");
  fprintf(stderr, "  -s           Format-2 state filtering mode (output format 0; one record per term).\n");
  fprintf(stderr, "  -H <N>       If no names are provided, take the first N samples from the input index.\n");
  fprintf(stderr, "  -T <N>       If no names are provided, take the last  N samples from the input index.\n");
  fprintf(stderr, "  -h           Show this help message.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Notes:\n");
  fprintf(stderr, "  * -H/-T only apply when you did NOT provide an explicit name list.\n");
  fprintf(stderr, "  * -T requires an index (same as default sample subsetting).\n");
  fprintf(stderr, "  * In -s mode, the input is expected to be a single fmt2 record; the output\n");
  fprintf(stderr, "    contains one fmt0 record per requested term/state.\n");
  fprintf(stderr, "\n");
  return 1;
}

int main_subset(int argc, char *argv[]) {

  int c0; char *fname_snames = NULL;
  int head = -1, tail = -1;
  char *fname_out = NULL;
  int filter_fmt2_states = 0;
  while ((c0 = getopt(argc, argv, "o:l:sH:T:h"))>=0) {
    switch (c0) {
    case 'o': fname_out = strdup(optarg); break;
    case 'l': fname_snames = strdup(optarg); break;
    case 's': filter_fmt2_states = 1; break;
    case 'H': head = atoi(optarg); break;
    case 'T': tail = atoi(optarg); break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c0);
    }
  }

  if (argc < optind + 1) { 
    usage(); 
    wzfatal("Please supply input file and output file.\n"); 
  }

  // input
  cfile_t cf = open_cfile(argv[optind]);
  char *fname_index = get_fname_index(argv[optind]);
  index_t *idx = loadIndex(fname_index);
  free(fname_index);
  optind++;

  // get sample names
  snames_t snames = {0};
  if (optind < argc) {          // sample names from command line
    for(int i = optind; i < argc; ++i) {
      snames.s = realloc(snames.s, (snames.n+1)*sizeof(const char*));
      snames.s[snames.n++] = strdup(argv[i]);
    }
  } else {                      // from a file list
    snames = loadSampleNames(fname_snames, 1);
  }

  if (filter_fmt2_states) {
    subset_fmt2_states(cf, snames, fname_out);
  } else {
    subset_samples(cf, idx, snames, fname_out, head, tail);
  }
  
  // clean up
  bgzf_close(cf.fh);
  cleanIndex(idx);
  cleanSampleNames(&snames);
  
  return 0;
}

