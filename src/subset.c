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

void subset_fmt2_states(cfile_t cf, snames_t snames, char *fname_out) {
  cdata_t c = read_cdata1(&cf);
  decompress2(&c);
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
  cdata_t c0 = {.n = c.n, .fmt = '0'}; // output data
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

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: yame subset [options] <in.cx> [<sample1> <sample2> ...]\n");
  fprintf(stderr, "If -o <out.cx>, an index will also be generated. Otherwise, output .cx to stdout without index.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -v        verbose\n");
  fprintf(stderr, "    -o        output cx file name. if missing, output to stdout without index.\n");
  fprintf(stderr, "    -l        Path to the sample list. Ignored if sample names are provided on the command line.\n");
  fprintf(stderr, "    -s        Filter format 2 <in.cx> instead of samples in files. Output format 0.\n");
  fprintf(stderr, "    -H [N]    Process N samples from the start of the list, where N is less than or equal to the total number of samples.\n");
  fprintf(stderr, "    -T [N]    Process N samples from the end of the list, where N is less than or equal to the total number of samples. Requires index.\n");
  fprintf(stderr, "    -h        This help\n");
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

