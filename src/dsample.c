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
#include <time.h>
#include "cfile.h"
#include "snames.h"

static int usage(void) {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: yame dsample [options] <in.cx> [out.cx]\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Downsample methylation data for format 3 or 6.\n");
  fprintf(stderr, "  - For format 3, downsampling masks by setting M=U=0.\n");
  fprintf(stderr, "  - For format 6, downsampling masks by clearing the universe bit.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -o [PATH] output .cx file name.\n");
  fprintf(stderr, "              If missing, write to stdout (no index will be written).\n");
  fprintf(stderr, "    -s [int]  seed for random sampling (default: current time).\n");
  fprintf(stderr, "    -N [int]  number of records to sample/keep per sample (default: 100).\n");
  fprintf(stderr, "              If N >= available records, all available records are kept.\n");
  fprintf(stderr, "    -r [int]  number of downsampled replicates per input sample (default: 1).\n");
  fprintf(stderr, "              Each replicate is independently re-sampled.\n");
  fprintf(stderr, "    -h        this help.\n");
  fprintf(stderr, "\n");

  return 1;
}

/**
 * Select K unique indices (filled before the function) using a
 * partial Fisher–Yates shuffle.
 *
 * We shuffle only the first K positions in-place; after this call,
 * array[0..K-1] contains K distinct indices in [0, N).
 *
 * The caller must provide an array of length N.
 */
static void fisher_yates_shuffle_select(uint64_t *array, uint64_t N, uint64_t K) {
  if (N == 0 || K == 0) return;
  if (K > N) K = N;  // guard: cannot select more than N elements

  // Partial Fisher–Yates: randomize the first K elements
  for (uint64_t i = 0; i < K; ++i) {
    /* Scale rand() into [0, N-i) using double arithmetic to avoid
       out-of-range j and reduce bias compared to naive modulo. */
    double r = (double)rand() / ((double)RAND_MAX + 1.0);
    uint64_t range = N - i;
    uint64_t offset = (uint64_t)(r * (double)range);
    uint64_t j = i + offset; // j in [i, N-1]

    // swap array[i] and array[j]
    uint64_t tmp = array[i];
    array[i] = array[j];
    array[j] = tmp;
  }
}

/**
 * Downsample a format-3 cdata_t.
 *
 * - Input  c: compressed or uncompressed fmt3.
 * - Output cout: decompressed fmt3 with some entries masked (M=U=0).
 * - N: desired number of non-zero records to keep (per sample).
 *
 * `indices` and `to_include` are workspace buffers:
 *   - indices: length >= cout->n
 *   - to_include: bitset with at least ceil(cout->n/8) bytes
 */
cdata_t dsample_fmt3(cdata_t *c, uint64_t N,
                     uint64_t *indices, uint8_t *to_include) {
  // copy and expand into cout
  // decompress(c, cout);
  assert(c->compressed == 0);
  cdata_t cout = cdata_duplicate(*c);

  uint64_t N_indices = 0;
  for (uint64_t i = 0; i < cout.n; ++i)
    if (f3_get_mu(&cout, i))
      N_indices++;

  // If we want >= all available entries, do nothing (keep all)
  if (N >= N_indices) return cout;

  // Collect indices of eligible sites
  uint64_t j = 0;
  for (uint64_t i = 0; i < cout.n; ++i)
    if (f3_get_mu(&cout, i))
      indices[j++] = i;

  // Randomly select N of them into the first N slots of indices[]
  fisher_yates_shuffle_select(indices, N_indices, N);

  // Build bitset of which indices to keep
  uint64_t nbytes = cout.n / 8 + 1;
  memset(to_include, 0, nbytes * sizeof(uint8_t));

  for (uint64_t k = 0; k < N; ++k)
    _FMT0_SET(to_include, indices[k]);

  // Mask out everything not selected
  for (uint64_t i = 0; i < cout.n; ++i) {
    if (f3_get_mu(&cout, i) && !_FMT0_IN_SET(to_include, i)) {
      // set M=U=0
      f3_set_mu(&cout, i, 0, 0);
    }
  }
  return cout;
}

/**
 * Downsample a format-6 cdata_t.
 *
 * - Input  c: compressed or uncompressed fmt6.
 * - Output cout: decompressed fmt6 with universe bits cleared for
 *   positions outside the sampled subset.
 *
 * `indices` and `to_include` are workspace buffers:
 *   - indices: length >= cout->n
 *   - to_include: bitset with at least ceil(cout->n/8) bytes
 */
cdata_t dsample_fmt6(cdata_t *c, uint64_t N,
                  uint64_t *indices, uint8_t *to_include) {
  /* decompress(c, cout); */
  assert(c->compressed == 0);
  cdata_t cout = cdata_duplicate(*c);

  // count eligible entries in universe
  uint64_t N_indices = 0;
  for (uint64_t i = 0; i < cout.n; ++i)
    if (FMT6_IN_UNI(cout, i))
      N_indices++;

  if (N >= N_indices) return cout; // keep all

  uint64_t j = 0;
  for (uint64_t i = 0; i < cout.n; ++i)
    if (FMT6_IN_UNI(cout, i))
      indices[j++] = i;

  fisher_yates_shuffle_select(indices, N_indices, N);

  // bitset of which universe positions to keep
  uint64_t nbytes = cout.n / 8 + 1;
  memset(to_include, 0, nbytes * sizeof(uint8_t));
  for (uint64_t k = 0; k < N; ++k)
    _FMT0_SET(to_include, indices[k]);

  for (uint64_t i = 0; i < cout.n; ++i) {
    if (FMT6_IN_UNI(cout, i) && !_FMT0_IN_SET(to_include, i)) {
      // mark as NA outside the sampled subset
      FMT6_SET_NA(cout, i);
    }
  }
  return cout;
}

/**
 * Write an index for the downsampled output, taking into account replicates.
 *
 * - fname:      input .cx filename (used to load existing index, if any).
 * - fname_out:  output .cx filename; if NULL, no index is written.
 * - n_in:       number of input samples read.
 * - n_rep:      number of replicates per input sample.
 *
 * When an input index exists, its keys are used as the base names; otherwise
 * we fall back to 0-based numeric sample IDs. Replicates are suffixed with
 * `-0`, `-1`, ..., `-(n_rep-1)` when n_rep > 1.
 */
static void write_index_with_rep(char *fname,
                                 char *fname_out,
                                 int64_t n_in,
                                 int n_rep) {
  
  if (!fname_out) return; // cannot index stdout

  // load input index (if any)
  int npairs = 0;
  index_pair_t *pairs = NULL;
  char *fname_index = get_fname_index(fname);
  index_t *idx = loadIndex(fname_index);
  if (idx) pairs = index_pairs(idx, &npairs);
  free(fname_index);

  // sanity: if we have an index, its size should match n_in
  if (idx && npairs != n_in) {
    fprintf(stderr,
            "[Warning] write_index_with_rep: input index has %d entries, "
            "but n_in = %" PRId64 ". Using min of the two.\n",
            npairs, n_in);
    fflush(stderr);
    if (npairs < n_in) n_in = npairs;
  }

  // open output CX for reading to recover BGZF offsets
  cfile_t cf2 = open_cfile(fname_out);
  index_t *idx2 = kh_init(index);

  int64_t addr = bgzf_tell(cf2.fh);
  cdata_t c_tmp = (cdata_t){0};

  for (int64_t i = 0; i < n_in; ++i) {
    for (int j = 0; j < n_rep; ++j) {

      if (!read_cdata2(&cf2, &c_tmp)) {
        fprintf(stderr,
                "[Error] write_index_with_rep: data is shorter than "
                "expected (%" PRId64 " records, n_rep=%d).\n",
                n_in, n_rep);
        fflush(stderr);
        exit(1);
      }

      // decide the base name
      const char *base = NULL;
      char buf[64]; // 64 characters is enough for n_in

      if (pairs) {
        base = pairs[i].key;           // input index key
      } else {
        // fallback: use i if input index is missing
        snprintf(buf, sizeof(buf), "%" PRId64, i);
        base = buf;
      }

      // construct the final key string
      char *s = NULL;
      if (n_rep == 1) {
        int n = snprintf(NULL, 0, "%s", base);
        s = (char *)malloc(n + 1);
        snprintf(s, n + 1, "%s", base);
      } else {
        int n = snprintf(NULL, 0, "%s-%d", base, j);
        s = (char *)malloc(n + 1);
        snprintf(s, n + 1, "%s-%d", base, j);
      }

      insert_index(idx2, s, addr);

      // addr for the *next* record
      addr = bgzf_tell(cf2.fh);
    }
  }

  free_cdata(&c_tmp);
  bgzf_close(cf2.fh);

  // write output index file
  char *fname_index2 = get_fname_index(fname_out);
  FILE *out = fopen(fname_index2, "w");
  if (!out) {
    fprintf(stderr,
            "[Error] write_index_with_rep: cannot open index file %s\n",
            fname_index2);
    fflush(stderr);
    exit(1);
  }
  writeIndex(out, idx2);
  fclose(out);
  free(fname_index2);

  cleanIndex(idx2); // free the keys too
  if (pairs) free(pairs);
  if (idx)   cleanIndex(idx);   // or freeIndex(idx) if needed
}

int main_dsample(int argc, char *argv[]) {

  int c;
  unsigned seed = (unsigned) time(NULL);
  uint64_t N = 100;
  int n_rep = 1;
  char *fname_out = NULL;
  while ((c = getopt(argc, argv, "o:r:s:N:h"))>=0) {
    switch (c) {
    case 'o': fname_out = strdup(optarg); break;
    case 'r': n_rep = atoi(optarg); break;
    case 's': seed = (unsigned) strtoul(optarg, NULL, 10); break;
    case 'N': N = strtoul(optarg, NULL, 10); break;
    case 'h': return usage(); break;
    default:
      usage();
      wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 1 > argc) {
    usage(); 
    wzfatal("Please supply input file.\n");
  }

  char *fname = argv[optind];
  BGZF* fp_out;
  if (fname_out) fp_out = bgzf_open2(fname_out, "wb");
  else fp_out = bgzf_dopen(fileno(stdout), "wb");

  if (fp_out == NULL) {
    fprintf(stderr, "[%s:%d] Error opening file for writing: %s\n",
            __func__, __LINE__, fname_out ? fname_out : "<stdout>");
    fflush(stderr);
    exit(1);
  }
  
  cfile_t cf = open_cfile(fname);
  srand(seed);
  
  uint64_t *indices = NULL;
  uint8_t *to_include = NULL;

  int64_t n_samples = 0;
  
  for (;;++n_samples) {
    cdata_t c_in = read_cdata1(&cf);
    if (c_in.n == 0) {
      free_cdata(&c_in);
      break;  // end-of-file
    }
    
    decompress2(&c_in);
    uint64_t n_positions = c_in.n;

    if (!indices) {
      indices = (uint64_t *) calloc(n_positions, sizeof(uint64_t));
      if (!indices) {
        wzfatal("Failed to allocate indices buffer.");
      }
      uint64_t nbytes = (n_positions + 7) / 8;
      to_include = (uint8_t *) calloc(nbytes, 1);
      if (!to_include) {
        wzfatal("Failed to allocate bitset buffer.");
      }
    }

    /* For each requested replicate, independently downsample and write. */
    for (int rep = 0; rep < n_rep; ++rep) {
      cdata_t c_out = {0};

      switch (c_in.fmt) {
      case '3':
        c_out = dsample_fmt3(&c_in, N, indices, to_include);
        break;
      case '6':
        c_out = dsample_fmt6(&c_in, N, indices, to_include);
        break;
      default:
        wzfatal("Format %d not recognized (only 3 and 6 are supported).\n", c_in.fmt);
      }

      fflush(stderr);
      if (!c_out.compressed) {
        cdata_compress(&c_out);
      }
      fflush(stderr);
      cdata_write1(fp_out, &c_out);
      free_cdata(&c_out);
    }
    free_cdata(&c_in);
  }

  free(indices);
  free(to_include);
  bgzf_close(cf.fh);
  bgzf_close(fp_out);

  /* Only attempt to write an index when we actually wrote to a regular file. */
  write_index_with_rep(fname, fname_out, n_samples, n_rep);
  if (fname_out) free(fname_out);
   
  return 0;
}


