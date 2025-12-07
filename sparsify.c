#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include "cfile.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: yame sparsify [options] <binary.cg>\n");
  fprintf(stderr, "Take format 6 data as input, masking data with NA.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -p [float] masking probability (default: 0.8, masking 80%% data).\n");
  fprintf(stderr, "    -r [int]   number of replicates per input (default: 1).\n");
  fprintf(stderr, "    -s [int]   seed for random sampling. (default: use time to seed).\n");
  fprintf(stderr, "    -o         output cx file name (format 6). if missing, output to stdout.\n");
  fprintf(stderr, "    -h         This help\n");
  fprintf(stderr, "\n");

  return 1;
}

/**
 * Fisher-Yates shuffle implementation.
 * It selects K unique indices from the N available indices.
 *
 * @param array The array of indices (0 to N-1).
 * @param N The total number of elements/indices (e.g., 1,000,000).
 * @param K The number of elements to select (e.g., 100,000).
 */
void fisher_yates_shuffle_select(uint64_t *array, uint64_t N, uint64_t K) {
    if (K > N) K = N; // Cannot select more than N elements

    // initialize
    uint64_t i;
    for (i = 0; i < N; ++i) array[i] = i;

    for (i = 0; i < K; i++) {
        // Generate a random index j between i and N-1 (inclusive)
        // Note: Using rand() here is fast because we only do K calls, not N.
        // It's also safer than (rand() % (N - i)) which can have bias if (N-i) is large.
        uint64_t max_rand = RAND_MAX;
        uint64_t r = rand();
        uint64_t j = i + r / (max_rand / (N - i) + 1);

        // Ensure we handle edge case where rand() returns RAND_MAX
        if (j >= N) j = N - 1; 

        // fprintf(stderr, "%ld\t%ld\n", i, j);
        // Swap the elements at array[i] and array[j]
        uint64_t temp = array[i];
        array[i] = array[j];
        array[j] = temp;
    }
}

int main_sparsify(int argc, char *argv[]) {

  int c; int n_rep = 1;
  char *fname_out = NULL;
  double p_mask = 0.8;
  int seed = 0;
  while ((c = getopt(argc, argv, "o:r:p:s:h"))>=0) {
    switch (c) {
    case 'o': fname_out = strdup(optarg); break;
    case 'r': n_rep = atoi(optarg); break;
    case 'p': p_mask = atof(optarg); break;
    case 's': seed = atoi(optarg); break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 1 > argc) {
    usage(); 
    wzfatal("Please supply input file.\n");
  }

  char *fname = argv[optind];

  BGZF *fp_out;
  if (fname_out) fp_out = bgzf_open2(fname_out, "w");
  else fp_out = bgzf_dopen(fileno(stdout), "w");
  if (fp_out == NULL) {
    fprintf(stderr, "Error opening file for writing: %s\n", fname_out);
    exit(1);
  }

  cfile_t cf = open_cfile(fname);
  char *fname_index = get_fname_index(fname);
  index_t *idx = loadIndex(fname_index);
  free(fname_index);
  if (seed > 0) srand(seed);
  else srand(time(NULL));

  int n_in = 0;
  uint8_t *mask = NULL;
  uint64_t k_selected = 0;
  uint64_t *indices_selected = NULL;
  
  while (1) {
    cdata_t c = read_cdata1(&cf);
    if (c.n == 0) break;
    n_in++;
    decompress2(&c);
    if (c.fmt != '6') {
      fprintf(stderr, "[%s:%d] Only format 6 (given %d) files are supported.\n", __func__, __LINE__, c.fmt);
      fflush(stderr);
      exit(1);
    }

    for (int j=0; j<n_rep; ++j) {
      cdata_t cout = {.fmt = '6', .n = c.n};
      cout.s = calloc((cout.n+3)/4, sizeof(uint8_t));
      memcpy(cout.s, c.s, (cout.n+3)/4);

      if (!mask) {
        mask = calloc(c.n, sizeof(uint8_t));
        indices_selected = calloc(c.n, sizeof(uint64_t));
        k_selected = floor((1-p_mask)*c.n);
      }

      fisher_yates_shuffle_select(indices_selected, c.n, k_selected);
      for (uint64_t ii=0; ii<cout.n; ++ii) mask[ii] = 1;
      for (uint64_t ii=0; ii<k_selected; ++ii) mask[indices_selected[ii]] = 0;
      for (uint64_t i=0; i<cout.n; ++i) {
        if (mask[i]) {
          FMT6_SET_NA(cout, i);
        }
      }
      cdata_compress(&cout);
      cdata_write1(fp_out, &cout);
      free_cdata(&cout);
    }
    free_cdata(&c);
  }

  if (mask) free(mask);
  if (indices_selected) free(indices_selected);

  if (fname_out) {              // output index
    int npairs = 0;
    index_pair_t *pairs = NULL;
    if (idx) pairs = index_pairs(idx, &npairs);
    
    cfile_t cf2 = open_cfile(fname_out);
    index_t *idx2 = kh_init(index);
    int64_t addr = bgzf_tell(cf2.fh);
    cdata_t c_tmp = {0};
    for (int i=0; i< n_in; ++i) {
      for (int j=0; j < n_rep; ++j) {
        if (!read_cdata2(&cf2, &c_tmp)) {
          fprintf(stderr, "[Error] Data is shorter than the sample name list.\n");
          fflush(stderr);
          exit(1);
        }
        if (n_rep == 1) {
          if (pairs) insert_index(idx2, pairs[i].key, addr);
          else {
            int n = snprintf(NULL, 0, "%d", i);
            char *s = (char*) malloc(n+1);
            snprintf(s, n+1, "%d", i);
            insert_index(idx2, s, addr);
          }
        } else {
          if (pairs) {          // use names in the input index
            int n = snprintf(NULL, 0, "%s-%d", pairs[i].key, j);
            char *s = (char*) malloc(n+1);
            snprintf(s, n+1, "%s-%d", pairs[i].key, j);
            insert_index(idx2, s, addr);
          } else {        // use i if input index is missing, e.g., from a pipe
            int n = snprintf(NULL, 0, "%d-%d", i, j);
            char *s = (char*) malloc(n+1);
            snprintf(s, n+1, "%d-%d", i, j);
            insert_index(idx2, s, addr);
          }
        }
        addr = bgzf_tell(cf2.fh);
      }
    }
    free_cdata(&c_tmp);

    char *fname_index2 = get_fname_index(fname_out);
    FILE *out = fopen(fname_index2, "w");
    writeIndex(out, idx2);
    fclose(out);
    free(fname_index2);
    bgzf_close(cf2.fh);
    freeIndex(idx2);
    free(pairs);
    cleanIndex(idx);
  }

  if (fname_out) free(fname_out);
  bgzf_close(fp_out);
  bgzf_close(cf.fh);

  return 0;
}
