#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include "cfile.h"
#include "snames.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: yame dsample [options] <in.cx> <out.cx>\n");
  fprintf(stderr, "Take format 6 and format 3 as input.\n");
  fprintf(stderr, "If format 3, this function mask by setting M=U=0.\n");
  fprintf(stderr, "If format 6, this function mask by setting universe bit.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -o        output cg file name. if missing, output to stdout.\n");
  fprintf(stderr, "    -s [int]  seed for random sampling (default: use time to seed) .\n");
  fprintf(stderr, "    -N [int]  number of records to sample/keep. Rest will be masked (default: 100).\n");
  fprintf(stderr, "    -r [int]  number of replicates per input (default: 1).\n");
  fprintf(stderr, "              When higher than available, capped to available.\n");
  fprintf(stderr, "    -h        This help\n");
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

void dsample_fmt3(cdata_t *c, cdata_t *cout, uint64_t N,
                  uint64_t *indices, uint8_t *to_include) {
  // copy and expand into cout
  decompress(c, cout);

  uint64_t N_indices = 0;
  for (uint64_t i = 0; i < cout->n; ++i)
    if (f3_get_mu(cout, i))
      N_indices++;

  // If we want >= all available entries, do nothing (keep all)
  if (N >= N_indices) return;

  // Collect indices of eligible sites
  uint64_t j = 0;
  for (uint64_t i = 0; i < cout->n; ++i)
    if (f3_get_mu(cout, i))
      indices[j++] = i;

  // Randomly select N of them into the first N slots of indices[]
  fisher_yates_shuffle_select(indices, N_indices, N);

  // Build bitset of which indices to keep
  uint64_t nbytes = cout->n / 8 + 1;
  memset(to_include, 0, nbytes * sizeof(uint8_t));

  for (uint64_t k = 0; k < N; ++k)
    _FMT0_SET(to_include, indices[k]);

  // Mask out everything not selected
  for (uint64_t i = 0; i < cout->n; ++i) {
    if (f3_get_mu(cout, i) && !_FMT0_IN_SET(to_include, i)) {
      // set M=U=0
      f3_set_mu(cout, i, 0, 0);
    }
  }
}

void dsample_fmt6(cdata_t *c, cdata_t *cout, uint64_t N,
                  uint64_t *indices, uint8_t *to_include) {
  decompress(c, cout);

  uint64_t N_indices = 0;
  for (uint64_t i = 0; i < cout->n; ++i)
    if (FMT6_IN_UNI(*cout, i))
      N_indices++;

  if (N >= N_indices) return; // keep all

  uint64_t j = 0;
  for (uint64_t i = 0; i < cout->n; ++i)
    if (FMT6_IN_UNI(*cout, i))
      indices[j++] = i;

  fisher_yates_shuffle_select(indices, N_indices, N);

  // bitset of which universe positions to keep
  uint64_t nbytes = cout->n / 8 + 1;
  memset(to_include, 0, nbytes * sizeof(uint8_t));
  for (uint64_t k = 0; k < N; ++k)
    _FMT0_SET(to_include, indices[k]);

  for (uint64_t i = 0; i < cout->n; ++i) {
    if (FMT6_IN_UNI(*cout, i) && !_FMT0_IN_SET(to_include, i)) {
      // mark as NA outside the sampled subset
      FMT6_SET_NA(*cout, i);
    }
  }
}

void write_index_with_rep(char *fname, char *fname_out, int64_t n_in, int n_rep) {
    if (!fname_out) return;

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

    freeIndex(idx2);
    if (pairs) free(pairs);
    if (idx)   cleanIndex(idx);   // or freeIndex(idx) if needed
}

int main_dsample(int argc, char *argv[]) {

  int c; char *fname_out = NULL; 
  unsigned seed = (unsigned int)time(NULL);
  uint64_t N = 100; int n_rep = 1;
  while ((c = getopt(argc, argv, "o:r:s:N:h"))>=0) {
    switch (c) {
    case 'o': fname_out = strdup(optarg); break;
    case 'r': n_rep = atoi(optarg); break;
    case 's': seed = strtoul(optarg, NULL, 10); break;
    case 'N': N = strtoul(optarg, NULL, 10); break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 1 > argc) {
    usage(); 
    wzfatal("Please supply input file.\n");
  }

  char *fname = argv[optind];
  if (argc >= optind + 2)
    fname_out = strdup(argv[optind+1]);
  
  BGZF* fp;
  if (fname_out) fp = bgzf_open2(fname_out, "wb");
  else fp = bgzf_dopen(fileno(stdout), "wb");
  if (fp == NULL) {
    fprintf(stderr, "[%s:%d] Error opening file for writing: %s\n", __func__, __LINE__, fname_out);
    fflush(stderr);
    exit(1);
  }
  
  cfile_t cf = open_cfile(fname);
  srand(seed);
  uint64_t *indices = NULL;
  uint8_t *to_include = NULL;

  uint64_t kq = 0;
  for (kq=0;;++kq) {
    cdata_t c = read_cdata1(&cf);
    if (c.n == 0) break;
    cdata_t c2 = {0};

    if (!indices) {
      indices = calloc(c.n, sizeof(uint64_t));
      to_include = calloc(c.n/8+1, 1);
    }

    for (int k=0; k<n_rep; ++k) {
      switch (c.fmt) {
      case '3': dsample_fmt3(&c, &c2, N, indices, to_include); break;
      case '6': dsample_fmt6(&c, &c2, N, indices, to_include); break;
      default: wzfatal("Format %d not recognized.\n", c.fmt);
      }
      if (!c2.compressed) cdata_compress(&c2);
      cdata_write1(fp, &c2);
    }
    
    free_cdata(&c); free_cdata(&c2);
  }
  free(indices); free(to_include);
  write_index_with_rep(fname, fname_out, kq, n_rep);
  
  bgzf_close(cf.fh);
  bgzf_close(fp);
  if (fname_out) free(fname_out);
  
  return 0;
}


