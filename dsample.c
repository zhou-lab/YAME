#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include "cfile.h"
#include "snames.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: yame dsample [options] <in.cx> <out.cx>\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -v        verbose\n");
  fprintf(stderr, "    -s [N]    seed for sampling.\n");
  fprintf(stderr, "    -N [N]    number of cg to sample.\n");
  fprintf(stderr, "    -h        This help\n");
  fprintf(stderr, "\n");

  return 1;
}

void Fisher_Yates_shuffle(uint64_t *array, uint64_t n) {
  for (uint64_t i = 0; i<n-1; ++i) {
    uint64_t j = i + rand() / (RAND_MAX / (n-i) + 1);
    uint64_t t = array[j];
    array[j] = array[i];
    array[i] = t;
  }
}

void downsample(cdata_t *c, uint64_t N) {
  cdata_t c2 = {0};
  decompress(c, &c2);
  uint64_t N_indices = 0;
  for (uint64_t i=0; i<c->n; ++i)
    if (f3_get_mu(c, i)) N_indices++;

  if (N >= N_indices) return; // more elements to sample than exist
  
  uint64_t *indices = calloc(N_indices, sizeof(uint64_t));
  uint64_t j = 0;
  for (uint64_t i = 0; i<c->n; ++i)
    if (f3_get_mu(c, i)) indices[j++] = i;

  Fisher_Yates_shuffle(indices, N_indices);

  uint8_t *to_include = calloc(c->n/8+1, 1);
  for (uint64_t j = 0; j<N; ++j)
    to_include[indices[j]>>3] |= (1<<(indices[j]&0x7));

  for (uint64_t i = 0; i<c->n; ++i) {
    if (f3_get_mu(c, i) && !(to_include[i>>3]&(1<<(i&0x7)))) {
      f3_set_mu(c, i, 0, 0);
    }
  }
}

int main_dsample(int argc, char *argv[]) {

  int c, verbose = 0;
  unsigned seed = (unsigned int)time(NULL);
  uint64_t N = 0;
  while ((c = getopt(argc, argv, "vo:c:h"))>=0) {
    switch (c) {
    case 'v': verbose = 1; break;
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
  char *fname_out = NULL;
  if (argc >= optind + 2)
    fname_out = strdup(argv[optind+1]);

  cfile_t cf = open_cfile(fname);
  srand(seed);
  for (uint64_t kq=0;;++kq) {
    cdata_t c = read_cdata1(&cf);
    if (c.n == 0) break;
    downsample(&c, N);
    cdata_write(fname_out, &c, "wb", verbose);
  }
  
  bgzf_close(cf.fh);
  if (fname_out) free(fname_out);
  
  return 0;
}


