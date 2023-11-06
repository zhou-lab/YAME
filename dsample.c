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
  fprintf(stderr, "    -N [N]    number of records to sample. (default: 100).\n");
  fprintf(stderr, "              When higher than available, capped to available.\n");
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

cdata_t* downsample(cdata_t *c0, cdata_t *c, uint64_t N) {

  if (c0->fmt != '3') {
    fprintf(stderr, "[%s:%d] Input must be of format 3.\n", __func__, __LINE__);
    fflush(stderr);
    exit(1);
  }
  
  decompress(c0, c);
  uint64_t N_indices = 0;
  for (uint64_t i=0; i<c->n; ++i)
    if (f3_get_mu(c, i)) N_indices++;

  if (N >= N_indices) return c; // more elements to sample than exist
  
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
  free(indices); free(to_include);
  return c;
}

int main_dsample(int argc, char *argv[]) {

  int c, verbose = 0;
  unsigned seed = (unsigned int)time(NULL);
  uint64_t N = 100;
  while ((c = getopt(argc, argv, "vs:N:h"))>=0) {
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
    cdata_t c2 = {0};
    downsample(&c, &c2, N);
    cdata_write(fname_out, &c2, "wb", verbose);
    free_cdata(&c); free_cdata(&c2);
  }
  
  bgzf_close(cf.fh);
  if (fname_out) free(fname_out);
  
  return 0;
}


