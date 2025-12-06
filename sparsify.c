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
  fprintf(stderr, "    -o         output cx file name (format 6). if missing, output to stdout.\n");
  fprintf(stderr, "    -h         This help\n");
  fprintf(stderr, "\n");

  return 1;
}

int main_sparsify(int argc, char *argv[]) {

  int c;
  char *fname_out = NULL;
  double p_mask = 0.8;
  while ((c = getopt(argc, argv, "o:p:h"))>=0) {
    switch (c) {
    case 'o': fname_out = strdup(optarg); break;
    case 'p': p_mask = atof(optarg); break;
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
  srand(time(NULL));
  
  while (1) {
    cdata_t c = read_cdata1(&cf);
    if (c.n == 0) break;
    decompress2(&c);
    if (c.fmt != '6') {
      fprintf(stderr, "[%s:%d] Only format 6 (given %d) files are supported.\n", __func__, __LINE__, c.fmt);
      fflush(stderr);
      exit(1);
    }

    cdata_t cout = {.fmt = '6', .n = c.n};
    cout.s = calloc((cout.n+3)/4, sizeof(uint8_t));
    memcpy(cout.s, c.s, (cout.n+3)/4);
    for (uint64_t i=0; i<cout.n; ++i) {
      // Calculate a random float R in the range [0.0, 1.0)
      double R = (double)rand() / (RAND_MAX + 1.0);
      if (R < p_mask) { // THIS IS PROBABLY WRONG
        FMT6_SET_NA(cout, i);
      }
    }
    cdata_compress(&cout);
    cdata_write1(fp_out, &cout);
    free_cdata(&cout);
    free_cdata(&c);
  }

  if (idx && fname_out) {              // output index
    int npairs = 0;
    index_pair_t *pairs = index_pairs(idx, &npairs);
    
    cfile_t cf2 = open_cfile(fname_out);
    index_t *idx2 = kh_init(index);
    int64_t addr = bgzf_tell(cf2.fh);
    cdata_t c_tmp = {0};
    for (int i=0; i< npairs; ++i) {
      if (!read_cdata2(&cf2, &c_tmp)) {
        fprintf(stderr, "[Error] Data is shorter than the sample name list.\n");
        fflush(stderr);
        exit(1);
      }
      insert_index(idx2, pairs[i].key, addr);
      addr = bgzf_tell(cf2.fh);
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
