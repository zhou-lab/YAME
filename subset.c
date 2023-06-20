#include "cgfile.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: kycg subset [options] <in.cg> <out.cg>\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -v        verbose\n");
  fprintf(stderr, "    -l        Path to the sample list. Ignored if sample names are provided on the command line.\n");
  fprintf(stderr, "    -H [N]    Process N samples from the start of the list, where N is less than or equal to the total number of samples.\n");
  fprintf(stderr, "    -T [N]    Process N samples from the end of the list, where N is less than or equal to the total number of samples. Requires index.\n");
  fprintf(stderr, "    -n        sample name list\n");
  fprintf(stderr, "    -h        This help\n");
  fprintf(stderr, "\n");

  return 1;
}

int main_subset(int argc, char *argv[]) {

  int c; char *fname_snames = NULL;
  int head = -1, tail = -1;
  while ((c = getopt(argc, argv, "l:H:T:h"))>=0) {
    switch (c) {
    case 'l': fname_snames = strdup(optarg); break;
    case 'H': head = atoi(optarg); break;
    case 'T': tail = atoi(optarg); break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (argc < optind + 2) { 
    usage(); 
    wzfatal("Please supply input file and output file.\n"); 
  }

  // input
  cgfile_t cgf = open_cgfile(argv[optind]);
  char *fname_index = get_fname_index(argv[optind]);
  index_t *idx = loadIndex(fname_index);
  free(fname_index);
  optind++;

  // output
  char *fname_out = strdup(argv[optind++]);

  snames_t snames = {0};
  if (optind + 1 < argc) {      // The requested sample names from command line
    for(int i = optind + 1; i < argc; ++i) {
      snames.array = realloc(snames.array, (snames.n+1));
      snames.array[snames.n++] = strdup(argv[i]);
    }
  } else {                      // from a file list
    snames = loadSampleNames(fname_snames, 1);
  }

  // check if we have index
  if (!idx) {
    fprintf(stderr, "Error, the cg file needs indexing for subsetting.\n");
    fflush(stderr);
    exit(1);
  }

  // read in the cgs
  if (snames.n == 0) {
    int npairs = 0;
    index_pair_t *pairs = index_pairs(idx, &npairs);
    if (tail > 0) {
      if (tail > npairs) tail = npairs;
      for (int i=0; i<tail; ++i) {
        snames.array = realloc(snames.array, (snames.n+1)*sizeof(const char*));
        snames.array[snames.n++] = strdup(pairs[npairs-tail+i].key);
      }
    } else {
      if (head < 1) head = 1;     // default to head 1
      if (head > npairs) head = npairs;
      for (int i=0; i<head; ++i) {
        snames.array = realloc(snames.array, (snames.n+1)*sizeof(const char*));
        snames.array[snames.n++] = strdup(pairs[i].key);
      }
    }
    free(pairs);
  }

  // output
  BGZF *fp = bgzf_open2(fname_out, "w");
  if (fp == NULL) {
    fprintf(stderr, "Error opening file for writing: %s\n", fname_out);
    exit(1);
  }
  cgdata_t cg = {0};
  for (int i=0; i<snames.n; ++i) {
    int64_t index = getIndex(idx, snames.array[i]);
    assert(index >= 0);
    assert(bgzf_seek(cgf.fh, index, SEEK_SET) == 0);
    read_cg2(&cgf, &cg);
    if (cg.n <= 0) break;
    cgdata_write1(fp, &cg);
  }
  bgzf_close(fp);               // output

  // output index
  cgfile_t cgf2 = open_cgfile(fname_out);
  index_t *idx2 = kh_init(index);
  int64_t addr = bgzf_tell(cgf2.fh);
  for (int i=0; i< snames.n; ++i) {
    if (!read_cg2(&cgf2, &cg)) {
      fprintf(stderr, "[Error] Data is shorter than the sample name list.\n");
      fflush(stderr);
      exit(1);
    }
    insert_index(idx2, snames.array[i], addr);
    addr = bgzf_tell(cgf2.fh);
  }
  free(cg.s);

  char *fname_index2 = get_fname_index(fname_out);
  FILE *out = fopen(fname_index2, "w");
  writeIndex(out, idx2);
  fclose(out);
  free(fname_index2);

  // clean up
  free(fname_out);
  bgzf_close(cgf.fh);
  bgzf_close(cgf2.fh);

  cleanIndex(idx);
  freeIndex(idx2);
  cleanSampleNames(&snames);
  
  return 0;
}
