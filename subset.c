#include "cgfile.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: kycg subset [options] <in.cg> [<sample1>, ...]\n");
  fprintf(stderr, "If -o <out.cg>, an index will also be generated. Otherwise, output .cg to stdout without index.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -v        verbose\n");
  fprintf(stderr, "    -o        output cg file name. if missing, output to stdout without index.\n");
  fprintf(stderr, "    -l        Path to the sample list. Ignored if sample names are provided on the command line.\n");
  fprintf(stderr, "    -H [N]    Process N samples from the start of the list, where N is less than or equal to the total number of samples.\n");
  fprintf(stderr, "    -T [N]    Process N samples from the end of the list, where N is less than or equal to the total number of samples. Requires index.\n");
  fprintf(stderr, "    -h        This help\n");
  fprintf(stderr, "\n");

  return 1;
}

int main_subset(int argc, char *argv[]) {

  int c; char *fname_snames = NULL;
  int head = -1, tail = -1;
  char *fname_out = NULL;
  while ((c = getopt(argc, argv, "o:l:H:T:h"))>=0) {
    switch (c) {
    case 'o': fname_out = strdup(optarg); break;
    case 'l': fname_snames = strdup(optarg); break;
    case 'H': head = atoi(optarg); break;
    case 'T': tail = atoi(optarg); break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (argc < optind + 1) { 
    usage(); 
    wzfatal("Please supply input file and output file.\n"); 
  }

  // input
  cgfile_t cgf = open_cgfile(argv[optind]);
  char *fname_index = get_fname_index(argv[optind]);
  index_t *idx = loadIndex(fname_index);
  free(fname_index);
  optind++;

  // get sample names
  snames_t snames = {0};
  if (optind < argc) {      // sample names from command line
    for(int i = optind; i < argc; ++i) {
      snames.s = realloc(snames.s, (snames.n+1));
      snames.s[snames.n++] = strdup(argv[i]);
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
      if (head < 1) head = 1;     // default to head 1
      if (head > npairs) head = npairs;
      for (int i=0; i<head; ++i) {
        snames.s = realloc(snames.s, (snames.n+1)*sizeof(const char*));
        snames.s[snames.n++] = strdup(pairs[i].key);
      }
    }
    free(pairs);
  }

  // output
  BGZF *fp;
  if (fname_out) fp = bgzf_open2(fname_out, "w");
  else fp = bgzf_dopen(fileno(stdout), "w");
  
  if (fp == NULL) {
    fprintf(stderr, "Error opening file for writing: %s\n", fname_out);
    exit(1);
  }
  cgdata_t cg = {0};
  for (int i=0; i<snames.n; ++i) {
    int64_t index = getIndex(idx, snames.s[i]);
    assert(index >= 0);
    assert(bgzf_seek(cgf.fh, index, SEEK_SET) == 0);
    read_cg2(&cgf, &cg);
    if (cg.n <= 0) break;
    cgdata_write1(fp, &cg);
  }
  bgzf_close(fp);               // output

  if (fname_out) {
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
      insert_index(idx2, snames.s[i], addr);
      addr = bgzf_tell(cgf2.fh);
    }
    free(cg.s);

    char *fname_index2 = get_fname_index(fname_out);
    FILE *out = fopen(fname_index2, "w");
    writeIndex(out, idx2);
    fclose(out);
    free(fname_index2);
    free(fname_out);
    bgzf_close(cgf2.fh);
    freeIndex(idx2);
  }
  
  // clean up
  bgzf_close(cgf.fh);
  cleanIndex(idx);
  cleanSampleNames(&snames);
  
  return 0;
}
