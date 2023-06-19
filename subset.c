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

  int c, verbose = 0; char *fname_snames = NULL;
  while ((c = getopt(argc, argv, "l:H:T:hv"))>=0) {
    switch (c) {
    case 'v': verbose = 1; break;
    case 'l': fname_snames = strdup(optarg); break;
    case 'H': head = atoi(optarg); break;
    case 'T': tail = atoi(optarg); break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 1 > argc) { 
    usage(); 
    wzfatal("Please supply input file.\n"); 
  }

  cgfile_t cgf = open_cgfile(argv[optind]);
  char *fname_index = get_fname_index(argv[optind]);
  index_t *idx = loadIndex(fname_index);

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
  if (!idx && (snames.n > 0 || tail > 0)) {
    fprintf(stderr, "Error, the cg file needs indexing for random sample access.\n");
    fflush(stderr);
    exit(1);
  }

  // read in the cgs
  if (snames.n > 0) {
    cgdata_t cg = {0};
    for (int i=0; i<snames->n; ++i) {
      int64_t index = getIndex(idx, snames->array[i]);
      assert(index >= 0);
      assert(bgzf_seek(cgf.fh, index, SEEK_SET) == 0);
      read_cg2(cgf, &cg);
      if (cg.n <= 0) break;
      cgdata_write1(cgf.fh, &cg);
    }
    free(cg.s);
  } else if (head > 0) {
    cgdata_t cg = {0};
    for (int i=0; i<head; ++i) {
      read_cg2(cgf, &cg);
      if (cg.n <= 0) break;
      cgdata_write1(cgf.fh, &cg);
    }
    free(cg.s);
  } else if (tail > 0) {
    int npairs = 0;
    index_pair_t *pairs = index_pairs(idx, &npairs);
    if (tail > npairs) tail = npairs;
    cgdata_t cg = {0};
    for (int i=0; i<tail; ++i) {
      int64_t index = pairs[npairs-tail+i].value;
      assert(index >= 0);
      assert(bgzf_seek(cgf.fh, index, SEEK_SET) == 0);
      read_cg2(cgf, &cg);
      if (cg.n <= 0) break;
      cgdata_write1(cgf.fh, &cg);
    }
    free(cg.s);
  } else {
    cgdata_t cg = {0};
    read_cg2(cgf, &cg);
    if (cg.n <= 0) break;
    cgdata_write1(cgf.fh, &cg);
    free(cg.s);
  }

  // clean up
  bgzf_close(cgf.fh);
  if (idx) destroyIndex(idx);
  cleanSampleNames(&snames);
  
  return 0;
}
