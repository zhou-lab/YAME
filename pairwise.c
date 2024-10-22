#include "cfile.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: yame pairwise [options] <MU1.cx> (<MU2.cx>)\n");
  fprintf(stderr, "Return a format 6 set that represent differential methylation between MU1 and MU2.\n");
  fprintf(stderr, "If MU2 is not given, use the top 2 samples in MU1.cx.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -o        output cx file name. if missing, output to stdout without index.\n");
  fprintf(stderr, "    -H        1: higher meth level in sample 1 than 2 (default).\n");
  fprintf(stderr, "              2: higher meth level in sample 2 than 1.\n");
  fprintf(stderr, "              others: diff levels, i.e., 1 and 2 combined.\n");
  fprintf(stderr, "    -c        minimum coverage (default: 1)\n");
  fprintf(stderr, "    -h        This help\n");
  fprintf(stderr, "\n");

  return 1;
}

int main_pairwise(int argc, char *argv[]) {

  int c; int direc = 1;
  int64_t min_coverage = 1; char *fname_out = NULL;
  while ((c = getopt(argc, argv, "o:c:H:h"))>=0) {
    switch (c) {
    case 'o': fname_out = strdup(optarg); break;
    case 'c': min_coverage = atoi(optarg); break;
    case 'H': direc = atoi(optarg); break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }
  if (min_coverage < 1) min_coverage = 1;

  if (optind + 1 > argc) { 
    usage(); 
    wzfatal("Please supply input file.\n");
  }

  cfile_t cf1 = open_cfile(argv[optind++]);
  cdata_t c1 = read_cdata1(&cf1);
  cdata_t c2 = {0};
  if (optind >= argc) {
    c2 = read_cdata1(&cf1);
  } else {
    cfile_t cf2 = open_cfile(argv[optind++]);
    c2 = read_cdata1(&cf2);
    bgzf_close(cf2.fh);
  }
  bgzf_close(cf1.fh);

  decompress2(&c1);
  decompress2(&c2);

  if (c1.n != c2.n) wzfatal("Two inputs have different dimensions: %"PRIu64" vs %"PRIu64"\n", c1.n, c2.n);

  cdata_t c_out = {.fmt = '6', .n = c1.n };
  c_out.s = calloc((c_out.n+3)/4, sizeof(uint8_t));
  for (uint64_t i=0; i<c1.n; ++i) {
    uint64_t mu1 = f3_get_mu(&c1, i);
    uint64_t mu2 = f3_get_mu(&c2, i);
    if (MU2cov(mu1) >= (uint64_t) min_coverage &&
        MU2cov(mu2) >= (uint64_t) min_coverage) {
      if (direc == 1) {
        if (MU2beta(mu1) > MU2beta(mu2)) FMT6_SET1(c_out, i);
        else FMT6_SET0(c_out, i);
      } else if (direc == 2) {
        if (MU2beta(mu1) < MU2beta(mu2)) FMT6_SET1(c_out, i);
        else FMT6_SET0(c_out, i);
      } else {
        if (MU2beta(mu1) != MU2beta(mu2)) FMT6_SET1(c_out, i);
        else FMT6_SET0(c_out, i);
      }
    }
  }
  free_cdata(&c1);
  free_cdata(&c2);
  
  cdata_compress(&c_out);
  BGZF *fp_out;
  if (fname_out) fp_out = bgzf_open2(fname_out, "w");
  else fp_out = bgzf_dopen(fileno(stdout), "w");
  if (fp_out == NULL) {
    fprintf(stderr, "Error opening file for writing: %s\n", fname_out);
    exit(1);
  }
  cdata_write1(fp_out, &c_out);
  bgzf_close(fp_out);
  free_cdata(&c_out);
  
  return 0;
}
