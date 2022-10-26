#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdio.h>
#include "kycg.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: kycg pack [options] <in.bed> <out.cg>\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -f        format (0,1,2,3,4,a)\n");
  fprintf(stderr, "    -v        verbose\n");
  fprintf(stderr, "    -h        This help\n");
  fprintf(stderr, "\n");

  return 1;
}

int main_unpack(int argc, char *argv[]) {

  int c; int verbose = 0;
  while ((c = getopt(argc, argv, "vh"))>=0) {
    switch (c) {
    case 'v': verbose = 1; break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 1 > argc) { 
    usage(); 
    wzfatal("Please supply input file.\n"); 
  }
  
  cgdata_t *cg = read_cg(argv[optind]);
  switch (cg->fmt) {
  case '0': {
    uint64_t i;
    for(i=0; i<cg->n; ++i) {
      fputc(((cg->s[i>>3]>>(i&0x7))&0x1)+'0', stdout); fputc('\n', stdout);
    }
    break;
  }
  case '1': {
    cgdata_t *cg2 = fmt1_decompress(cg);
    uint64_t i;
    for(i=0; i<cg2->n; ++i) {
      fputc(cg2->s[i], stdout); fputc('\n', stdout);
    }
    free_cgdata(cg2);
    break;
  }
  case '3': {
    cgdata_t *cg2 = fmt3_decompress(cg);
    uint64_t i; uint64_t *s = (uint64_t*) cg2->s;
    for(i=0; i<cg2->n; ++i)
      fprintf(stdout, "%"PRIu64"\t%"PRIu64"\n",s[i]>>32,s[i]<<32>>32);
    free_cgdata(cg2);
    break;
  }
  case '4': {
    cgdata_t *cg2 = fmt4_decompress(cg);
    uint64_t i;
    for(i=0; i<cg2->n; ++i) {
      if (((float_t*) (cg2->s))[i]<0) {
        fputs("NA\n", stdout);
      } else {
        fprintf(stdout, "%1.3f\n", ((float_t*) (cg2->s))[i]);
      }
    }
    free_cgdata(cg2);
    break;
  }
  case '5': {
    cgdata_t *cg2 = fmt5_decompress(cg);
    uint64_t i;
    for(i=0; i<cg2->n; ++i) {
      fputc(cg2->s[i]+'0', stdout); fputc('\n', stdout);
    }
    free_cgdata(cg2);
    break;
  }
  default: usage(); wzfatal("Unrecognized format: %c.\n", cg->fmt);
  }

  free_cgdata(cg);
  return 0;
}
