#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdio.h>
#include "kycg.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: kycg unpack [options] <in.cg>\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -v        verbose\n");
  fprintf(stderr, "    -h        This help\n");
  fprintf(stderr, "\n");

  return 1;
}

static void print_cg1(cgdata_t *cg, uint64_t i) {
  switch (cg->fmt) {
  case '0': {
      fputc(((cg->s[i>>3]>>(i&0x7))&0x1)+'0', stdout);
      break;
  }
  case '1': {
    fputc(cg->s[i], stdout);
    break;
  }
  case '3': {
    uint64_t *s = (uint64_t*) cg->s;
    fprintf(stdout, "%"PRIu64"\t%"PRIu64"",s[i]>>32,s[i]<<32>>32);
    break;
  }
  case '4': {
    if (((float_t*) (cg->s))[i]<0) {
      fputs("NA", stdout);
    } else {
      fprintf(stdout, "%1.3f", ((float_t*) (cg->s))[i]);
    }
    break;
  }
  case '5': {
    fputc(cg->s[i]+'0', stdout);
    break;
  }
  default: usage(); wzfatal("Unrecognized format: %c.\n", cg->fmt);
  }
}

static void print_cg(cgdata_t *cg) {

  cgdata_t cg2 = decompress(cg);
  uint64_t i;
  for (i=0; i<cg2.n; ++i) {
    print_cg1(&cg2, i); fputc('\n', stdout);
  }
  free(cg2.s);
}

static void print_cgs(cgdata_v *cgs) {
  uint64_t i,k; cgdata_v *cgs_d = init_cgdata_v(cgs->size);
  for (k=0; k<cgs->size; ++k)
    *next_ref_cgdata_v(cgs_d) = decompress(ref_cgdata_v(cgs, k));
  for (i=0; i<ref_cgdata_v(cgs_d,0)->n; ++i) {
    for (k=0; k<cgs->size; ++k) {
      if(k) fputc('\t', stdout);
      print_cg1(ref_cgdata_v(cgs_d,k),i);
    }
    fputc('\n', stdout);
  }
}

int main_unpack(int argc, char *argv[]) {

  int c, verbose = 0, read_all = 0;
  while ((c = getopt(argc, argv, "avh"))>=0) {
    switch (c) {
    case 'v': verbose = 1; break;
    case 'a': read_all = 1; break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 1 > argc) { 
    usage(); 
    wzfatal("Please supply input file.\n"); 
  }

  cgfile_t cgf = open_cgfile(argv[optind]);
  if (read_all) {
    cgdata_v *cgs = read_cg_all(&cgf);
    print_cgs(cgs);
    uint64_t i;
    for (i=0; i<cgs->size; ++i) free(ref_cgdata_v(cgs,i)->s);
    free_cgdata_v(cgs);
  } else {
    cgdata_t cg = read_cg(&cgf);
    print_cg(&cg);
    free(cg.s);
  }
  gzclose(cgf.fh);
  
  return 0;
}
