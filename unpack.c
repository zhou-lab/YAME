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

  cgdata_t expanded = {0};
  decompress(cg, &expanded);
  uint64_t i;
  for (i=0; i<expanded.n; ++i) {
    print_cg1(&expanded, i); fputc('\n', stdout);
  }
  free(expanded.s);
}

static void print_cgs_chunk(cgdata_v *cgs, uint64_t s) {
  uint64_t i,k,m;
  cgdata_t expanded = {0};
  decompress(ref_cgdata_v(cgs, 0), &expanded);
  uint64_t n = expanded.n;
  cgdata_t *sliced = calloc(cgs->size, sizeof(cgdata_t));
  for (m=0; m <= n/s; ++m) {
    for (k=0; k<cgs->size; ++k) {
      decompress(ref_cgdata_v(cgs, k), &expanded);
      slice(&expanded, m*s, (m+1)*s-1, &sliced[k]);
    }
    for (i=0; i<sliced[k].n; ++i) {
      for (k=0; k<cgs->size; ++k) {
        if(k) fputc('\t', stdout);
        print_cg1(&sliced[k],i);
      }
      fputc('\n', stdout);
    }
  }

  for (k=0; k<cgs->size; ++k) free(sliced[k].s);
  free(expanded.s); free(sliced);
}

static void print_cgs(cgdata_v *cgs) {
  uint64_t i,k; cgdata_v *cgs_d = init_cgdata_v(cgs->size);
  cgdata_t *expanded;
  for (k=0; k<cgs->size; ++k) {
    expanded = next_ref_cgdata_v(cgs_d);
    decompress(ref_cgdata_v(cgs, k), expanded);
  }
  for (i=0; i<ref_cgdata_v(cgs_d,0)->n; ++i) {
    for (k=0; k<cgs->size; ++k) {
      if(k) fputc('\t', stdout);
      print_cg1(ref_cgdata_v(cgs_d,k),i);
    }
    fputc('\n', stdout);
  }
}

int main_unpack(int argc, char *argv[]) {

  int c, verbose = 0, read_all = 0, chunk = 0;
  uint64_t chunk_size = 1000000;
  while ((c = getopt(argc, argv, "s:avh"))>=0) {
    switch (c) {
    case 'c': chunk = 1; break;
    case 's': chunk_size = atoi(optarg); break;
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
    if (chunk) print_cgs_chunk(cgs, chunk_size);
    else print_cgs(cgs);
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
