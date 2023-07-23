#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdio.h>
#include "cgfile.h"
#include "vector.h"
#include "snames.h"

static int usage() {
  fprintf(stderr, "\nUsage: kycg unpack [options] <in.cg> [[sample 1], [sample 2], ...]\n\n");

  fprintf(stderr, "Options:\n");

  fprintf(stderr, "    -a        Process all samples\n");
  fprintf(stderr, "    -l        Path to the sample list. Ignored if sample names are provided on the command line.\n");
  fprintf(stderr, "    -H [N]    Process N samples from the start of the list, where N is less than or equal to the total number of samples.\n");
  fprintf(stderr, "    -T [N]    Process N samples from the end of the list, where N is less than or equal to the total number of samples. Requires index.\n");
  fprintf(stderr, "    -f [N]    Display format. Options are:\n");
  fprintf(stderr, "                   N == 0: Compound MU\n");
  fprintf(stderr, "                   N <  0: M<tab>U\n");
  fprintf(stderr, "                   N >  0: Fraction (with number for the min coverage)\n");
  fprintf(stderr, "    -c        Enable chunk process\n");
  fprintf(stderr, "    -s        Specify chunk size (default is 1M)\n");
  fprintf(stderr, "    -h        Display this help message\n\n");

  return 1;
}


static void print_cg1(cgdata_t *cg, uint64_t i, int printfmt3) {
  switch (cg->fmt) {
  case '0': {
      fputc(((cg->s[i>>3]>>(i&0x7))&0x1)+'0', stdout);
      break;
  }
  case '1': {
    fputc(cg->s[i], stdout);
    break;
  }
  case '2': {
    if (!cg->aux) fmt2_set_keys(cg);
    keys_t keys = *((keys_t*) cg->aux);
    uint64_t *data = (uint64_t*) fmt2_get_data(cg);
    assert(data[i] < keys.n);
    fprintf(stdout, "%s", keys.s[data[i]]);
    break;
  }
  case '3': {
    uint64_t *s = (uint64_t*) cg->s;
    if (printfmt3 == 0)
      fprintf(stdout, "%"PRIu32"", compressMU32(s[i]>>32, s[i]<<32>>32));
    else if (printfmt3 < 0)
      fprintf(stdout, "%"PRIu64"\t%"PRIu64"",s[i]>>32,s[i]<<32>>32);
    else {
      uint64_t M = s[i]>>32;
      uint64_t U = s[i]<<32>>32;
      if ((M==0 && U==0) || (M+U) < (uint32_t) printfmt3) fputs("NA", stdout);
      else fprintf(stdout, "%1.3f", (double) M/(M+U));
    }
      
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
    if (cg->s[i] == 2) {
      fputs("NA", stdout);
    } else {
      fputc(cg->s[i]+'0', stdout);
    }
    break;
  }
  case '6': {
    uint64_t *s = (uint64_t*) cg->s;
    fprintf(stdout, "%"PRIu64, s[i]);
    break;
  }
  default: usage(); wzfatal("Unrecognized format: %c.\n", cg->fmt);
  }
}

static void print_cgs_chunk(cgdata_v *cgs, uint64_t s, int printfmt3) {
  uint64_t i,m, k, kn = cgs->size;
  cgdata_t expanded = {0};
  decompress(ref_cgdata_v(cgs, 0), &expanded);
  uint64_t n = expanded.n;
  cgdata_t *sliced = calloc(kn, sizeof(cgdata_t));
  for (m=0; m <= n/s; ++m) {
    for (k=0; k<kn; ++k) {
      decompress(ref_cgdata_v(cgs, k), &expanded);
      slice(&expanded, m*s, (m+1)*s-1, &sliced[k]);
    }
    for (i=0; i<sliced[0].n; ++i) {
      for (k=0; k<kn; ++k) {
        if(k) fputc('\t', stdout);
        print_cg1(&sliced[k],i,printfmt3);
      }
      fputc('\n', stdout);
    }
  }

  for (k=0; k<kn; ++k) free(sliced[k].s);
  free(expanded.s); free(sliced);
}

static void print_cgs(cgdata_v *cgs, int printfmt3) {
  uint64_t i, k, kn = cgs->size;
  cgdata_t *expanded = calloc(kn, sizeof(cgdata_t));
  for (k=0; k<kn; ++k) {
    decompress(ref_cgdata_v(cgs, k), expanded+k);
  }
  for (i=0; i<expanded[0].n; ++i) {
    for (k=0; k<kn; ++k) {
      if(k) fputc('\t', stdout);
      print_cg1(expanded+k, i, printfmt3);
    }
    fputc('\n', stdout);
  }
  for (k=0; k<kn; ++k) free_cgdata(&expanded[k]);
  free(expanded);
}

int main_unpack(int argc, char *argv[]) {

  int c, read_all = 0, chunk = 0;
  int printfmt3 = 0;
  uint64_t chunk_size = 1000000; char *fname_snames = NULL;
  int head = -1, tail = -1;
  while ((c = getopt(argc, argv, "cs:H:T:f:ah"))>=0) {
    switch (c) {
    case 'c': chunk = 1; break;
    case 's': chunk_size = atoi(optarg); break;
    case 'l': fname_snames = strdup(optarg); break;
    case 'H': head = atoi(optarg); break;
    case 'T': tail = atoi(optarg); break;
    case 'a': read_all = 1; break;
    case 'f': printfmt3 = atoi(optarg); break;
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
  if (fname_index) free(fname_index);

  snames_t snames = {0};
  if (optind + 1 < argc) {      // The requested sample names from command line
    for(int i = optind + 1; i < argc; ++i) {
      snames.s = realloc(snames.s, (snames.n+1));
      snames.s[snames.n++] = strdup(argv[i]);
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
  cgdata_v *cgs = NULL;
  if (snames.n > 0) {
    cgs = read_cgs_with_snames(&cgf, idx, &snames);
  } else if (read_all) {
    cgs = read_cgs_all(&cgf);
  } else if (head > 0) {
    cgs = read_cgs_from_head(&cgf, head);
  } else if (tail > 0) {
    cgs = read_cgs_from_tail(&cgf, idx, tail);
  } else {
    cgs = read_cgs_from_head(&cgf, 1);
  }

  // output the cgs
  if (chunk) print_cgs_chunk(cgs, chunk_size, printfmt3);
  else print_cgs(cgs, printfmt3);

  // clean up
  for (uint64_t i=0; i<cgs->size; ++i) free_cgdata(ref_cgdata_v(cgs,i));
  free_cgdata_v(cgs);
  bgzf_close(cgf.fh);
  if (idx) cleanIndex(idx);
  cleanSampleNames(&snames);
  
  return 0;
}
