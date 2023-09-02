#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdio.h>
#include "cfile.h"
#include "vector.h"
#include "snames.h"

static int usage() {
  fprintf(stderr, "\nUsage: yame unpack [options] <in.cx> [[sample 1], [sample 2], ...]\n\n");

  fprintf(stderr, "Options:\n");

  fprintf(stderr, "    -a        Process all samples\n");
  fprintf(stderr, "    -C        Output column names\n");
  fprintf(stderr, "    -l        Path to the sample list. Ignored if sample names are provided on the command line.\n");
  fprintf(stderr, "    -H [N]    Process N samples from the start of the list, where N is less than or equal to the total number of samples.\n");
  fprintf(stderr, "    -T [N]    Process N samples from the end of the list, where N is less than or equal to the total number of samples. Requires index.\n");
  fprintf(stderr, "    -f [N]    Display format. Options are:\n");
  fprintf(stderr, "                   N == 0: Compound MU\n");
  fprintf(stderr, "                   N <  0: M<tab>U\n");
  fprintf(stderr, "                   N >  0: Fraction (with number for the min coverage)\n");
  fprintf(stderr, "    -c        Enable chunk process\n");
  fprintf(stderr, "    -u [int]  number of bytes for each unit data while inflated. Lower number is more memory efficient but could be lossier. Can only be 1-8. 0 means this will be inferred from data.\n");
  fprintf(stderr, "    -s        Specify chunk size (default is 1M)\n");
  fprintf(stderr, "    -h        Display this help message\n\n");

  return 1;
}

static void print_cdata1(cdata_t *c, uint64_t i, int f3_fmt) {
  switch (c->fmt) {
  case '0': {
      fputc(((c->s[i>>3]>>(i&0x7))&0x1)+'0', stdout);
      break;
  }
  case '1': {
    fputc(c->s[i], stdout);
    break;
  }
  case '2': {
    fprintf(stdout, "%s", f2_get_string(c, i));
    break;
  }
  case '3': {
    uint64_t mu = f3_get_mu(c, i);
    if (f3_fmt == 0)
      fprintf(stdout, "%"PRIu64"", mu);
    else if (f3_fmt < 0)
      fprintf(stdout, "%"PRIu64"\t%"PRIu64"",mu>>32, mu<<32>>32);
    else {
      uint64_t M = mu>>32;
      uint64_t U = mu<<32>>32;
      if ((M==0 && U==0) || (M+U) < (uint64_t) f3_fmt) fputs("NA", stdout);
      else fprintf(stdout, "%1.3f", (double) M/(M+U));
    }
    break;
  }
  case '4': {
    if (((float_t*) (c->s))[i]<0) {
      fputs("NA", stdout);
    } else {
      fprintf(stdout, "%1.3f", ((float_t*) (c->s))[i]);
    }
    break;
  }
  case '5': {
    if (c->s[i] == 2) {
      fputs("NA", stdout);
    } else {
      fputc(c->s[i]+'0', stdout);
    }
    break;
  }
  case '6': {
    uint64_t *s = (uint64_t*) c->s;
    fprintf(stdout, "%"PRIu64, s[i]);
    break;
  }
  default: usage(); wzfatal("Unrecognized format: %c.\n", c->fmt);
  }
}

static void print_cdata_chunk(cdata_v *cs, uint64_t s, int f3_fmt) {
  uint64_t i,m, k, kn = cs->size;
  cdata_t expanded = {0};
  decompress(ref_cdata_v(cs, 0), &expanded);
  uint64_t n = expanded.n;
  cdata_t *sliced = calloc(kn, sizeof(cdata_t));
  for (m=0; m <= n/s; ++m) {
    for (k=0; k<kn; ++k) {
      decompress(ref_cdata_v(cs, k), &expanded);
      slice(&expanded, m*s, (m+1)*s-1, &sliced[k]);
    }
    for (i=0; i<sliced[0].n; ++i) {
      for (k=0; k<kn; ++k) {
        if(k) fputc('\t', stdout);
        print_cdata1(&sliced[k],i,f3_fmt);
      }
      fputc('\n', stdout);
    }
  }

  for (k=0; k<kn; ++k) free(sliced[k].s);
  free(expanded.s); free(sliced);
}

static void print_cdata(cdata_v *cs, int f3_fmt) {
  uint64_t i, k, kn = cs->size;
  cdata_t *expanded = calloc(kn, sizeof(cdata_t));
  for (k=0; k<kn; ++k) {
    decompress(ref_cdata_v(cs, k), expanded+k);
  }
  for (i=0; i<expanded[0].n; ++i) {
    for (k=0; k<kn; ++k) {
      if(k) fputc('\t', stdout);
      print_cdata1(expanded+k, i, f3_fmt);
    }
    fputc('\n', stdout);
  }
  for (k=0; k<kn; ++k) free_cdata(&expanded[k]);
  free(expanded);
}

int main_unpack(int argc, char *argv[]) {

  int c, read_all = 0, chunk = 0;
  int f3_fmt = 0;
  uint64_t chunk_size = 1000000; char *fname_snames = NULL;
  int head = -1, tail = -1;
  uint8_t unit = 0; // default: auto-inferred
  int print_column_names = 0;
  while ((c = getopt(argc, argv, "cs:l:H:T:f:u:Cah"))>=0) {
    switch (c) {
    case 'c': chunk = 1; break;
    case 's': chunk_size = atoi(optarg); break;
    case 'l': fname_snames = strdup(optarg); break;
    case 'H': head = atoi(optarg); break;
    case 'T': tail = atoi(optarg); break;
    case 'u': unit = atoi(optarg); break;
    case 'C': print_column_names = 1; break;
    case 'a': read_all = 1; break;
    case 'f': f3_fmt = atoi(optarg); break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 1 > argc) { 
    usage(); 
    wzfatal("Please supply input file.\n"); 
  }

  char *fname_in = strdup(argv[optind]);
  cfile_t cf = open_cfile(fname_in);
  char *fname_index = get_fname_index(fname_in);
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
  if ((tail > 0 && !idx) || (snames.n > 0 && strcmp(fname_in, "-") != 0 && !idx)) {
    fprintf(stderr, "Error, the cx file needs indexing for random sample access.\n");
    fflush(stderr);
    exit(1);
  }
  
  // read in the cdata
  cdata_v *cs = NULL;
  if (idx && snames.n > 0) {
    cs = read_cdata_with_snames(&cf, idx, &snames);
  } else if (read_all) {
    cs = read_cdata_all(&cf);
  } else if (head > 0) {
    cs = read_cdata_from_head(&cf, head);
  } else if (tail > 0) {
    cs = read_cdata_from_tail(&cf, idx, tail);
  } else {
    cs = read_cdata_from_head(&cf, 1);
  }

  // set unit size
  if (unit > 8 || ((unit&0x1) && unit != 1 && unit)) {
    fprintf(stderr, "[%s:%d] Unit size (%u) can only be 1,2,4,6,8.\n", __func__, __LINE__, unit);
    fflush(stderr);
    exit(1);
  }
  for (uint64_t i=0; i<cs->size; ++i) ref_cdata_v(cs, i)->unit = unit;

  // output headers
  if (print_column_names) {
    if (!snames.n) {
      fprintf(stderr, "[%s:%d] Error, index file is missing for printing sample names.\n", __func__, __LINE__);
      fflush(stderr);
      exit(1);
    }
    for (int i=0; i<snames.n; ++i) {
      if (i) fputc('\t', stdout);
      fputs(snames.s[i], stdout);
    }
    fputc('\n', stdout);
  }

  // output the cs
  if (chunk) print_cdata_chunk(cs, chunk_size, f3_fmt);
  else print_cdata(cs, f3_fmt);

  // clean up
  for (uint64_t i=0; i<cs->size; ++i) free_cdata(ref_cdata_v(cs,i));
  free_cdata_v(cs);
  bgzf_close(cf.fh);
  free(fname_in);
  if (idx) cleanIndex(idx);
  cleanSampleNames(&snames);
  
  return 0;
}
