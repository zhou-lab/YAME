#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdio.h>
#include "wzmisc.h"
#include "wzbed.h"
#include "kycg.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: kycg overlap [options] <query.cg> <feature.cg>\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -u        Universe .cg file (optional).\n");
  fprintf(stderr, "    -h        This help.\n");
  fprintf(stderr, "Returns four numbers, nu, nf, nq, nfq;\n");
  fprintf(stderr, "\n");

  return 1;
}

static void bit_mask(uint8_t *s, uint8_t *mask, size_t n) {
  size_t i;
  for (i=0; i<(n>>3)+1; ++i) s[i] &= mask[i];
}

static size_t bit_count(cgdata_t cg) {

  /* create a look-up table */
  int byte2cnt[256]; int p;
  for (p=0; p<256; ++p) {
    unsigned char q = p; int ii, cnt = 0;
    for (ii=0; ii<8; ++ii) { if (q&1) cnt++; q>>=1; }
    byte2cnt[p] = cnt;
  }
  
  size_t i,k,m = 0;
  for (i=0; i<(cg.n>>3); ++i) m += byte2cnt[cg.s[i]];
  for (k=0; k<(cg.n&0x7); ++k) m += (cg.s[i]>>k) & 0x1;
  return m;
}

/* static unsigned char* read_vec(FILE *fh, int64_t *n) { */

/*   uint64_t ncols; wzfread(&ncols, sizeof(uint64_t), 1, fh); */
/*   uint8_t fmt; unsigned char *s = NULL; */
/*   wzfread(&fmt, 1, 1, fh); */
/*   if (fmt == 0) {               /\* bit vector *\/ */
/*     wzfread(n, sizeof(int64_t), 1, fh); */
/*     s = malloc(((*n)>>3)+1); */
/*     wzfread(s, 1, ((*n)>>3)+1, fh); */
/*     return s; */
    
/*   } else if (fmt == 1) {        /\* run length encoding *\/ */
/*     int64_t n_rle; wzfread(&n_rle, sizeof(int64_t), 1, fh); */
/*     unsigned char *s_rle = malloc(n_rle); */
/*     wzfread(s_rle, 1, n_rle, fh); */
/*     int i; */
/*     *n=0; */
/*     for (i=0; i<n_rle/3; ++i) */
/*       (*n) += *((uint16_t*) (s_rle+i*3+1)); */
/*     s = calloc(((*n)>>3)+1, 1); size_t sum; */
/*     for (i=0, sum=0; i<n_rle/3; ++i) { */
/*       uint16_t l = *((uint16_t*) (s_rle+i*3+1)); */
/*       if (s_rle[i*3] == 1) { */
/*         size_t j; */
/*         for(j=sum; j<sum+l; ++j) { */
/*           s[j>>3] |= (1<<(j&0x7)); */
/*         } */
/*       } */
/*       sum += l; */
/*     } */
/*     free(s_rle); */

/*     return s; */
    
/*   } else { */
/*     wzfatal("Unknown format code: %d", fmt); */
/*   } */
/*   return s; */
/* } */

static cgdata_t read_cg2fmt0(cgfile_t *cgf) {
  cgdata_t cg = read_cg(cgf);
  switch (cg.fmt) {
  case '0': break;
  case '1': {
    /* int64_t n_rle; wzfread(&n_rle, sizeof(int64_t), 1, fh); */
    /* unsigned char *s_rle = malloc(n_rle); */
    /* wzfread(s_rle, 1, n_rle, fh); */
    uint8_t *s_rle = cg.s; uint64_t n_rle = cg.n;
    uint64_t i;
    cg.n=0;
    for (i=0; i<n_rle/3; ++i) cg.n += *((uint16_t*) (s_rle+i*3+1));
    cg.s = calloc((cg.n>>3)+1, 1); size_t sum;
    for (i=0, sum=0; i<n_rle/3; ++i) {
      uint16_t l = *((uint16_t*) (s_rle+i*3+1));
      if (s_rle[i*3] == '1') {
        size_t j;
        for(j=sum; j<sum+l; ++j) {
          cg.s[j>>3] |= (1<<(j&0x7));
        }
      }
      sum += l;
    }
    free(s_rle);
    break;
  }
  default: wzfatal("Format %c unsupported.\n", cg.fmt);
  }
  return cg;
}

/* The design, first 10 bytes are uint64_t (length) + uint16_t (0=vec; 1=rle) */
int main_overlap(int argc, char *argv[]) {
  int c; char *upath = NULL;
  while ((c = getopt(argc, argv, "u:h"))>=0) {
    switch (c) {
    case 'u': upath = strdup(optarg); break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 2 > argc) { 
    usage(); 
    wzfatal("Please supply input and output file.\n"); 
  }

  cgfile_t cgf_fea = open_cgfile(argv[optind++]);
  cgfile_t cgf_qry = open_cgfile(argv[optind++]);
  cgdata_t cg_fea = read_cg2fmt0(&cgf_fea);
  cgdata_t cg_qry = read_cg2fmt0(&cgf_qry);
  if (cg_fea.n != cg_qry.n) {
    wzfatal("Feature (%"PRId64") and query (%"PRId64") has different lengths.", cg_fea.n, cg_qry.n);
  }
  bgzf_close(cgf_fea.fh);
  bgzf_close(cgf_qry.fh);

  size_t m_uni = 0;
  if (upath != NULL) {
    cgfile_t cgf_uni = open_cgfile(upath);
    cgdata_t cg_uni = read_cg2fmt0(&cgf_uni);
    if (cg_uni.n != cg_fea.n || cg_uni.n != cg_qry.n) {
      wzfatal("Universe has different length.");
    }
    bit_mask(cg_fea.s, cg_uni.s, cg_uni.n);
    bit_mask(cg_qry.s, cg_uni.s, cg_uni.n);
    m_uni = bit_count(cg_uni);
    free(cg_uni.s); free(upath);
    bgzf_close(cgf_uni.fh);
  } else { m_uni = cg_fea.n; }
  
  /* FILE *fh_fea = fopen(argv[optind++], "rb"); int64_t n_fea; */
  /* FILE *fh_qry = fopen(argv[optind++], "rb"); int64_t n_qry; */
  /* unsigned char *s_fea = read_vec(fh_fea, &n_fea); */
  /* unsigned char *s_qry = read_vec(fh_qry, &n_qry); */
  /* fclose(fh_fea); */
  /* fclose(fh_qry); */
  /* if (n_fea != n_qry) { */
  /*   wzfatal("Feature (%"PRId64") and query (%"PRId64") has different lengths.", n_fea, n_qry); } */

  /* size_t m_uni = 0; */
  /* if (upath != NULL) { */
  /*   FILE *fh_uni = fopen(upath, "rb"); int64_t n_uni; */
  /*   uint8_t *s_uni = read_vec(fh_uni, &n_uni); */
  /*   if (n_uni != n_fea || n_uni != n_qry) { wzfatal("Universe has different length."); } */
  /*   bit_mask(s_fea, s_uni, n_uni); */
  /*   bit_mask(s_qry, s_uni, n_uni); */
  /*   m_uni = bit_count(s_uni, n_uni); */
  /*   free(s_uni); */
  /*   fclose(fh_uni); free(upath); */
  /* } else { m_uni = n_fea; } */

  size_t nf = bit_count(cg_fea);
  size_t nq = bit_count(cg_qry);
  bit_mask(cg_fea.s, cg_qry.s, cg_fea.n);
  size_t nfq = bit_count(cg_fea);
  
  /* int64_t i; int64_t nfq=0, nf=0, nq=0; */
  /* for (i=0; i<(n_fea>>3); ++i) { */
  /*   nfq += byte2cnt[s_fea[i] & s_qry[i]]; */
  /*   nf  += byte2cnt[s_fea[i]]; */
  /*   nq  += byte2cnt[s_qry[i]]; */
  /* } */
  /* size_t k; */
  /* for (k=0; k<(n_fea&0x7); ++k) { */
  /*   nfq += ((s_fea[i] & s_qry[i])>>k)&0x1; */
  /*   nf  += ((s_fea[i])>>k)&0x1; */
  /*   nq  += ((s_qry[i])>>k)&0x1; */
  /* } */
  
  fprintf(stdout, "%zu\t%zu\t%zu\t%zu\n", m_uni, nf, nq, nfq);
  free(cg_fea.s);
  free(cg_qry.s);
  
  /* free(s_fea); */
  /* free(s_qry); */
  return 0;
}
