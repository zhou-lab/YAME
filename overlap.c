#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdio.h>
#include "wzmisc.h"
#include "wzbed.h"
#include "cgfile.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: kycg overlap [options] <feature.cg> <query.cg>\n");
  fprintf(stderr, "Query or feature can be a multi-sample set.\n");
  fprintf(stderr, "If feature.cg is unseekable, only first sample will be used.\n");
  fprintf(stderr, "If the input is format 3: test overlap of M+U.\n");
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

/* The design, first 10 bytes are uint64_t (length) + uint16_t (0=vec; 1=rle) */
int main_overlap(int argc, char *argv[]) {
  int c; char *upath = NULL; int print_header=0;
  while ((c = getopt(argc, argv, "u:Hh"))>=0) {
    switch (c) {
    case 'u': upath = strdup(optarg); break;
    case 'H': print_header = 1; break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 2 > argc) { 
    usage(); 
    wzfatal("Please supply input and output file.\n"); 
  }

  /* import universe if provided */
  cgdata_t cg_uni = {0};
  size_t m_uni = 0;
  if (upath != NULL) {
    cgfile_t cgf_uni = open_cgfile(upath);
    cg_uni = read_cg(&cgf_uni);
    m_uni = bit_count(cg_uni);
    free(upath);
    bgzf_close(cgf_uni.fh);
  }

  char *fname_fea = argv[optind++];
  cgfile_t cgf_fea = open_cgfile(fname_fea);
  snames_t snames_fea = loadSampleNamesFromIndex(fname_fea);
  int unseekable = bgzf_seek(cgf_fea.fh, 0, SEEK_SET);
  char *fname_qry = argv[optind];
  cgfile_t cgf_qry = open_cgfile(fname_qry);
  snames_t snames_qry = loadSampleNamesFromIndex(fname_qry);

  cgdata_t cg_fea = {0};
  if (unseekable) {             /* only the first cg */
    cg_fea = read_cg(&cgf_fea);
    convertToFmt0(&cg_fea);
  }
  if (print_header)
    fputs("Query\tFeature\tN_universe\tN_feature\tN_query\tN_overlap\n", stdout);
  
  for (uint64_t kq=0;;++kq) {
    cgdata_t cg_qry = read_cg(&cgf_qry);
    if (cg_qry.n == 0) break;
    convertToFmt0(&cg_qry);

    if (cg_uni.n > 0) {
      assert(cg_uni.n == cg_qry.n);
      bit_mask(cg_qry.s, cg_uni.s, cg_uni.n);
    } else { m_uni = cg_qry.n; }
    size_t nf = bit_count(cg_qry);

    if (unseekable) {
      assert(cg_qry.n == cg_fea.n);
      if (cg_uni.n > 0) bit_mask(cg_fea.s, cg_uni.s, cg_uni.n);
      size_t nq = bit_count(cg_fea);
      bit_mask(cg_qry.s, cg_fea.s, cg_qry.n);
      size_t nfq = bit_count(cg_qry);

      if (snames_qry.n) { fputs(snames_qry.s[kq], stdout); fputc('\t', stdout); }
      else fprintf(stdout, "%"PRIu64"\t", kq+1);
      if (snames_fea.n) { fputs(snames_fea.s[0], stdout); fputc('\t', stdout); }
      else fputs("1\t", stdout);
      fprintf(stdout, "%zu\t%zu\t%zu\t%zu\n", m_uni, nf, nq, nfq);
    } else {
      assert(bgzf_seek(cgf_fea.fh, 0, SEEK_SET)==0);
      for (uint64_t kf=0;;++kf) {
        cgdata_t cg_fea = read_cg(&cgf_fea);
        if (cg_fea.n == 0) break;
        convertToFmt0(&cg_fea);
        
        assert(cg_qry.n == cg_fea.n);
        if (cg_uni.n > 0) bit_mask(cg_fea.s, cg_uni.s, cg_uni.n);
        size_t nq = bit_count(cg_fea);
        bit_mask(cg_fea.s, cg_qry.s, cg_qry.n);
        size_t nfq = bit_count(cg_fea);

        if (snames_qry.n) { fputs(snames_qry.s[kq], stdout); fputc('\t', stdout); }
        else fprintf(stdout, "%"PRIu64"\t", kq+1);
        if (snames_fea.n) { fputs(snames_fea.s[kf], stdout); fputc('\t', stdout); }
        else fprintf(stdout, "%"PRIu64"\t", kf+1);
        fprintf(stdout, "%zu\t%zu\t%zu\t%zu\n", m_uni, nf, nq, nfq);

        free(cg_fea.s);
      }
    }
    free(cg_qry.s);
  }
  if (unseekable) free(cg_fea.s);
  if (cg_uni.n > 0) free(cg_uni.s);
  bgzf_close(cgf_qry.fh);
  bgzf_close(cgf_fea.fh);
  cleanSampleNames2(snames_qry);
  cleanSampleNames2(snames_fea);
  
  return 0;
}
