#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdio.h>
#include "wzmisc.h"
#include "wzbed.h"
#include "cfile.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: yame overlap [options] <feature.cx> <query.cm>\n");
  fprintf(stderr, "Query or feature can be a multi-sample set.\n");
  fprintf(stderr, "If feature.cx is unseekable, only first sample will be used.\n");
  fprintf(stderr, "If the input is format 3: test overlap of M+U.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -u        Universe .cb file (optional).\n");
  fprintf(stderr, "    -h        This help.\n");
  fprintf(stderr, "Returns four numbers, nu, nf, nq, nfq;\n");
  fprintf(stderr, "\n");

  return 1;
}

static void bit_mask(uint8_t *s, uint8_t *mask, size_t n) {
  size_t i;
  for (i=0; i<(n>>3)+1; ++i) s[i] &= mask[i];
}

static size_t bit_count(cdata_t c) {

  /* create a look-up table */
  int byte2cnt[256]; int p;
  for (p=0; p<256; ++p) {
    unsigned char q = p; int ii, cnt = 0;
    for (ii=0; ii<8; ++ii) { if (q&1) cnt++; q>>=1; }
    byte2cnt[p] = cnt;
  }
  
  size_t i,k,m = 0;
  for (i=0; i<(c.n>>3); ++i) m += byte2cnt[c.s[i]];
  for (k=0; k<(c.n&0x7); ++k) m += (c.s[i]>>k) & 0x1;
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
  cdata_t c_uni = {0};
  size_t m_uni = 0;
  if (upath != NULL) {
    cfile_t cf_uni = open_cfile(upath);
    c_uni = read_cdata1(&cf_uni);
    m_uni = bit_count(c_uni);
    free(upath);
    bgzf_close(cf_uni.fh);
  }

  char *fname_fea = argv[optind++];
  cfile_t cf_fea = open_cfile(fname_fea);
  snames_t snames_fea = loadSampleNamesFromIndex(fname_fea);
  int unseekable = bgzf_seek(cf_fea.fh, 0, SEEK_SET);
  char *fname_qry = argv[optind];
  cfile_t cf_qry = open_cfile(fname_qry);
  snames_t snames_qry = loadSampleNamesFromIndex(fname_qry);

  cdata_t c_fea = {0};
  if (unseekable) {             /* only the first cdata */
    c_fea = read_cdata1(&cf_fea);
    if (c_fea.fmt == '2') decompress2(&c_fea);
    else convertToFmt0(&c_fea);
  }
  if (print_header)
    fputs("Query\tFeature\tN_universe\tN_feature\tN_query\tN_overlap\n", stdout);
  
  for (uint64_t kq=0;;++kq) {
    cdata_t c_qry = read_cdata1(&cf_qry);
    if (c_qry.n == 0) break;
    convertToFmt0(&c_qry);

    if (c_uni.n > 0) {
      assert(c_uni.n == c_qry.n);
      bit_mask(c_qry.s, c_uni.s, c_uni.n);
    } else { m_uni = c_qry.n; }
    size_t nf = bit_count(c_qry);

    if (unseekable) {
      assert(c_qry.n == c_fea.n);
      if (c_uni.n > 0) bit_mask(c_fea.s, c_uni.s, c_uni.n);
      size_t nq = bit_count(c_fea);
      bit_mask(c_qry.s, c_fea.s, c_qry.n);
      size_t nfq = bit_count(c_qry);

      if (snames_qry.n) { fputs(snames_qry.s[kq], stdout); fputc('\t', stdout); }
      else fprintf(stdout, "%"PRIu64"\t", kq+1);
      if (snames_fea.n) { fputs(snames_fea.s[0], stdout); fputc('\t', stdout); }
      else fputs("1\t", stdout);
      fprintf(stdout, "%zu\t%zu\t%zu\t%zu\n", m_uni, nf, nq, nfq);
    } else {
      assert(bgzf_seek(cf_fea.fh, 0, SEEK_SET)==0);
      for (uint64_t kf=0;;++kf) {
        cdata_t c_fea = read_cdata1(&cf_fea);
        if (c_fea.n == 0) break;
        if (c_fea.fmt == '2') decompress2(&c_fea);
        else convertToFmt0(&c_fea);
        
        assert(c_qry.n == c_fea.n);
        if (c_uni.n > 0) bit_mask(c_fea.s, c_uni.s, c_uni.n);
        size_t nq = bit_count(c_fea);
        bit_mask(c_fea.s, c_qry.s, c_qry.n);
        size_t nfq = bit_count(c_fea);

        if (snames_qry.n) { fputs(snames_qry.s[kq], stdout); fputc('\t', stdout); }
        else fprintf(stdout, "%"PRIu64"\t", kq+1);
        if (snames_fea.n) { fputs(snames_fea.s[kf], stdout); fputc('\t', stdout); }
        else fprintf(stdout, "%"PRIu64"\t", kf+1);
        fprintf(stdout, "%zu\t%zu\t%zu\t%zu\n", m_uni, nf, nq, nfq);

        free(c_fea.s);
      }
    }
    free(c_qry.s);
  }
  if (unseekable) free(c_fea.s);
  if (c_uni.n > 0) free(c_uni.s);
  bgzf_close(cf_qry.fh);
  bgzf_close(cf_fea.fh);
  cleanSampleNames2(snames_qry);
  cleanSampleNames2(snames_fea);
  
  return 0;
}
