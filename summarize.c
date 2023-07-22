#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdio.h>
#include "wzmisc.h"
#include "wzbed.h"
#include "cgfile.h"
#include "snames.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: kycg summarize [options] <query.cg>\n");
  fprintf(stderr, "Query should be of format 3, can be a multi-sample set.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -m        Binary .cg file, can be multi-sample.\n");
  fprintf(stderr, "    -H        Print header?.\n");
  fprintf(stderr, "    -s        Sample list provided to override the query index file.\n");
  fprintf(stderr, "    -h        This help.\n");
  fprintf(stderr, "\n");

  return 1;
}

#define MU2beta(mu) (double) ((mu)>>32) / (((mu)>>32) + ((mu)&0xffffffff))

typedef struct stats_t {
  double mean_beta;
  uint64_t n_u;                 // universe
  uint64_t n_q;                 // query
  uint64_t n_m;                 // mask
  uint64_t n_o;                 // overlap
} stats_t;

static stats_t summarize1(cgdata_t cg, cgdata_t cg_mask) {
  assert(cg.compressed == 0);
  assert(cg.fmt == '3');
  uint64_t *s = (uint64_t*) cg.s;
  uint64_t sum = 0;

  stats_t stats = {0};
  stats.n_u = cg.n;
  if (cg_mask.n) {
    assert(cg_mask.n == cg.n);
    for (uint64_t i=0; i<cg.n; ++i) {
      if (s[i]) stats.n_q++;
      if (cg_mask.s[i>>3]&(1<<(i&0x7))) {
        stats.n_m++;
        if (s[i]) {
          sum += MU2beta(s[i]);
          stats.n_o++;
        }}}
  } else {
    for (uint64_t i=0; i<cg.n; ++i) {
      if (s[i]) {
        sum += MU2beta(s[i]);
        stats.n_o++;
        stats.n_q++;
      }}
  }
  stats.mean_beta = (double) sum / stats.n_o;
  return stats;
}

static void format_stats(stats_t st) {
  fprintf(stdout,
          "%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%f\n",
          st.n_u, st.n_q, st.n_m, st.n_o, st.mean_beta);
}

/* The design, first 10 bytes are uint64_t (length) + uint16_t (0=vec; 1=rle) */
int main_summarize(int argc, char *argv[]) {
  int c; int print_header = 0;
  char *fname_mask = NULL;
  char *fname_snames = NULL;
  while ((c = getopt(argc, argv, "m:Hs:h"))>=0) {
    switch (c) {
    case 'm': fname_mask = strdup(optarg); break;
    case 'H': print_header = 1; break;
    case 's': fname_snames = strdup(optarg); break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 1 > argc) { 
    usage(); 
    wzfatal("Please supply input file.\n"); 
  }

  cgfile_t cgf_mask; int unseekable = 0; cgdata_t cg_mask = {0};
  snames_t snames_mask = {0};
  if (fname_mask) {
    cgf_mask = open_cgfile(fname_mask);
    unseekable = bgzf_seek(cgf_mask.fh, 0, SEEK_SET);
    if (unseekable) {             /* only the first cg */
      cg_mask = read_cg(&cgf_mask);
      convertToFmt0(&cg_mask);
    }
    snames_mask = loadSampleNamesFromIndex(fname_mask);
  }
  char *fname_qry = argv[optind];
  cgfile_t cgf_qry = open_cgfile(fname_qry);
  snames_t snames_qry = {0};
  
  if (fname_snames) snames_qry = loadSampleNames(fname_snames, 1);
  else snames_qry = loadSampleNamesFromIndex(fname_qry);
  
  if (print_header)
    fputs("Query\tMask\tN_univ\tN_query\tN_mask\tN_overlap\tbeta\n", stdout);
  
  cgdata_t cg_qry_inflated = {0};
  for (uint64_t kq=0;;++kq) {
    cgdata_t cg_qry = read_cg(&cgf_qry);
    if (cg_qry.n == 0) break;
    decompress(&cg_qry, &cg_qry_inflated);

    if (fname_mask) {           /* apply any mask? */
      if (unseekable) {         /* mask is unseekable */
        stats_t st = summarize1(cg_qry_inflated, cg_mask);
        if (snames_qry.n) { fputs(snames_qry.s[kq], stdout); fputc('\t', stdout); }
        else fprintf(stdout, "%"PRIu64"\t", kq+1);
        fputs("mask\t", stdout);
        format_stats(st);
      } else {                  /* mask is seekable */
        assert(bgzf_seek(cgf_mask.fh, 0, SEEK_SET)==0);
        for (uint64_t km=0;;++km) {
          cg_mask = read_cg(&cgf_mask);
          if (cg_mask.n == 0) break;
          convertToFmt0(&cg_mask);
          stats_t st = summarize1(cg_qry_inflated, cg_mask);

          if (snames_qry.n) { fputs(snames_qry.s[kq], stdout); fputc('\t', stdout); }
          else fprintf(stdout, "%"PRIu64"\t", kq+1);
          if (snames_mask.n) { fputs(snames_mask.s[km], stdout); fputc('\t', stdout); }
          else fprintf(stdout, "%"PRIu64"\t", km+1);
          format_stats(st);
          free(cg_mask.s); cg_mask.s = 0;
        }
      }
    } else {                    /* whole dataset summary if missing mask */
      stats_t st = summarize1(cg_qry_inflated, cg_mask);
      if (snames_qry.n) { fputs(snames_qry.s[kq], stdout); fputc('\t', stdout); }
      else fprintf(stdout, "%"PRIu64"\t", kq+1);
      fputs("global\t", stdout);
      format_stats(st);
    }
    free(cg_qry.s); cg_qry.s = NULL;
  }
  if (cg_qry_inflated.s) free(cg_qry_inflated.s);
  if (cg_mask.s) free(cg_mask.s);
  bgzf_close(cgf_qry.fh);
  if (fname_mask) bgzf_close(cgf_mask.fh);
  cleanSampleNames2(snames_qry);
  cleanSampleNames2(snames_mask);

  /* char *fname_fea = argv[optind++]; */
  /* cgfile_t cgf_fea = open_cgfile(fname_fea); int n_fea=0; */
  /* index_pair_t *idx_pairs_fea = load_index_pairs(fname_fea, &n_fea); */
  /* int unseekable = bgzf_seek(cgf_fea.fh, 0, SEEK_SET); */
  /* char *fname_qry = argv[optind]; */
  /* cgfile_t cgf_qry = open_cgfile(fname_qry); int n_qry=0; */
  /* index_pair_t *idx_pairs_qry = load_index_pairs(fname_qry, &n_qry); */

  
  
  
  /* for (uint64_t kq=0;;++kq) { */
  /*   cgdata_t cg_qry = read_cg(&cgf_qry); */
  /*   if (cg_qry.n == 0) break; */
  /*   convertToFmt0(&cg_qry); */

  /*   if (cg_uni.n > 0) { */
  /*     assert(cg_uni.n == cg_qry.n); */
  /*     bit_mask(cg_qry.s, cg_uni.s, cg_uni.n); */
  /*   } else { m_uni = cg_qry.n; } */
  /*   size_t nf = bit_count(cg_qry); */

  /*   if (unseekable) { */
  /*     assert(cg_qry.n == cg_fea.n); */
  /*     if (cg_uni.n > 0) bit_mask(cg_fea.s, cg_uni.s, cg_uni.n); */
  /*     size_t nq = bit_count(cg_fea); */
  /*     bit_mask(cg_qry.s, cg_fea.s, cg_qry.n); */
  /*     size_t nfq = bit_count(cg_qry); */

  /*     if (idx_pairs_qry) { fputs(idx_pairs_qry[kq].key, stdout); fputc('\t', stdout); } */
  /*     else fprintf(stdout, "%"PRIu64"\t", kq+1); */
  /*     if (idx_pairs_fea) { fputs(idx_pairs_fea[0].key, stdout); fputc('\t', stdout); } */
  /*     else fputs("1\t", stdout); */
  /*     fprintf(stdout, "%zu\t%zu\t%zu\t%zu\n", m_uni, nf, nq, nfq); */
  /*   } else { */
  /*     assert(bgzf_seek(cgf_fea.fh, 0, SEEK_SET)==0); */
  /*     for (uint64_t kf=0;;++kf) { */
  /*       cgdata_t cg_fea = read_cg(&cgf_fea); */
  /*       if (cg_fea.n == 0) break; */
  /*       convertToFmt0(&cg_fea); */
        
  /*       assert(cg_qry.n == cg_fea.n); */
  /*       if (cg_uni.n > 0) bit_mask(cg_fea.s, cg_uni.s, cg_uni.n); */
  /*       size_t nq = bit_count(cg_fea); */
  /*       bit_mask(cg_fea.s, cg_qry.s, cg_qry.n); */
  /*       size_t nfq = bit_count(cg_fea); */

  /*       if (idx_pairs_qry) { fputs(idx_pairs_qry[kq].key, stdout); fputc('\t', stdout); } */
  /*       else fprintf(stdout, "%"PRIu64"\t", kq+1); */
  /*       if (idx_pairs_fea) { fputs(idx_pairs_fea[kf].key, stdout); fputc('\t', stdout); } */
  /*       else fprintf(stdout, "%"PRIu64"\t", kf+1); */
  /*       fprintf(stdout, "%zu\t%zu\t%zu\t%zu\n", m_uni, nf, nq, nfq); */

  /*       free(cg_fea.s); */
  /*     } */
  /*   } */
  /*   free(cg_qry.s); */
  /* } */
  /* if (unseekable) free(cg_fea.s); */
  /* if (cg_uni.n > 0) free(cg_uni.s); */

  /* bgzf_close(cgf_fea.fh); */
  /* if (idx_pairs_fea) clean_index_pairs(idx_pairs_fea, n_fea); */
  /* if (idx_pairs_qry) clean_index_pairs(idx_pairs_qry, n_qry); */
  
  return 0;
}
