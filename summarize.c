#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdio.h>
#include "wzmisc.h"
#include "wzbed.h"
#include "cgfile.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: kycg summarize [options] <query.cg>\n");
  fprintf(stderr, "Query can be a multi-sample set.\n");
  fprintf(stderr, "If the input is format 3: test summarize of M+U.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -m        binary .cg file, can be multi-sample.\n");
  fprintf(stderr, "    -H        print header.\n");
  fprintf(stderr, "    -h        This help.\n");
  fprintf(stderr, "Returns four numbers, nu, nf, nq, nfq;\n");
  fprintf(stderr, "\n");

  return 1;
}

#define MU2beta(mu) (double) ((mu)>>32) / (((mu)>>32) + ((mu)&0xffffffff))

static double meanBetas(cgdata_t cg, uint64_t *n_nonzero, cgdata_t cg_mask) {
  assert(cg.compressed == 0);
  assert(cg.fmt == '3');
  uint64_t *s = (uint64_t*) cg.s;
  uint64_t sum = 0; uint64_t n = 0;

  if (cg_mask.n) {
    for (uint64_t i=0; i<cg.n; ++i) {
      if (s[i] && cg_mask.s[i>>3]&(1<<(i&0x7))) { sum += MU2beta(s[i]); n++; }}
  } else {
    for (uint64_t i=0; i<cg.n; ++i) { if (s[i]) { sum += MU2beta(s[i]); n++; }}
  }
  *n_nonzero = n;
  return ((double)sum)/n;
}

/* The design, first 10 bytes are uint64_t (length) + uint16_t (0=vec; 1=rle) */
int main_summarize(int argc, char *argv[]) {
  int c; char *fname_mask = NULL; int print_header = 0;
  while ((c = getopt(argc, argv, "Hm:h"))>=0) {
    switch (c) {
    case 'm': fname_mask = strdup(optarg); break;
    case 'H': print_header = 1; break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 1 > argc) { 
    usage(); 
    wzfatal("Please supply input file.\n"); 
  }

  cgfile_t cgf_mask; int unseekable = 0; cgdata_t cg_mask = {0};
  index_pair_t *idx_pairs_mask = NULL; int n_mask = 0;
  if (fname_mask) {
    cgf_mask = open_cgfile(fname_mask);
    unseekable = bgzf_seek(cgf_mask.fh, 0, SEEK_SET);
    if (unseekable) {             /* only the first cg */
      cg_mask = read_cg(&cgf_mask);
      convertToFmt0(&cg_mask);
    }
    idx_pairs_mask = load_index_pairs(fname_mask, &n_mask);
  }
  char *fname_qry = argv[optind];
  cgfile_t cgf_qry = open_cgfile(fname_qry); int n_qry = 0;
  index_pair_t *idx_pairs_qry = load_index_pairs(fname_qry, &n_qry);
  if (print_header) fputs("Query\tFeature\tN_nonzero\tbeta\n", stdout);

  uint64_t n_nonzero = 0; double mean_betas = 0.0;
  cgdata_t cg_qry_inflated = {0};
  for (uint64_t kq=0;;++kq) {
    cgdata_t cg_qry = read_cg(&cgf_qry);
    if (cg_qry.n == 0) break;
    decompress(&cg_qry, &cg_qry_inflated);

    if (fname_mask) {           /* apply any mask? */
      if (unseekable) {         /* mask is unseekable */
        assert(cg_mask.n == cg_qry_inflated.n);
        mean_betas = meanBetas(cg_qry_inflated, &n_nonzero, cg_mask);
        if (idx_pairs_qry) { fputs(idx_pairs_qry[kq].key, stdout); fputc('\t', stdout); }
        else fprintf(stdout, "%"PRIu64"\t", kq+1);
        fprintf(stdout, "Masked\t%"PRIu64"\t%f\n", n_nonzero, mean_betas);
      } else {                  /* mask is seekable */
        assert(bgzf_seek(cgf_mask.fh, 0, SEEK_SET)==0);
        for (uint64_t km=0;;++km) {
          cg_mask = read_cg(&cgf_mask);
          if (cg_mask.n == 0) break;
          convertToFmt0(&cg_mask);
          assert(cg_mask.n == cg_qry_inflated.n);
          mean_betas = meanBetas(cg_qry_inflated, &n_nonzero, cg_mask);

          if (idx_pairs_qry) { fputs(idx_pairs_qry[kq].key, stdout); fputc('\t', stdout); }
          else fprintf(stdout, "%"PRIu64"\t", kq+1);
          if (idx_pairs_mask) { fputs(idx_pairs_mask[km].key, stdout); fputc('\t', stdout); }
          else fprintf(stdout, "%"PRIu64"\t", km+1);
          fprintf(stdout, "%"PRIu64"\t%f\n", n_nonzero, mean_betas);
          free(cg_mask.s); cg_mask.s = 0;
        }
      }
    } else {
      mean_betas = meanBetas(cg_qry_inflated, &n_nonzero, cg_mask);
      if (idx_pairs_qry) { fputs(idx_pairs_qry[kq].key, stdout); fputc('\t', stdout); }
      else fprintf(stdout, "%"PRIu64"\t", kq+1);
      fprintf(stdout, "global\t%"PRIu64"\t%f\n", n_nonzero, mean_betas);
    }
    free(cg_qry.s); cg_qry.s = NULL;
  }
  if (cg_qry_inflated.s) free(cg_qry_inflated.s);
  if (cg_mask.s) free(cg_mask.s);
  bgzf_close(cgf_qry.fh);
  if (fname_mask) bgzf_close(cgf_mask.fh);
  if (idx_pairs_qry) clean_index_pairs(idx_pairs_qry, n_qry);
  if (idx_pairs_mask) clean_index_pairs(idx_pairs_mask, n_mask);

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
