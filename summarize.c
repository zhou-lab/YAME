#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdio.h>
#include "wzmisc.h"
#include "wzbed.h"
#include "cgfile.h"
#include "snames.h"
#include "kstring.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: yame summarize [options] <query.cg>\n");
  fprintf(stderr, "Query should be of format 3, can be a multi-sample set.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -m        Binary (.cb) file, can be multi-sample.\n");
  fprintf(stderr, "    -M        All masks will be loaded to memory. This save disk IO.\n");
  fprintf(stderr, "    -H        Print header?.\n");
  fprintf(stderr, "    -s        Sample list provided to override the query index file.\n");
  fprintf(stderr, "    -h        This help.\n");
  fprintf(stderr, "\n");

  return 1;
}

#define MU2beta(mu) (double) ((mu)>>32) / (((mu)>>32) + ((mu)&0xffffffff))

typedef struct stats_t {
  double mean_beta;
  double sum_beta;
  uint64_t n_u;                 // universe
  uint64_t n_q;                 // query
  uint64_t n_m;                 // mask
  uint64_t n_o;                 // overlap
} stats_t;

static void summarize1(cgdata_t cg, cgdata_t cg_mask, uint64_t *n_st, stats_t **st) {
  if (cg_mask.fmt == '2') {     // state mask
    if (cg_mask.n != cg.n) {
      fprintf(stderr, "[%s:%d] mask (N=%"PRIu64") and query (N=%"PRIu64") are not of the same length.\n", __func__, __LINE__, cg_mask.n, cg.n);
      fflush(stderr);
      exit(1);
    }
    *n_st = fmt2_get_keys_n(&cg_mask);
    (*st) = realloc((*st), sizeof(stats_t)*(*n_st));
    memset(*st, 0, sizeof(stats_t)*(*n_st));
    
    uint64_t nq=0;
    for (uint64_t i=0; i<cg.n; ++i) {
      uint64_t index = f2_unpack_uint64(&cg_mask, i);
      uint64_t mu = f3_unpack_mu(&cg, i);
      if (index >= (*n_st)) {
        fprintf(stderr, "[%s:%d] State data is corrupted.\n", __func__, __LINE__);
        fflush(stderr);
        exit(1);
      }
      if (mu) {
        (*st)[index].sum_beta += MU2beta(mu);
        (*st)[index].n_o++;
        nq++;
      }
      (*st)[index].n_m++;
    }

    for (uint64_t k=0; k < (*n_st); ++k) {
      (*st)[k].n_q = nq;
      (*st)[k].n_u = cg.n;
      (*st)[k].mean_beta = (*st)[k].sum_beta / (*st)[k].n_o;
    }
  } else {                      // binary mask
    if (cg.fmt != '3') {
      fprintf(stderr, "[%s:%d] Query is not sequencing data.\n", __func__, __LINE__);
      fflush(stderr);
      exit(1);
    }
    *n_st = 1;
    (*st) = realloc((*st), sizeof(stats_t)*(*n_st));
    stats_t *st1 = *st;
    memset(st1, 0, sizeof(stats_t));

    st1->n_u = cg.n;
    if (cg_mask.n) {
      if (cg_mask.n != cg.n) {
        fprintf(stderr, "[%s:%d] mask (N=%"PRIu64") and query (N=%"PRIu64") are not of the same length.\n", __func__, __LINE__, cg_mask.n, cg.n);
        fflush(stderr);
        exit(1);
      }
      for (uint64_t i=0; i<cg.n; ++i) {
        uint64_t mu = f3_unpack_mu(&cg, i);
        if (mu) st1->n_q++;
        if (cg_mask.s[i>>3]&(1<<(i&0x7))) {
          st1->n_m++;
          if (mu) {
            st1->sum_beta += MU2beta(mu);
            st1->n_o++;
          }}}
    } else {
      for (uint64_t i=0; i<cg.n; ++i) {
        uint64_t mu = f3_unpack_mu(&cg, i);
        if (mu) {
          st1->sum_beta += MU2beta(mu);
          st1->n_o++;
          st1->n_q++;
        }}
    }
    st1->mean_beta = (double) st1->sum_beta / st1->n_o;
  }
}

static void format_stats(stats_t *st, uint64_t n_st, char *sq, char *sm, cgdata_t *cg_mask) {
  if (cg_mask->fmt == '2') {
    if (!cg_mask->aux) fmt2_set_aux(cg_mask);
    f2_aux_t *aux = (f2_aux_t*) cg_mask->aux;
    for (uint64_t i=0; i<n_st; ++i) {
      fprintf(
        stdout,
        "%s\t%s-%s\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%f\n",
        sq, sm, aux->keys[i],
        st[i].n_u, st[i].n_q, st[i].n_m, st[i].n_o, st[i].mean_beta);
    }
  } else {
    for (uint64_t i=0; i<n_st; ++i) {
      fprintf(
        stdout,
        "%s\t%s\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%f\n",
        sq, sm, st[i].n_u, st[i].n_q, st[i].n_m, st[i].n_o, st[i].mean_beta);
    }
  }
}

/* The design, first 10 bytes are uint64_t (length) + uint16_t (0=vec; 1=rle) */
int main_summarize(int argc, char *argv[]) {
  int c; int print_header = 0; int in_memory = 0;
  char *fname_mask = NULL;
  char *fname_snames = NULL;
  while ((c = getopt(argc, argv, "m:MHs:h"))>=0) {
    switch (c) {
    case 'm': fname_mask = strdup(optarg); break;
    case 'M': in_memory = 1; break;
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
      if (cg_mask.fmt == '2') decompress2(&cg_mask);
      else convertToFmt0(&cg_mask);

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

  cgdata_t *cg_masks = NULL; uint64_t cg_masks_n = 0;
  if (in_memory) {              /* load in-memory masks */
    cg_masks = calloc(1, sizeof(cgdata_t));
    cg_masks_n = 0;
    for (;;++cg_masks_n) {
      cg_mask = read_cg(&cgf_mask);
      if (cg_mask.n == 0) break;
      if (cg_mask.fmt == '2') decompress2(&cg_mask);
      else convertToFmt0(&cg_mask);
      cg_masks = realloc(cg_masks, (cg_masks_n+1)*sizeof(cgdata_t));
      cg_masks[cg_masks_n] = cg_mask;
    }
  }

  uint64_t n_st = 0; stats_t *st = NULL;
  cgdata_t cg_qry_inflated = {0};
  for (uint64_t kq=0;;++kq) {
    cgdata_t cg_qry = read_cg(&cgf_qry);
    if (cg_qry.n == 0) break;
    kstring_t sq = {0};
    if (snames_qry.n) kputs(snames_qry.s[kq], &sq);
    else ksprintf(&sq, "%"PRIu64"\t", kq+1);
    decompress(&cg_qry, &cg_qry_inflated);

    if (fname_mask) {           /* apply any mask? */
      if (in_memory) {          /* in memory mask */
        for (uint64_t km=0;km<cg_masks_n;++km) {
          cg_mask = cg_masks[km];
          summarize1(cg_qry_inflated, cg_mask, &n_st, &st);
          kstring_t sm = {0};
          if (snames_mask.n) kputs(snames_mask.s[km], &sm);
          else ksprintf(&sm, "%"PRIu64"", km+1);
          format_stats(st, n_st, sq.s, sm.s, &cg_mask);
          free(sm.s);
        }
      } else if (unseekable) {  /* mask is unseekable */
        summarize1(cg_qry_inflated, cg_mask, &n_st, &st);
        kstring_t sm = {0};
        kputs("mask", &sm);
        format_stats(st, n_st, sq.s, sm.s, &cg_mask);
        free(sm.s);
      } else {                  /* mask is seekable */
        if (bgzf_seek(cgf_mask.fh, 0, SEEK_SET)!=0) {
          fprintf(stderr, "[%s:%d] Cannot seek mask.\n", __func__, __LINE__);
          fflush(stderr);
          exit(1);
        }
        for (uint64_t km=0;;++km) {
          cg_mask = read_cg(&cgf_mask);
          if (cg_mask.n == 0) break;
          if (cg_mask.fmt == '2') decompress2(&cg_mask);
          else convertToFmt0(&cg_mask);
          summarize1(cg_qry_inflated, cg_mask, &n_st, &st);

          kstring_t sm = {0};
          if (snames_mask.n) kputs(snames_mask.s[km], &sm);
          else ksprintf(&sm, "%"PRIu64"", km+1);
          format_stats(st, n_st, sq.s, sm.s, &cg_mask);
          free(sm.s);
          free_cgdata(&cg_mask);
        }
      }
    } else {                    /* whole dataset summary if missing mask */
      summarize1(cg_qry_inflated, cg_mask, &n_st, &st);
      kstring_t sm = {0};
      kputs("global", &sm);
      format_stats(st, n_st, sq.s, sm.s, &cg_mask);
      free(sm.s);
    }
    free(sq.s);
    free(cg_qry.s); cg_qry.s = NULL;
  }
  if (cg_qry_inflated.s) free(cg_qry_inflated.s);
  if (in_memory) {
    for (uint64_t i=0; i<cg_masks_n; ++i) free_cgdata(&cg_masks[i]);
    free(cg_masks);
  } else if (cg_mask.s) free_cgdata(&cg_mask);
  bgzf_close(cgf_qry.fh);
  if (fname_mask) bgzf_close(cgf_mask.fh);
  cleanSampleNames2(snames_qry);
  cleanSampleNames2(snames_mask);
  
  return 0;
}
