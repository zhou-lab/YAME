#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdio.h>
#include "wzmisc.h"
#include "wzbed.h"
#include "cfile.h"
#include "snames.h"
#include "kstring.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: yame summary [options] <query.cm>\n");
  fprintf(stderr, "Query should be of format 0,1,2,3, can be a multi-sample set.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -m        Binary (.cb) file, can be multi-sample.\n");
  fprintf(stderr, "              If '-', the whole sample will bed kept in memory, same as -M.\n");
  fprintf(stderr, "    -M        All masks will be loaded to memory. This save disk IO.\n");
  fprintf(stderr, "    -H        Suppress header printing.\n");
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
  char* sm;                     // mask name
  char* sq;                     // query name
} stats_t;

static stats_t* summarize1_fmt1(cdata_t c, cdata_t c_mask, uint64_t *n_st, char *sm, char *sq) {

  stats_t *st = NULL;
  if (c_mask.n == 0) {          // no mask

    *n_st = 1;
    st = calloc(1, sizeof(stats_t));
    st[0].n_u = c.n;
    st[0].n_q = bit_count(c);
    st[0].sm = strdup(sm);
    st[0].sq = strdup(sq);
    
  } else if (c_mask.fmt <= '1') { // binary mask

    if (c_mask.n != c.n) {
      fprintf(stderr, "[%s:%d] mask (N=%"PRIu64") and query (N=%"PRIu64") are of different lengths.\n", __func__, __LINE__, c_mask.n, c.n);
      fflush(stderr);
      exit(1);
    }
    
    *n_st = 1;
    st = calloc(1, sizeof(stats_t));
    st[0].n_u = c.n;
    st[0].n_q = bit_count(c);
    st[0].n_m = bit_count(c_mask);
    cdata_t tmp = {0};
    tmp.s = malloc((c.n>>3)+1); tmp.n = c.n;
    memcpy(tmp.s, c.s, (c.n>>3)+1);
    bit_mask(tmp.s, c_mask.s, c_mask.n);
    st[0].n_o = bit_count(tmp);
    free(tmp.s);
    st[0].sm = strdup(sm);
    st[0].sq = strdup(sq);
    st[0].mean_beta = -1;

  } else if (c_mask.fmt == '2') { // state mask

    if (c_mask.n != c.n) {
      fprintf(stderr, "[%s:%d] mask (N=%"PRIu64") and query (N=%"PRIu64") are of different lengths.\n", __func__, __LINE__, c_mask.n, c.n);
      fflush(stderr);
      exit(1);
    }
    if (!c_mask.aux) fmt2_set_aux(&c_mask);
    f2_aux_t *aux = (f2_aux_t*) c_mask.aux;
    *n_st = aux->nk;
    st = calloc((*n_st), sizeof(stats_t));
    uint64_t nq=0;
    for (uint64_t i=0; i<c.n; ++i) {
      uint64_t index = f2_get_uint64(&c_mask, i);
      if (index >= (*n_st)) {
        fprintf(stderr, "[%s:%d] State data is corrupted.\n", __func__, __LINE__);
        fflush(stderr);
        exit(1);
      }
      if (c.s[i>>3]&(1<<(i&0x7))) {
        st[index].n_o++;
        nq++;
      }
      st[index].n_m++;
    }
    for (uint64_t k=0; k < (*n_st); ++k) {
      st[k].n_q = nq;
      st[k].n_u = c.n;
      st[k].sum_beta = -1;
      st[k].mean_beta = -1;
      kstring_t tmp = {0};
      ksprintf(&tmp, "%s-%s", sm, aux->keys[k]);
      st[k].sm = tmp.s;
      st[k].sq = strdup(sq);
    }
    
  } else {                      // other masks
    fprintf(stderr, "[%s:%d] Mask format %c unsupported.\n", __func__, __LINE__, c_mask.fmt);
    fflush(stderr);
    exit(1);
  }
  return st;
}

static stats_t* summarize1_fmt2(cdata_t c, cdata_t c_mask, uint64_t *n_st, char *sm, char *sq) {

  stats_t *st = NULL;
  if (c_mask.n == 0) {          // no mask
    
    if (!c.aux) fmt2_set_aux(&c);
    f2_aux_t *aux = (f2_aux_t*) c.aux;
    *n_st = aux->nk;
    uint64_t *cnts = calloc(aux->nk, sizeof(uint64_t));
    for (uint64_t i=0; i<c.n; ++i) cnts[f2_get_uint64(&c, i)]++;
    st = calloc(aux->nk, sizeof(stats_t));
    for (uint64_t k=0; k<aux->nk; ++k) {
      st[k].n_u = c.n;
      st[k].n_q = cnts[k];
      st[k].n_m = 0;
      st[k].n_o = 0;
      st[k].mean_beta = -1;
      st[k].sum_beta = -1;
      st[k].sm = strdup(sm);
      kstring_t tmp = {0};
      ksprintf(&tmp, "%s-%s", sq, aux->keys[k]);
      st[k].sq = tmp.s;
    }
    free(cnts);
    
  } else if (c_mask.fmt <= '1') { // binary mask

    if (!c.aux) fmt2_set_aux(&c);
    f2_aux_t *aux = (f2_aux_t*) c.aux;
    *n_st = aux->nk;
    uint64_t *cnts = calloc(aux->nk, sizeof(uint64_t));
    uint64_t *cnts_q = calloc(aux->nk, sizeof(uint64_t));
    uint64_t n_m = 0;
    for (uint64_t i=0; i<c.n; ++i) {
      if (c_mask.s[i>>3]&(1<<(i&0x7))) {
        n_m++;
        cnts[f2_get_uint64(&c, i)]++;
      }
      cnts_q[f2_get_uint64(&c, i)]++;
    }
    st = calloc(aux->nk, sizeof(stats_t));
    for (uint64_t k=0; k<aux->nk; ++k) {
      st[k].n_u = c.n;
      st[k].n_q = cnts_q[k];
      st[k].n_o = cnts[k];
      st[k].n_m = n_m;
      st[k].mean_beta = -1;
      st[k].sum_beta = -1;
      st[k].sm = strdup(sm);
      kstring_t tmp = {0};
      ksprintf(&tmp, "%s-%s", sq, aux->keys[k]);
      st[k].sq = tmp.s;
    }
    free(cnts);
    
  } else if (c_mask.fmt == '2') { // state mask

    if (c_mask.n != c.n) {
      fprintf(stderr, "[%s:%d] mask (N=%"PRIu64") and query (N=%"PRIu64") are of different lengths.\n", __func__, __LINE__, c_mask.n, c.n);
      fflush(stderr);
      exit(1);
    }

    if (!c_mask.aux) fmt2_set_aux(&c_mask);
    f2_aux_t *aux_m = (f2_aux_t*) c_mask.aux;

    if (!c.aux) fmt2_set_aux(&c);
    f2_aux_t *aux_q = (f2_aux_t*) c.aux;

    *n_st = aux_m->nk * aux_q->nk;
    st = calloc((*n_st), sizeof(stats_t));
    uint64_t *nq = calloc(aux_q->nk, sizeof(uint64_t));
    uint64_t *nm = calloc(aux_m->nk, sizeof(uint64_t));
    for (uint64_t i=0; i<c.n; ++i) {
      uint64_t im = f2_get_uint64(&c_mask, i);
      uint64_t iq = f2_get_uint64(&c, i);
      st[im * aux_q->nk + iq].n_o++;
      nq[iq]++;
      nm[im]++;
    }
    for (uint64_t im=0; im<aux_m->nk; ++im) {
      for (uint64_t iq=0; iq<aux_q->nk; ++iq) {
        stats_t *st1 = &st[im * aux_q->nk + iq];
        st1->n_o++;
        st1->n_u = c.n;
        st1->n_q = nq[iq];
        st1->n_m = nm[im];
        st1->mean_beta = -1;
        kstring_t tmp = {0};
        ksprintf(&tmp, "%s-%s", sm, aux_m->keys[im]);
        st1->sm = tmp.s;
        memset(&tmp, 0, sizeof(kstring_t));
        ksprintf(&tmp, "%s-%s", sq, aux_q->keys[iq]);
        st1->sq = tmp.s;
      }
    }
    free(nq); free(nm);
    
  } else {                      // other masks
    fprintf(stderr, "[%s:%d] Mask format %c unsupported.\n", __func__, __LINE__, c_mask.fmt);
    fflush(stderr);
    exit(1);
  }
  return st;
}

static stats_t* summarize1_fmt3(cdata_t c, cdata_t c_mask, uint64_t *n_st, char *sm, char *sq) {

  stats_t *st = NULL;
  if (c_mask.n == 0) {            // no mask
    
    *n_st = 1;
    st = calloc(1, sizeof(stats_t));
    st[0].n_u = c.n;
    for (uint64_t i=0; i<c.n; ++i) {
      uint64_t mu = f3_get_mu(&c, i);
      if (mu) {
        st[0].sum_beta += MU2beta(mu);
        st[0].n_o++;
        st[0].n_q++;
      }}
    st[0].sm = strdup(sm);
    st[0].sq = strdup(sq);
    
  } else if (c_mask.fmt <= '1') { // binary mask
    
    *n_st = 1;
    st = calloc(1, sizeof(stats_t));
    st[0].n_u = c.n;
    if (c_mask.n != c.n) {
      fprintf(stderr, "[%s:%d] mask (N=%"PRIu64") and query (N=%"PRIu64") are of different lengths.\n", __func__, __LINE__, c_mask.n, c.n);
      fflush(stderr);
      exit(1);
    }
    for (uint64_t i=0; i<c.n; ++i) {
      uint64_t mu = f3_get_mu(&c, i);
      if (mu) st[0].n_q++;
      if (c_mask.s[i>>3]&(1<<(i&0x7))) {
        st[0].n_m++;
        if (mu) {
          st[0].sum_beta += MU2beta(mu);
          st[0].n_o++;
        }}}
    st[0].mean_beta = (double) st[0].sum_beta / st[0].n_o;
    st[0].sm = strdup(sm);
    st[0].sq = strdup(sq);
    
  } else if (c_mask.fmt == '2') { // state mask
    
    if (c_mask.n != c.n) {
      fprintf(stderr, "[%s:%d] mask (N=%"PRIu64") and query (N=%"PRIu64") are of different lengths.\n", __func__, __LINE__, c_mask.n, c.n);
      fflush(stderr);
      exit(1);
    }
    if (!c_mask.aux) fmt2_set_aux(&c_mask);
    f2_aux_t *aux = (f2_aux_t*) c_mask.aux;
    *n_st = aux->nk;
    st = calloc((*n_st), sizeof(stats_t));
    uint64_t nq=0;
    for (uint64_t i=0; i<c.n; ++i) {
      uint64_t index = f2_get_uint64(&c_mask, i);
      uint64_t mu = f3_get_mu(&c, i);
      if (index >= (*n_st)) {
        fprintf(stderr, "[%s:%d] State data is corrupted.\n", __func__, __LINE__);
        fflush(stderr);
        exit(1);
      }
      if (mu) {
        st[index].sum_beta += MU2beta(mu);
        st[index].n_o++;
        nq++;
      }
      st[index].n_m++;
    }
    for (uint64_t k=0; k < (*n_st); ++k) {
      st[k].n_q = nq;
      st[k].n_u = c.n;
      st[k].mean_beta = st[k].sum_beta / st[k].n_o;
      kstring_t tmp = {0};
      ksprintf(&tmp, "%s-%s", sm, aux->keys[k]);
      st[k].sm = tmp.s;
      st[k].sq = strdup(sq);
    }
    
  } else {                      // other masks
    fprintf(stderr, "[%s:%d] Mask format %c unsupported.\n", __func__, __LINE__, c_mask.fmt);
    fflush(stderr);
    exit(1);
  }
  return st;
}

static stats_t* summarize1(cdata_t c, cdata_t c_mask, uint64_t *n_st, char *sm, char *sq) {
  if (c.fmt == '3') {
    return summarize1_fmt3(c, c_mask, n_st, sm, sq);
  } else if (c.fmt == '2') {
    return summarize1_fmt2(c, c_mask, n_st, sm, sq);
  } else if (c.fmt == '0') {
    return summarize1_fmt1(c, c_mask, n_st, sm, sq);
  } else {
    fprintf(stderr, "[%s:%d] Query format %c unsupported.\n", __func__, __LINE__, c.fmt);
    fflush(stderr);
    exit(1);
  }
}

static void format_stats_and_clean(stats_t *st, uint64_t n_st, const char *fname_mask, const char *fname_qry) {
  const char *fmask = "NA";
  for (uint64_t i=0; i<n_st; ++i) {
    stats_t s = st[i];
    char *odds_ratio = NULL;
    if (fname_mask) {
      double n_mm = s.n_u - s.n_q - s.n_m + s.n_o;
      double n_mp = s.n_q - s.n_o;
      double n_pm = s.n_m - s.n_o;
      kstring_t tmp = {0};
      ksprintf(&tmp, "%1.2f", log2(n_mm*s.n_o / (n_mp*n_pm)));
      odds_ratio = tmp.s;
      fmask = fname_mask;
    } else {
      kstring_t tmp = {0};
      kputs("NA", &tmp);
      odds_ratio = tmp.s;
      fmask = "NA";
    }
    fprintf(stdout,
      "%s\t%s\t%s\t%s\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%s\t",
      fname_qry, s.sq, fmask, s.sm, s.n_u, s.n_q, s.n_m, s.n_o, odds_ratio);
    if (s.mean_beta < 0) {
      fputs("NA", stdout);
    } else {
      fprintf(stdout, "%1.3f", s.mean_beta);
    }
    fputc('\n', stdout);
    free(odds_ratio);
  }

  if (n_st) {
    for (uint64_t i=0; i<n_st; ++i) {
      free(st[i].sm);
      free(st[i].sq);
    }
    free(st);
  }
}

static void prepare_query(cdata_t *c) {
  if (c->fmt < '2') {
    convertToFmt0(c);
  } else {
    decompress2(c);
  }
}

static void prepare_mask(cdata_t *c) {
  if (c->fmt < '2') {
    convertToFmt0(c);
  } else {
    decompress2(c);
  }
}

/* The design, first 10 bytes are uint64_t (length) + uint16_t (0=vec; 1=rle) */
int main_summary(int argc, char *argv[]) {
  int c; int print_header = 1; int in_memory = 0;
  char *fname_mask = NULL;
  char *fname_snames = NULL;
  while ((c = getopt(argc, argv, "m:MHs:h"))>=0) {
    switch (c) {
    case 'm': fname_mask = strdup(optarg); break;
    case 'M': in_memory = 1; break;
    case 'H': print_header = 0; break;
    case 's': fname_snames = strdup(optarg); break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 1 > argc) { 
    usage(); 
    wzfatal("Please supply input file.\n"); 
  }

  cfile_t cf_mask; int unseekable = 0; cdata_t c_mask = {0};
  snames_t snames_mask = {0};
  if (fname_mask) {
    cf_mask = open_cfile(fname_mask);
    unseekable = bgzf_seek(cf_mask.fh, 0, SEEK_SET);
    if (unseekable) {             /* only the first c */
      c_mask = read_cdata1(&cf_mask);
      prepare_mask(&c_mask);
    }
    snames_mask = loadSampleNamesFromIndex(fname_mask);
  }
 
  if (print_header)
    fputs("QFile\tQuery\tMFile\tMask\tN_univ\tN_query\tN_mask\tN_overlap\tLog2OddsRatio\tBeta\n", stdout);

  cdata_t *c_masks = NULL; uint64_t c_masks_n = 0;
  if (in_memory) {              /* load in-memory masks */
    c_masks = calloc(1, sizeof(cdata_t));
    c_masks_n = 0;
    for (;;++c_masks_n) {
      c_mask = read_cdata1(&cf_mask);
      if (c_mask.n == 0) break;
      prepare_mask(&c_mask);
      c_masks = realloc(c_masks, (c_masks_n+1)*sizeof(cdata_t));
      c_masks[c_masks_n] = c_mask;
    }
  }

  for (int j = optind; j < argc; ++j) {
    
    char *fname_qry = argv[j];
    cfile_t cf_qry = open_cfile(fname_qry);
    snames_t snames_qry = {0};
  
    if (fname_snames) snames_qry = loadSampleNames(fname_snames, 1);
    else snames_qry = loadSampleNamesFromIndex(fname_qry);

    for (uint64_t kq=0;;++kq) {
      cdata_t c_qry = read_cdata1(&cf_qry);
      if (c_qry.n == 0) break;
      kstring_t sq = {0};
      if (snames_qry.n) kputs(snames_qry.s[kq], &sq);
      else ksprintf(&sq, "%"PRIu64"", kq+1);
      prepare_query(&c_qry);

      if (fname_mask) {           /* apply any mask? */
        if (in_memory || unseekable) { /* in memory mask if unseekable */
          for (uint64_t km=0;km<c_masks_n;++km) {
            c_mask = c_masks[km];
            kstring_t sm = {0};
            if (snames_mask.n) kputs(snames_mask.s[km], &sm);
            else ksprintf(&sm, "%"PRIu64"", km+1);
            uint64_t n_st = 0;
            stats_t *st = summarize1(c_qry, c_mask, &n_st, sm.s, sq.s);
            format_stats_and_clean(st, n_st, fname_mask, fname_qry);
            free(sm.s);
          }
        } else {                  /* mask is seekable */
          if (bgzf_seek(cf_mask.fh, 0, SEEK_SET)!=0) {
            fprintf(stderr, "[%s:%d] Cannot seek mask.\n", __func__, __LINE__);
            fflush(stderr);
            exit(1);
          }
          for (uint64_t km=0;;++km) {
            c_mask = read_cdata1(&cf_mask);
            if (c_mask.n == 0) break;
            prepare_mask(&c_mask);

            kstring_t sm = {0};
            if (snames_mask.n) kputs(snames_mask.s[km], &sm);
            else ksprintf(&sm, "%"PRIu64"", km+1);
            uint64_t n_st = 0;
            stats_t *st = summarize1(c_qry, c_mask, &n_st, sm.s, sq.s);
            format_stats_and_clean(st, n_st, fname_mask, fname_qry);
            free(sm.s);
            free_cdata(&c_mask);
          }
        }
      } else {                    /* whole dataset summary if missing mask */
        kstring_t sm = {0};
        kputs("global", &sm);
        uint64_t n_st = 0;
        stats_t *st = summarize1(c_qry, c_mask, &n_st, sm.s, sq.s);
        format_stats_and_clean(st, n_st, fname_mask, fname_qry);
        free(sm.s);
      }
      free(sq.s);
      free(c_qry.s); c_qry.s = NULL;
    }
    if (in_memory) {
      for (uint64_t i=0; i<c_masks_n; ++i) free_cdata(&c_masks[i]);
      free(c_masks);
    } else if (c_mask.s) free_cdata(&c_mask);
    bgzf_close(cf_qry.fh);
    cleanSampleNames2(snames_qry);
    if (fname_snames) free(fname_snames);
  }
  if (fname_mask) bgzf_close(cf_mask.fh);
  if (fname_mask) free(fname_mask);
  cleanSampleNames2(snames_mask);
  
  return 0;
}
