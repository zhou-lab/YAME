// SPDX-License-Identifier: AGPL-3.0-or-later
/**
 * This file is part of YAME.
 *
 * Copyright (C) 2021-present Wanding Zhou
 *
 * YAME is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * YAME is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with YAME.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdio.h>
#include "wzmisc.h"
#include "wzbed.h"
#include "cfile.h"
#include "snames.h"
#include "summary.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: yame summary [options] <query.cm>\n");
  fprintf(stderr, "Query should be of format 0,1,2,3, can be a multi-sample set.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -m        Mask feature (.cx) file, can be multi-sample.\n");
  fprintf(stderr, "              If '-', the whole sample will bed kept in memory, same as -M.\n");
  fprintf(stderr, "    -M        All masks will be loaded to memory. This save disk IO.\n");
  fprintf(stderr, "    -u        Optional universe set as a .cx file. If given, the masks and queries are both subset.\n");
  fprintf(stderr, "    -H        Suppress header printing.\n");
  fprintf(stderr, "    -q        The backup query file name if the query file name is '-'.\n");
  fprintf(stderr, "    -F        Use full feature/query file name instead of base name.\n");
  fprintf(stderr, "    -T        State features always show section names.\n");
  fprintf(stderr, "    -s        Sample list provided to override the query index file. Only applies to the first query.\n");
  fprintf(stderr, "    -h        This help.\n");
  fprintf(stderr, "\n");

  return 1;
}

static stats_t* summarize1_queryfmt0(
  cdata_t *c, cdata_t *c_mask, uint64_t *n_st, char *sm, char *sq, config_t *config) {

  stats_t *st = NULL;
  if (c_mask->n == 0) {          // no mask
    
    *n_st = 1;
    st = calloc(1, sizeof(stats_t));
    st[0].n_u = c->n;
    st[0].n_q = bit_count(c[0]);
    st[0].sm = strdup(sm);
    st[0].sq = strdup(sq);
    
  } else if (c_mask->fmt <= '1') { // binary mask

    if (c_mask->n != c->n) {
      fprintf(stderr, "[%s:%d] mask (N=%"PRIu64") and query (N=%"PRIu64") are of different lengths.\n", __func__, __LINE__, c_mask->n, c->n);
      fflush(stderr);
      exit(1);
    }
    
    *n_st = 1;
    st = calloc(1, sizeof(stats_t));
    st[0].n_u = c->n;
    st[0].n_q = bit_count(c[0]);
    st[0].n_m = bit_count(c_mask[0]);
    cdata_t tmp = {0};
    tmp.s = malloc((c->n>>3)+1); tmp.n = c->n;
    memcpy(tmp.s, c->s, (c->n>>3)+1);
    for (uint64_t i=0; i<(tmp.n>>3)+1; ++i) tmp.s[i] &= c_mask->s[i];
    st[0].n_o = bit_count(tmp);
    free(tmp.s);
    st[0].sm = strdup(sm);
    st[0].sq = strdup(sq);

  } else if (c_mask->fmt == '6') { // binary mask with universe

    if (c_mask->n != c->n) {
      fprintf(stderr, "[%s:%d] mask (N=%"PRIu64") and query (N=%"PRIu64") are of different lengths.\n", __func__, __LINE__, c_mask->n, c->n);
      fflush(stderr);
      exit(1);
    }
    
    *n_st = 1;
    stats_t st1 = {0};
    for (uint64_t i=0; i<c->n; ++i) {
      if (FMT6_IN_UNI(*c_mask, i)) {
        st1.n_u++;
        int in_q = FMT0_IN_SET(*c, i);
        int in_m = FMT6_IN_SET(*c_mask, i);
        if (in_q) st1.n_q++;
        if (in_m) st1.n_m++;
        if (in_q && in_m) st1.n_o++;
      }
    }
    st = calloc(1, sizeof(stats_t));
    st[0] = st1;
    st[0].sm = strdup(sm);
    st[0].sq = strdup(sq);

  } else if (c_mask->fmt == '2') { // state mask

    if (c_mask->n != c->n) {
      fprintf(stderr, "[%s:%d] mask (N=%"PRIu64") and query (N=%"PRIu64") are of different lengths.\n", __func__, __LINE__, c_mask->n, c->n);
      fflush(stderr);
      exit(1);
    }
    if (!c_mask->aux) fmt2_set_aux(c_mask);
    f2_aux_t *aux = (f2_aux_t*) c_mask->aux;
    *n_st = aux->nk;
    st = calloc((*n_st), sizeof(stats_t));
    uint64_t nq=0;
    for (uint64_t i=0; i<c->n; ++i) {
      uint64_t index = f2_get_uint64(c_mask, i);
      if (index >= (*n_st)) {
        fprintf(stderr, "[%s:%d] State data is corrupted.\n", __func__, __LINE__);
        fflush(stderr);
        exit(1);
      }
      if (FMT0_IN_SET(*c, i)) {
        st[index].n_o++;
        nq++;
      }
      st[index].n_m++;
    }
    for (uint64_t k=0; k < (*n_st); ++k) {
      st[k].n_q = nq;
      st[k].n_u = c->n;
      if (config->section_name) {
        kstring_t tmp = {0};
        ksprintf(&tmp, "%s-%s", sm, aux->keys[k]);
        st[k].sm = tmp.s;
      } else {
        st[k].sm = strdup(aux->keys[k]);
      }
      st[k].sq = strdup(sq);
    }
    
  } else {                      // other masks
    fprintf(stderr, "[%s:%d] Mask format %c unsupported.\n", __func__, __LINE__, c_mask->fmt);
    fflush(stderr);
    exit(1);
  }
  return st;
}

static stats_t* summarize1_queryfmt2(
  cdata_t *c, cdata_t *c_mask, uint64_t *n_st, char *sm, char *sq, config_t *config) {

  stats_t *st = NULL;
  if (c_mask->n == 0) {          // no mask
    
    if (!c->aux) fmt2_set_aux(c);
    f2_aux_t *aux = (f2_aux_t*) c->aux;
    *n_st = aux->nk;
    uint64_t *cnts = calloc(aux->nk, sizeof(uint64_t));
    for (uint64_t i=0; i<c->n; ++i) cnts[f2_get_uint64(c, i)]++;
    st = calloc(aux->nk, sizeof(stats_t));
    for (uint64_t k=0; k<aux->nk; ++k) {
      st[k].n_u = c->n;
      st[k].n_q = cnts[k];
      st[k].n_m = 0;
      st[k].n_o = 0;
      st[k].sm = strdup(sm);
      if (config->section_name) {
        kstring_t tmp = {0};
        ksprintf(&tmp, "%s-%s", sq, aux->keys[k]);
        st[k].sq = tmp.s;
      } else {
        st[k].sq = strdup(aux->keys[k]);
      }
    }
    free(cnts);
    
  } else if (c_mask->fmt <= '1') { // binary mask

    if (!c->aux) fmt2_set_aux(c);
    f2_aux_t *aux = (f2_aux_t*) c->aux;
    *n_st = aux->nk;
    uint64_t *cnts = calloc(aux->nk, sizeof(uint64_t));
    uint64_t *cnts_q = calloc(aux->nk, sizeof(uint64_t));
    uint64_t n_m = 0;
    for (uint64_t i=0; i<c->n; ++i) {
      if (FMT0_IN_SET(*c_mask, i)) {
        n_m++;
        cnts[f2_get_uint64(c, i)]++;
      }
      cnts_q[f2_get_uint64(c, i)]++;
    }
    st = calloc(aux->nk, sizeof(stats_t));
    for (uint64_t k=0; k<aux->nk; ++k) {
      st[k].n_u = c->n;
      st[k].n_q = cnts_q[k];
      st[k].n_o = cnts[k];
      st[k].n_m = n_m;
      st[k].sm = strdup(sm);
      kstring_t tmp = {0};
      ksprintf(&tmp, "%s-%s", sq, aux->keys[k]);
      st[k].sq = tmp.s;
    }
    free(cnts);

  } else if (c_mask->fmt == '6') { // binary mask with universe

    if (!c->aux) fmt2_set_aux(c);
    f2_aux_t *aux = (f2_aux_t*) c->aux;
    *n_st = aux->nk;
    uint64_t *cnts = calloc(aux->nk, sizeof(uint64_t));
    uint64_t *cnts_q = calloc(aux->nk, sizeof(uint64_t));
    uint64_t n_m = 0;
    for (uint64_t i=0; i<c->n; ++i) {
      if (FMT6_IN_UNI(*c_mask,i) && FMT6_IN_SET(*c_mask, i)) {
        n_m++;
        cnts[f2_get_uint64(c, i)]++;
      }
      cnts_q[f2_get_uint64(c, i)]++;
    }
    st = calloc(aux->nk, sizeof(stats_t));
    for (uint64_t k=0; k<aux->nk; ++k) {
      st[k].n_u = c->n;
      st[k].n_q = cnts_q[k];
      st[k].n_o = cnts[k];
      st[k].n_m = n_m;
      st[k].sm = strdup(sm);
      kstring_t tmp = {0};
      ksprintf(&tmp, "%s-%s", sq, aux->keys[k]);
      st[k].sq = tmp.s;
    }
    free(cnts);
    
  } else if (c_mask->fmt == '2') { // state mask

    if (c_mask->n != c->n) {
      fprintf(stderr, "[%s:%d] mask (N=%"PRIu64") and query (N=%"PRIu64") are of different lengths.\n", __func__, __LINE__, c_mask->n, c->n);
      fflush(stderr);
      exit(1);
    }

    if (!c_mask->aux) fmt2_set_aux(c_mask);
    f2_aux_t *aux_m = (f2_aux_t*) c_mask->aux;

    if (!c->aux) fmt2_set_aux(c);
    f2_aux_t *aux_q = (f2_aux_t*) c->aux;

    *n_st = aux_m->nk * aux_q->nk;
    st = calloc((*n_st), sizeof(stats_t));
    uint64_t *nq = calloc(aux_q->nk, sizeof(uint64_t));
    uint64_t *nm = calloc(aux_m->nk, sizeof(uint64_t));
    for (uint64_t i=0; i<c->n; ++i) {
      uint64_t im = f2_get_uint64(c_mask, i);
      uint64_t iq = f2_get_uint64(c, i);
      st[im * aux_q->nk + iq].n_o++;
      nq[iq]++;
      nm[im]++;
    }
    for (uint64_t im=0; im<aux_m->nk; ++im) {
      for (uint64_t iq=0; iq<aux_q->nk; ++iq) {
        stats_t *st1 = &st[im * aux_q->nk + iq];
        st1->n_o++;
        st1->n_u = c->n;
        st1->n_q = nq[iq];
        st1->n_m = nm[im];
        if (config->section_name) {
          kstring_t tmp = {0};
          ksprintf(&tmp, "%s-%s", sm, aux_m->keys[im]);
          st1->sm = tmp.s;
        } else {
          st1->sm = strdup(aux_m->keys[im]);
        }
        if (config->section_name) {
          kstring_t tmp = {0};
          ksprintf(&tmp, "%s-%s", sq, aux_q->keys[iq]);
          st1->sq = tmp.s;
        } else {
          st1->sq = strdup(aux_q->keys[iq]);
        }
      }
    }
    free(nq); free(nm);
    
  } else {                      // other masks
    fprintf(stderr, "[%s:%d] Mask format %c unsupported.\n", __func__, __LINE__, c_mask->fmt);
    fflush(stderr);
    exit(1);
  }
  return st;
}

static stats_t* summarize1_queryfmt3(
  cdata_t *c, cdata_t *c_mask, uint64_t *n_st, char *sm, char *sq, config_t *config) {

  stats_t *st = NULL;
  if (c_mask->n == 0) {            // no mask
    
    *n_st = 1;
    st = calloc(1, sizeof(stats_t));
    st[0].n_u = c->n;
    double sum_beta = 0.0;
    for (uint64_t i=0; i<c->n; ++i) {
      uint64_t mu = f3_get_mu(c, i);
      if (mu) {
        st[0].sum_depth += MU2cov(mu);
        sum_beta += MU2beta(mu);
        st[0].n_o++;
        st[0].n_q++;
      }}
    st[0].sm = strdup(sm);
    st[0].sq = strdup(sq);
    st[0].beta = sum_beta / st[0].n_o; // may have Inf
    
  } else if (c_mask->fmt <= '1') { // binary mask
    
    *n_st = 1;
    st = calloc(1, sizeof(stats_t));
    st[0].n_u = c->n;
    if (c_mask->n != c->n) {
      fprintf(stderr, "[%s:%d] mask (N=%"PRIu64") and query (N=%"PRIu64") are of different lengths.\n", __func__, __LINE__, c_mask->n, c->n);
      fflush(stderr);
      exit(1);
    }
    double sum_beta = 0.0;
    for (uint64_t i=0; i<c->n; ++i) {
      uint64_t mu = f3_get_mu(c, i);
      if (mu) st[0].n_q++;
      if (FMT0_IN_SET(*c_mask, i)) {
        st[0].n_m++;
        if (mu) {
          st[0].sum_depth += MU2cov(mu);
          st[0].sum_beta += MU2beta(mu);
          st[0].n_o++;
        }}}
    st[0].sm = strdup(sm);
    st[0].sq = strdup(sq);
    st[0].beta = sum_beta / st[0].n_o; // may have Inf when n_o == 0

  } else if (c_mask->fmt == '6') { // binary mask with universe
    
    *n_st = 1;
    st = calloc(1, sizeof(stats_t));
    st[0].n_u = c->n;
    if (c_mask->n != c->n) {
      fprintf(stderr, "[%s:%d] mask (N=%"PRIu64") and query (N=%"PRIu64") are of different lengths.\n", __func__, __LINE__, c_mask->n, c->n);
      fflush(stderr);
      exit(1);
    }
    double sum_beta = 0.0;
    for (uint64_t i=0; i<c->n; ++i) {
      uint64_t mu = f3_get_mu(c, i);
      if (mu) st[0].n_q++;
      if (FMT6_IN_UNI(*c_mask, i) && FMT6_IN_SET(*c_mask, i)) {
        st[0].n_m++;
        if (mu) {
          st[0].sum_depth += MU2cov(mu);
          sum_beta += MU2beta(mu);
          st[0].n_o++;
        }}}
    st[0].sm = strdup(sm);
    st[0].sq = strdup(sq);
    st[0].beta = sum_beta / st[0].n_o; // may have Inf when n_o == 0
    
  } else if (c_mask->fmt == '2') { // state mask
    
    if (c_mask->n != c->n) {
      fprintf(stderr, "[%s:%d] mask (N=%"PRIu64") and query (N=%"PRIu64") are of different lengths.\n", __func__, __LINE__, c_mask->n, c->n);
      fflush(stderr);
      exit(1);
    }
    if (!c_mask->aux) fmt2_set_aux(c_mask);
    f2_aux_t *aux = (f2_aux_t*) c_mask->aux;
    *n_st = aux->nk;
    st = calloc((*n_st), sizeof(stats_t));
    uint64_t nq=0;
    for (uint64_t i=0; i<c->n; ++i) {
      uint64_t index = f2_get_uint64(c_mask, i);
      uint64_t mu = f3_get_mu(c, i);
      if (index >= (*n_st)) {
        fprintf(stderr, "[%s:%d] State data is corrupted.\n", __func__, __LINE__);
        fflush(stderr);
        exit(1);
      }
      if (mu) {
        st[index].sum_depth += MU2cov(mu);
        st[index].sum_beta += MU2beta(mu);
        st[index].n_o++;
        nq++;
      }
      st[index].n_m++;
    }
    for (uint64_t k=0; k < (*n_st); ++k) {
      st[k].n_q = nq;
      st[k].n_u = c->n;
      st[k].beta = st[k].sum_beta / st[k].n_o;
      if (config->section_name) {
        kstring_t tmp = {0};
        ksprintf(&tmp, "%s-%s", sm, aux->keys[k]);
        st[k].sm = tmp.s;
      } else {
        st[k].sm = strdup(aux->keys[k]);
      }
      st[k].sq = strdup(sq);
    }
    
  } else {                      // other masks
    fprintf(stderr, "[%s:%d] Mask format %c unsupported.\n", __func__, __LINE__, c_mask->fmt);
    fflush(stderr);
    exit(1);
  }
  return st;
}

static stats_t* summarize1_queryfmt6(
  cdata_t *c, cdata_t *c_mask, uint64_t *n_st, char *sm, char *sq, config_t *config) {

  stats_t *st = NULL;
  if (c_mask->n == 0) {          // no mask
    
    *n_st = 1;
    st = calloc(1, sizeof(stats_t));
    for (uint64_t i=0; i<c->n; ++i) {
      if (FMT6_IN_UNI(*c,i)) {
        st[0].n_u++;
        if (FMT6_IN_SET(*c,i)) st[0].n_q++;
      }
    }
    st[0].sm = strdup(sm);
    st[0].sq = strdup(sq);
    st[0].beta = (double) st[0].n_q / st[0].n_u;
    
  } else if (c_mask->fmt <= '1') { // binary mask

    if (c_mask->n != c->n) {
      fprintf(stderr, "[%s:%d] mask (N=%"PRIu64") and query (N=%"PRIu64") are of different lengths.\n", __func__, __LINE__, c_mask->n, c->n);
      fflush(stderr);
      exit(1);
    }
    
    *n_st = 1;
    st = calloc(1, sizeof(stats_t));
    for (size_t i=0; i<c->n; ++i) {
      if (FMT6_IN_UNI(*c,i)) {
        st[0].n_u++;
        int in_q = FMT6_IN_SET(*c,i);
        int in_m = FMT0_IN_SET(*c_mask,i);
        if (in_q) st[0].n_q++;
        if (in_m) st[0].n_m++;
        if (in_q && in_m) st[0].n_o++;
      }
    }
    st[0].sm = strdup(sm);
    st[0].sq = strdup(sq);
    st[0].beta = (double) st[0].n_o / st[0].n_m;

  } else if (c_mask->fmt == '6') { // binary mask with universe

    if (c_mask->n != c->n) {
      fprintf(stderr, "[%s:%d] mask (N=%"PRIu64") and query (N=%"PRIu64") are of different lengths.\n", __func__, __LINE__, c_mask->n, c->n);
      fflush(stderr);
      exit(1);
    }
    
    *n_st = 1;
    st = calloc(1, sizeof(stats_t));
    for (size_t i=0; i<c->n; ++i) {
      if (FMT6_IN_UNI(*c,i) && FMT6_IN_UNI(*c_mask, i)) {
        st[0].n_u++;
        int in_q = FMT6_IN_SET(*c,i);
        int in_m = FMT6_IN_SET(*c_mask,i);
        if (in_q) st[0].n_q++;
        if (in_m) st[0].n_m++;
        if (in_q && in_m) st[0].n_o++;
      }
    }
    st[0].sm = strdup(sm);
    st[0].sq = strdup(sq);
    st[0].beta = (double) st[0].n_o / st[0].n_m;

  } else if (c_mask->fmt == '2') { // state mask

    if (c_mask->n != c->n) {
      fprintf(stderr, "[%s:%d] mask (N=%"PRIu64") and query (N=%"PRIu64") are of different lengths.\n", __func__, __LINE__, c_mask->n, c->n);
      fflush(stderr);
      exit(1);
    }
    if (!c_mask->aux) fmt2_set_aux(c_mask);
    f2_aux_t *aux = (f2_aux_t*) c_mask->aux;
    *n_st = aux->nk;
    st = calloc((*n_st), sizeof(stats_t));
    uint64_t nq = 0, nu = 0;
    for (uint64_t i=0; i<c->n; ++i) {
      uint64_t index = f2_get_uint64(c_mask, i);
      if (index >= (*n_st)) {
        fprintf(stderr, "[%s:%d] State data is corrupted.\n", __func__, __LINE__);
        fflush(stderr);
        exit(1);
      }
      if (FMT6_IN_UNI(*c,i)) {
        nu++;
        if (FMT6_IN_SET(*c,i)) {
          nq++;
          st[index].n_o++;
        }
        st[index].n_m++;
      }
    }
    for (uint64_t k=0; k < (*n_st); ++k) {
      st[k].n_q = nq;
      st[k].n_u = nu;
      if (config->section_name) {
        kstring_t tmp = {0};
        ksprintf(&tmp, "%s-%s", sm, aux->keys[k]);
        st[k].sm = tmp.s;
      } else {
        st[k].sm = strdup(aux->keys[k]);
      }
      st[k].sq = strdup(sq);
      st[k].beta = (double) st[k].n_o / st[k].n_m;
    }
    
  } else {                      // other masks
    fprintf(stderr, "[%s:%d] Mask format %c unsupported.\n", __func__, __LINE__, c_mask->fmt);
    fflush(stderr);
    exit(1);
  }
  return st;
}

stats_t* summarize1_queryfmt4(cdata_t *c, cdata_t *c_mask, uint64_t *n_st, char *sm, char *sq, config_t *config);

static stats_t* summarize1(cdata_t *c, cdata_t *c_mask, uint64_t *n_st, char *sm, char *sq, config_t *config) {
  if (c->fmt == '3') {
    return summarize1_queryfmt3(c, c_mask, n_st, sm, sq, config);
  } else if (c->fmt == '2') {
    return summarize1_queryfmt2(c, c_mask, n_st, sm, sq, config);
  } else if (c->fmt == '0' || c->fmt == '1') {
    return summarize1_queryfmt0(c, c_mask, n_st, sm, sq, config);
  } else if (c->fmt == '6') {
    return summarize1_queryfmt6(c, c_mask, n_st, sm, sq, config);
  } else if (c->fmt == '4') {
    return summarize1_queryfmt4(c, c_mask, n_st, sm, sq, config);
  } else {
    fprintf(stderr, "[%s:%d] Query format %c unsupported.\n", __func__, __LINE__, c->fmt);
    fflush(stderr);
    exit(1);
  }
}

static void format_stats_and_clean(stats_t *st, uint64_t n_st, const char *fname_qry, config_t *config) {
  const char *fmask = "NA";
  if (!config->full_name) fname_qry = get_basename(fname_qry);
  for (uint64_t i=0; i<n_st; ++i) {
    stats_t s = st[i];
    char *odds_ratio = NULL;
    if (config->fname_mask) {
      double n_mm = s.n_u - s.n_q - s.n_m + s.n_o;
      double n_mp = s.n_q - s.n_o;
      double n_pm = s.n_m - s.n_o;
      kstring_t tmp = {0};
      ksprintf(&tmp, "%1.2f", log2(n_mm*s.n_o / (n_mp*n_pm)));
      odds_ratio = tmp.s;
      if (config->full_name) fmask = config->fname_mask;
      else fmask = get_basename(config->fname_mask);
    } else {
      kstring_t tmp = {0};
      kputs("NA", &tmp);
      odds_ratio = tmp.s;
      fmask = "NA";
    }
    fprintf(stdout,
            "%s\t%s\t%s\t%s\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%s",
            fname_qry, s.sq, fmask, s.sm, s.n_u, s.n_q, s.n_m, s.n_o, odds_ratio);
    if (s.beta >=0) {
      fprintf(stdout, "\t%1.3f", s.beta);
    } else {
      fputs("\tNA", stdout);
    }
    if (s.sum_depth) {
      if (s.n_m) {
        fprintf(stdout, "\t%1.3f", (double) s.sum_depth / s.n_m);
      } else {
        fprintf(stdout, "\t%1.3f", (double) s.sum_depth / s.n_u);
      }
    } else {
      fputs("\tNA", stdout);
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

static void prepare_mask(cdata_t *c) {
  if (c->fmt < '2') {
    convertToFmt0(c);
  } else {
    decompress_in_situ(c);
  }
}

/* The design, first 10 bytes are uint64_t (length) + uint16_t (0=vec; 1=rle) */
int main_summary(int argc, char *argv[]) {
  int c;
  config_t config = {0};
  while ((c = getopt(argc, argv, "m:u:MHFTs:q:h"))>=0) {
    switch (c) {
    case 'm': config.fname_mask = strdup(optarg); break;
    case 'M': config.in_memory = 1; break;
    case 'H': config.no_header = 1; break;
    case 'F': config.full_name = 1; break;
    case 'T': config.section_name = 1; break;
    case 's': config.fname_snames = strdup(optarg); break;
    case 'q': config.fname_qry_stdin = optarg; break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 1 > argc) { 
    usage(); 
    wzfatal("Please supply input file.\n"); 
  }

  cfile_t cf_mask; int unseekable = 0;
  snames_t snames_mask = {0};
  cdata_t *c_masks = NULL; uint64_t c_masks_n = 0;
  if (config.fname_mask) {
    cf_mask = open_cfile(config.fname_mask);
    unseekable = bgzf_seek(cf_mask.fh, 0, SEEK_SET);
    snames_mask = loadSampleNamesFromIndex(config.fname_mask);
  }
  
  if (config.in_memory || unseekable) { /* load in-memory masks */
    c_masks = calloc(1, sizeof(cdata_t));
    c_masks_n = 0;
    for (;;++c_masks_n) {
      cdata_t c_mask = read_cdata1(&cf_mask);
      if (c_mask.n == 0) break;
      prepare_mask(&c_mask);
      c_masks = realloc(c_masks, (c_masks_n+1)*sizeof(cdata_t));
      c_masks[c_masks_n] = c_mask;
    }
  }
  
  if (!config.no_header) {
    fputs("QFile\tQuery\tMFile\tMask\tN_univ\tN_query\tN_mask\tN_overlap\tLog2OddsRatio\tBeta\tDepth\n", stdout);
  }

  for (int j = optind; j < argc; ++j) {
    char *fname_qry = argv[j];
    cfile_t cf_qry = open_cfile(fname_qry);
    snames_t snames_qry = {0};
  
    if (config.fname_snames) snames_qry = loadSampleNames(config.fname_snames, 1);
    else snames_qry = loadSampleNamesFromIndex(fname_qry);

    if (strcmp(fname_qry, "-")==0 && config.fname_qry_stdin)
      fname_qry = config.fname_qry_stdin;

    for (uint64_t kq=0;;++kq) {
      cdata_t c_qry = read_cdata1(&cf_qry);
      if (c_qry.n == 0) break;
      if (snames_qry.n && kq >= (unsigned) snames_qry.n) {
        fprintf(stderr, "[%s:%d] More data (N=%"PRIu64") found than specified in the index file (N=%d).\n", __func__, __LINE__, kq+1, snames_qry.n);
        fflush(stderr);
        exit(1);
      }
      if (c_qry.fmt == '7') { // skip format 7
        free_cdata(&c_qry); c_qry.s = NULL;
        continue;
      }
      kstring_t sq = {0};
      if (snames_qry.n) kputs(snames_qry.s[kq], &sq);
      else ksprintf(&sq, "%"PRIu64"", kq+1);
      prepare_mask(&c_qry);

      if (config.fname_mask) {   /* apply any mask? */
        if (c_masks_n) {        /* in memory or unseekable */
          for (uint64_t km=0;km<c_masks_n;++km) {
            cdata_t c_mask = c_masks[km];
            kstring_t sm = {0};
            if (snames_mask.n) kputs(snames_mask.s[km], &sm);
            else ksprintf(&sm, "%"PRIu64"", km+1);
            uint64_t n_st = 0;
            stats_t *st = summarize1(&c_qry, &c_mask, &n_st, sm.s, sq.s, &config);
            format_stats_and_clean(st, n_st, fname_qry, &config);
            free(sm.s);
          }
        } else {                /* mask is seekable */
          if (bgzf_seek(cf_mask.fh, 0, SEEK_SET)!=0) {
            fprintf(stderr, "[%s:%d] Cannot seek mask.\n", __func__, __LINE__);
            fflush(stderr);
            exit(1);
          }
          for (uint64_t km=0;;++km) {
            cdata_t c_mask = read_cdata1(&cf_mask);
            if (c_mask.n == 0) break;
            prepare_mask(&c_mask);

            kstring_t sm = {0};
            if (snames_mask.n) kputs(snames_mask.s[km], &sm);
            else ksprintf(&sm, "%"PRIu64"", km+1);
            uint64_t n_st = 0;
            stats_t *st = summarize1(&c_qry, &c_mask, &n_st, sm.s, sq.s, &config);
            format_stats_and_clean(st, n_st, fname_qry, &config);
            free(sm.s);
            free_cdata(&c_mask);
          }
        }
      } else {                  /* whole dataset summary if missing mask */
        kstring_t sm = {0}; cdata_t c_mask = {0};
        kputs("global", &sm);
        uint64_t n_st = 0;
        stats_t *st = summarize1(&c_qry, &c_mask, &n_st, sm.s, sq.s, &config);
        format_stats_and_clean(st, n_st, fname_qry, &config);
        free(sm.s);
      }
      free(sq.s);
      free_cdata(&c_qry); c_qry.s = NULL;
    }
    if (c_masks_n) {
      for (uint64_t i=0; i<c_masks_n; ++i) free_cdata(&c_masks[i]);
      free(c_masks);
    }
    bgzf_close(cf_qry.fh);
    cleanSampleNames2(snames_qry);
  }
  if (config.fname_snames) free(config.fname_snames);
  if (config.fname_mask) bgzf_close(cf_mask.fh);
  if (config.fname_mask) free(config.fname_mask);
  cleanSampleNames2(snames_mask);
  
  return 0;
}
