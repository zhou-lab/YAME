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

/**
 * yame summary
 * ============
 *
 * Goal
 * ----
 * Produce per-sample summary statistics for a query CX record, optionally
 * against one or more mask CX records. This is designed for fast “feature-set”
 * overlap/enrichment reporting on packed methylation formats. :contentReference[oaicite:1]{index=1}
 *
 * Inputs
 * ------
 * - Query file(s): one or more .cx streams. Each record is treated as one query
 *   sample. Query formats supported are routed by summarize1():
 *     fmt 0/1: binary presence/absence (fmt 1 is treated like binary query)
 *     fmt 2:   categorical/state data (keys + integer-coded states)
 *     fmt 3:   MU counts (depth and beta available)
 *     fmt 4:   float vector (summary handler may compute beta-like metrics)
 *     fmt 6:   2-bit packed (set + universe) summary
 *     fmt 7:   rows/coordinates (summary handler may be implemented separately)
 *
 * - Mask file (optional, -m): may contain multiple records (mask samples).
 *   For each query sample, we compute summaries against every mask sample.
 *
 * Mask loading strategy
 * ---------------------
 * If the mask BGZF stream is seekable, we re-seek to the beginning for each
 * query sample (lowest memory). If it is unseekable, or if -M is requested,
 * all mask records are read and prepared once into RAM. :contentReference[oaicite:2]{index=2}
 *
 * Preparation of query/mask records
 * ---------------------------------
 * prepare_mask() normalizes “mask-like operations”:
 *   - fmt 0/1 are converted to fmt 0 bitset (convertToFmt0)
 *   - fmt >= 2 are decompressed in-place (decompress2)
 * This ensures summarize1_* can assume consistent in-memory representation.
 *
 * What gets reported (one output row per reported category)
 * ---------------------------------------------------------
 * The output row schema is:
 *   QFile, Query, MFile, Mask, N_univ, N_query, N_mask, N_overlap,
 *   Log2OddsRatio, Beta, Depth
 *
 * Definitions:
 * - N_univ:   the universe size used for the test (may be constrained by fmt6 universe).
 * - N_query:  size of query set (or per-state count for fmt2 queries).
 * - N_mask:   size of mask set (or per-state count for fmt2 masks).
 * - N_overlap:|query ∩ mask| (or per-state overlap counts).
 * - Log2OddsRatio:
 *     log2( (n_mm * n_pp) / (n_mp * n_pm) ) where:
 *       n_pp = N_overlap
 *       n_mp = N_query - N_overlap
 *       n_pm = N_mask  - N_overlap
 *       n_mm = N_univ - N_query - N_mask + N_overlap
 *   If no mask is provided, this column prints NA.
 * - Beta:
 *   For MU-based summaries (fmt3 / fmt6), Beta is the mean methylation (or fraction)
 *   as implemented in the corresponding summarize1_queryfmt*.
 * - Depth:
 *   When available, reports mean depth over masked sites (or over universe when no mask).
 *
 * State (fmt2) naming
 * -------------------
 * For fmt2 data, rows may be emitted per state key. With -T, the state key is
 * appended/prefixed so outputs remain unique and interpretable even when
 * state names repeat across sections.
 *
 * Stdin naming
 * ------------
 * When the query filename is '-' (stdin), -q supplies a human-readable name
 * for the QFile column (otherwise it would be '-' and not traceable).
 */

static int usage(void) {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "  yame summary [options] <query.cx> [query2.cx ...]\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Purpose:\n");
  fprintf(stderr, "  Summarize a query feature set (or per-state composition) and optionally\n");
  fprintf(stderr, "  its overlap/enrichment against one or more masks.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Input:\n");
  fprintf(stderr, "  <query.cx> may contain one or multiple samples (records). Supported query\n");
  fprintf(stderr, "  formats: 0/1 (binary), 2 (state), 3 (MU counts), 4 (float),\n");
  fprintf(stderr, "           6 (set+universe), 7 (genomic coordinates).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Masking:\n");
  fprintf(stderr, "  -m <mask.cx>   Optional mask feature file (can be multi-sample).\n");
  fprintf(stderr, "                 If provided, every query sample is summarized against every\n");
  fprintf(stderr, "                 mask sample (cartesian product).\n");
  fprintf(stderr, "  -M             Load all masks into memory (faster when mask file is on slow IO).\n");
  fprintf(stderr, "                 Also auto-enabled when the mask stream is unseekable.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Naming / output formatting:\n");
  fprintf(stderr, "  -H             Suppress the header line.\n");
  fprintf(stderr, "  -F             Use full paths in QFile/MFile (default: basename only).\n");
  fprintf(stderr, "  -T             Always include section/state names in output labels when\n");
  fprintf(stderr, "                 summarizing format-2 (state) data.\n");
  fprintf(stderr, "  -s <list.txt>  Override query sample names using a plain-text list.\n");
  fprintf(stderr, "                 Only applies to the first query file.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Stdin helpers:\n");
  fprintf(stderr, "  -q <name>      Backup query file name used only when <query.cx> is '-'.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Other:\n");
  fprintf(stderr, "  -6             Treat format-6 query as 2bit quaternary than set/universe.\n");
  fprintf(stderr, "  -h             Show this help message.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Output columns:\n");
  fprintf(stderr, "  QFile  Query  MFile  Mask  N_univ  N_query  N_mask  N_overlap  Log2OddsRatio  Beta  Depth\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Notes:\n");
  fprintf(stderr, "  * For state masks (format 2), summary is emitted per state key (one row per key).\n");
  fprintf(stderr, "  * When no mask is given, Mask is reported as 'global'.\n");
  fprintf(stderr, "\n");
  return 1;
}

/* fprintf(stderr, "    -u        Optional universe set as a .cx file. If given, the masks and queries are both subset.\n"); */

stats_t* summarize1_queryfmt0(cdata_t *c, cdata_t *c_mask, uint64_t *n_st, char *sm, char *sq, config_t *config);
stats_t* summarize1_queryfmt2(cdata_t *c, cdata_t *c_mask, uint64_t *n_st, char *sm, char *sq, config_t *config);
stats_t* summarize1_queryfmt3(cdata_t *c, cdata_t *c_mask, uint64_t *n_st, char *sm, char *sq, config_t *config);
stats_t* summarize1_queryfmt4(cdata_t *c, cdata_t *c_mask, uint64_t *n_st, char *sm, char *sq, config_t *config);
stats_t* summarize1_queryfmt6(cdata_t *c, cdata_t *c_mask, uint64_t *n_st, char *sm, char *sq, config_t *config);
stats_t* summarize1_queryfmt7(cdata_t *c, cdata_t *c_mask, uint64_t *n_st, char *sm, char *sq, config_t *config);

static stats_t* summarize1(cdata_t *c, cdata_t *c_mask, uint64_t *n_st, char *sm, char *sq, config_t *config) {

  switch (c->fmt) {
  case '0': return summarize1_queryfmt0(c, c_mask, n_st, sm, sq, config);
  case '1': return summarize1_queryfmt0(c, c_mask, n_st, sm, sq, config);
  case '2': return summarize1_queryfmt2(c, c_mask, n_st, sm, sq, config);
  case '3': return summarize1_queryfmt3(c, c_mask, n_st, sm, sq, config);
  case '4': return summarize1_queryfmt4(c, c_mask, n_st, sm, sq, config);
  case '6': return summarize1_queryfmt6(c, c_mask, n_st, sm, sq, config);
  case '7': return summarize1_queryfmt7(c, c_mask, n_st, sm, sq, config);
  default: wzfatal("[%s:%d] Query format %c unsupported.\n", __func__, __LINE__, c->fmt);
  }
  stats_t *st = calloc(1, sizeof(stats_t));
  return st;
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
  while ((c = getopt(argc, argv, "m:u:MHFTs:6q:h"))>=0) {
    switch (c) {
    case 'm': config.fname_mask = strdup(optarg); break;
    case 'M': config.in_memory = 1; break;
    case '6': config.f6_as_2bit = 1; break;
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
      /* if (c_qry.fmt == '7') { // skip format 7 */
      /*   free_cdata(&c_qry); c_qry.s = NULL; */
      /*   continue; */
      /* } */
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
