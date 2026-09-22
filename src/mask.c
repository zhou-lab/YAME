// SPDX-License-Identifier: LicenseRef-CHOP-Academic-BSD-2-Clause
/**
 * This file is part of YAME.
 *
 * Copyright (C) 2021-present The Children's Hospital of Philadelphia
 *
 * Use of this software is available to academic and non-profit institutions
 * for research purposes under the 2-Clause BSD License; for use or transfers
 * to commercial entities, inquire with Dr. Wanding Zhou at zhouw3@chop.edu.
 * See the LICENSE file at the root of the repository for the full terms.
 */

#include <sys/stat.h>
#include "yame_ui.h"
#include <sys/types.h>
#include "cfile.h"
#include "assets.h"

static int usage() {
  yame_usage_head("yame mask [options] <in.cg> <mask>");
  yame_usage_text("Blank the sites the mask covers; every other site is left as it is.");
  yame_usage_text("`mask x.cg Blacklist.cm` REMOVES the blacklisted sites; -v keeps only");
  yame_usage_text("them. The mask may be format 0, 1, 3 (covered = M+U > 0) or 6 (covered");
  yame_usage_text("= in the universe), so a methylome can mask another: `mask truth.cg");
  yame_usage_text("query.cg` blanks the sites the query was shown. A bare name (CGI,");
  yame_usage_text("Blacklist) resolves in the store for the input's row space, as -m does");
  yame_usage_text("for summary. Row counts must match.");
  yame_usage_sec("Options:");
  yame_usage_opt("-o", "output cx file name. if missing, output to stdout without index.");
  yame_usage_opt("-c", "contextualize binary input to format 6 using '1's in mask.");
  yame_usage_cont("implicit for format 6 input (output is always format 6).");
  yame_usage_opt("-v", "invert the mask: blank the sites it does NOT cover.");
  yame_usage_opt("-h", "This help");
  fprintf(stderr, "\n");

  return 1;
}

void mask_fmt3(cdata_t *c, cdata_t c_mask, BGZF *fp_out) {
  for (uint64_t i=0; i<c->n; ++i) {
    if (FMT0_IN_SET(c_mask, i)) {
      f3_set_mu(c, i, 0, 0);
    }
  }
  cdata_compress(c);
  bgzf_flush(fp_out);           /* start this record on a block boundary */
  cdata_write1(fp_out, c);
}

void mask_fmt0(cdata_t *c, cdata_t c_mask, BGZF *fp_out) {
  for (uint64_t i=0; i<cdata_nbytes(c); ++i) {
    c->s[i] &= ~c_mask.s[i];
  }
  /* cdata_compress(&c); */
  bgzf_flush(fp_out);           /* start this record on a block boundary */
  cdata_write1(fp_out, c);
}

/* Blank the sites the mask covers -- the same sense as the fmt3 and fmt0
 * paths and as -h says. This branch used to do the opposite, keeping only
 * the covered sites, so `mask -o clean.cg x.cg Blacklist.cm` on a fmt6 file
 * kept the blacklist and blanked everything else (reported 2026-09-22). */
void mask_fmt6(cdata_t *c, cdata_t c_mask, BGZF *fp_out) {
  for (uint64_t i = 0; i < c->n; ++i) {
    if (FMT6_IN_UNI(*c, i) && FMT0_IN_SET(c_mask, i)) {
      FMT6_SET_NA(*c, i);
    }
  }
  cdata_compress(c);
  bgzf_flush(fp_out);           /* start this record on a block boundary */
  cdata_write1(fp_out, c);
}

void fmt0ContextualizeFmt6(cdata_t *c, cdata_t c_mask, BGZF *fp_out) {
  cdata_t c6 = {.fmt = '6', .n = c->n};
  c6.s = wzcalloc((c6.n+3)/4, sizeof(uint8_t));
  for (uint64_t i=0; i<c6.n; ++i) {
    if (FMT0_IN_SET(c_mask,i)) { // mask is used as universe, use -v to invert
      if (FMT0_IN_SET(*c, i)) FMT6_SET1(c6, i);
      else FMT6_SET0(c6, i);
    }
  }
  cdata_compress(&c6);
  bgzf_flush(fp_out);           /* start this record on a block boundary */
  cdata_write1(fp_out, &c6);
  free_cdata(&c6);
}

int main_mask(int argc, char *argv[]) {

  int c, reverse = 0, contextualize_to_fmt6 = 0;
  char *fname_out = NULL;
  while ((c = getopt(argc, argv, "o:cvh"))>=0) {
    switch (c) {
    case 'o': fname_out = wzstrdup(optarg); break;
    case 'c': contextualize_to_fmt6 = 1; break;
    case 'v': reverse = 1; break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 2 > argc) {
    usage(); 
    wzfatal("Please supply input file.\n");
  }

  char *fname = argv[optind];
  char *fname_mask = argv[optind+1];

  /* A bare name resolves in the store for the input's row space, the same
   * way summary -m does; an existing path is used as given. Before this,
   * `mask x.cg Blacklist` died with "Error opening file Blacklist" while
   * summary answered the same mistake with where to look. */
  char resolved[4096];
  const char *rname = NULL, *rfetch = NULL;
  {
    uint64_t rows = yame_ref_file_rows(fname);
    int st = yame_ref_resolve(fname_mask, rows, NULL, NULL, resolved,
                              sizeof(resolved), &rname, &rfetch);
    if (st != YAME_REF_OK) {
      yame_ref_explain_name(stderr, fname_mask, rows, st, rname, rfetch, "mask");
      return 1;
    }
    if (strcmp(resolved, fname_mask) != 0)
      fprintf(stderr, "[mask] %s -> %s\n", fname_mask, resolved);
  }

  cfile_t cf_mask = open_cfile(resolved);
  cdata_t c_mask = read_cdata1(&cf_mask);
  if (c_mask.fmt == '1' || c_mask.fmt == '3' || c_mask.fmt == '6') convertToFmt0(&c_mask);
  if (c_mask.fmt != '0')
    wzfatal("mask: %s is format %c; a mask must be format 0, 1, 3 or 6.\n",
            resolved, c_mask.fmt);
  if (reverse) {
    for (uint64_t i=0; i<cdata_nbytes(&c_mask); ++i) {
      c_mask.s[i] = ~(c_mask.s[i]);
    }
  }

  BGZF *fp_out;
  if (fname_out) fp_out = bgzf_open2(fname_out, "w");
  else fp_out = yame_bgzf_stdout("w", "mask");
  if (fp_out == NULL) {
    fprintf(stderr, "Error opening file for writing: %s\n", fname_out);
    exit(1);
  }

  cfile_t cf = open_cfile(fname);
  while (1) {
    cdata_t c = read_cdata1(&cf);
    if (c.n == 0) break;
    if (c.fmt == '1') convertToFmt0(&c);
    decompress_in_situ(&c);
    if (c.n != c_mask.n) {
      fprintf(stderr, "[%s:%d] mask (n=%"PRIu64") and query (N=%"PRIu64") are of different lengths.\n", __func__, __LINE__, c_mask.n, c.n);
      fflush(stderr);
      exit(1);
    }

    if (c.fmt == '3') {
      mask_fmt3(&c, c_mask, fp_out);
    } else if (c.fmt == '0') {
      if (contextualize_to_fmt6) {
        fmt0ContextualizeFmt6(&c, c_mask, fp_out);
      } else {
        mask_fmt0(&c, c_mask, fp_out);
      }
    } else if (c.fmt == '6') {
      mask_fmt6(&c, c_mask, fp_out);
    } else {
      fprintf(stderr, "[%s:%d] Only format %d files are supported.\n", __func__, __LINE__, c.fmt);
      fflush(stderr);
      exit(1);
    }
    free_cdata(&c);
  }

  if (fname_out) free(fname_out);
  bgzf_close(fp_out);
  free_cdata(&c_mask);
  bgzf_close(cf.fh);
  bgzf_close(cf_mask.fh);
  return 0;
}
