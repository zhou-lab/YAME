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

#include <time.h>
#include "yame_ui.h"
#include "cfile.h"

static int usage(void) {
  yame_usage_head("yame perturb [options] <in.cx>");
  yame_usage_text("Randomly flip 0/1 bits for format 0 and format 6.");
  yame_usage_text("- For format 0, each set bit (1) or unset bit (0) is independently");
  yame_usage_text("  flipped with probability p.");
  yame_usage_text("- For format 6, only in-universe sites are eligible; their set bit");
  yame_usage_text("  is flipped with probability p. NA sites are left unchanged.");
  yame_usage_sec("Options:");
  yame_usage_opt("-s [int]", "random seed (default: current time).");
  yame_usage_opt("-p [float]", "fraction of CpGs to flip, in [0,1] (default: 0.05).");
  yame_usage_opt("-o [PATH]", "output .cx file (default: stdout, no index written).");
  yame_usage_opt("-h", "this help.");
  fprintf(stderr, "\n");

  return 1;
}

static void perturb_fmt0(cdata_t *c, double p) {
  for (uint64_t i = 0; i < c->n; ++i) {
    if ((double)rand() / ((double)RAND_MAX + 1.0) < p) {
      /* flip bit i */
      c->s[i >> 3] ^= (1u << (i & 0x7));
    }
  }
}

static void perturb_fmt6(cdata_t *c, double p) {
  for (uint64_t i = 0; i < c->n; ++i) {
    if (!FMT6_IN_UNI(*c, i)) continue;
    if ((double)rand() / ((double)RAND_MAX + 1.0) < p) {
      if (FMT6_IN_SET(*c, i)) FMT6_SET0(*c, i);
      else FMT6_SET1(*c, i);
    }
  }
}

int main_perturb(int argc, char *argv[]) {
  int c;
  unsigned seed = (unsigned)time(NULL);
  double p = 0.05;
  char *fname_out = NULL;

  while ((c = getopt(argc, argv, "s:p:o:h")) >= 0) {
    switch (c) {
    case 's': seed = (unsigned)strtoul(optarg, NULL, 10); break;
    case 'p': p = atof(optarg); break;
    case 'o': fname_out = strdup(optarg); break;
    case 'h': return usage();
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (p < 0.0 || p > 1.0) wzfatal("Flip probability -p must be in [0, 1].\n");
  if (optind + 1 > argc) { usage(); wzfatal("Please supply input file.\n"); }

  char *fname = argv[optind];
  srand(seed);

  BGZF *fp_out;
  if (fname_out) fp_out = bgzf_open2(fname_out, "wb");
  else           fp_out = bgzf_dopen(fileno(stdout), "wb");
  if (!fp_out) {
    fprintf(stderr, "[%s:%d] Error opening output: %s\n",
            __func__, __LINE__, fname_out ? fname_out : "<stdout>");
    exit(1);
  }

  cfile_t cf = open_cfile(fname);
  for (;;) {
    cdata_t cin = read_cdata1(&cf);
    if (cin.n == 0) { free_cdata(&cin); break; }
    decompress_in_situ(&cin);

    switch (cin.fmt) {
    case '0': perturb_fmt0(&cin, p); break;
    case '6': perturb_fmt6(&cin, p); break;
    default:
      wzfatal("Format '%c' not supported by perturb (only 0 and 6).\n", cin.fmt);
    }

    cdata_compress(&cin);
    bgzf_flush(fp_out);         /* start this record on a block boundary */
    cdata_write1(fp_out, &cin);
    free_cdata(&cin);
  }

  bgzf_close(cf.fh);
  bgzf_close(fp_out);
  if (fname_out) free(fname_out);
  return 0;
}
