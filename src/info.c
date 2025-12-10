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

#include "cfile.h"
#include "kstring.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: yame info [options] <in.cx>\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -1        Report one record per file.\n");
  fprintf(stderr, "    -h        This help\n");
  fprintf(stderr, "\n");

  return 1;
}

static void cdata_length(cdata_t *c, uint64_t *n, uint8_t *u) {
  switch(c->fmt) {
  case '7': {
    *n = fmt7_data_length(c);
    *u = 1;
    break; }
  default: {
    cdata_t inflated = decompress(*c);
    *n = inflated.n;
    *u = inflated.unit;
    free(inflated.s);
    break; }
  }
}

int main_info(int argc, char *argv[]) {

  int c; int report1 = 0;
  while ((c = getopt(argc, argv, "1hv"))>=0) {
    switch (c) {
    case '1': report1 = 1; break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 1 > argc) { 
    usage(); 
    wzfatal("Please supply input file.\n"); 
  }

  fprintf(stdout, "File\tSample\tNcol\tNrow\tFormat\tUnitBytes\tKeys\n");
  for (int j = optind; j < argc; ++j) {
    char *fname_in = argv[j];
    cfile_t cf = open_cfile(fname_in);
    snames_t snames = loadSampleNamesFromIndex(fname_in);
    int i = 0;
    for (i=0; ; ++i) {
      cdata_t c = read_cdata1(&cf);
      if (c.n == 0) break;
      if (snames.n && i >= snames.n) {
        fprintf(stderr, "[%s:%d] More data (N=%d) found than specified in the index file (N=%d).\n", __func__, __LINE__, i+1, snames.n);
        fflush(stderr);
        exit(1);
      }

      fprintf(stdout, "%s\t", fname_in);
      if (snames.n) {
        fputs(snames.s[i], stdout);
      } else {
        fprintf(stdout, "%d", i+1);
      }
      if (snames.n) fprintf(stdout, "\t%d", snames.n);
      else fputs("\tNA", stdout);
      uint64_t length = 0; uint8_t unit = 0;
      cdata_length(&c, &length, &unit);
      fprintf(stdout, "\t%"PRIu64"\t%c\t%u\t", length, c.fmt, unit);

      kstring_t tmp = {0};
      if (c.fmt == '2') {
        uint64_t nk = fmt2_get_keys_n(&c);
        ksprintf(&tmp, "N=%"PRIu64"|", nk);
        const char *s = (char*) c.s;
        for (uint64_t k=0; k<nk; ++k) {
          if (k) kputc(',', &tmp);
          kputs(s, &tmp);
          s += (strlen(s)+1);
          if (k+1 < nk && strlen(tmp.s) > 25) {
            kputs(",...", &tmp);
            break;
          }
        }
      } else {
        kputs("NA", &tmp);
      }
      fputs(tmp.s, stdout);
      free(tmp.s);
      fputc('\n', stdout);
      
      free_cdata(&c);
      if (report1) break;
    }
    cleanSampleNames2(snames);
    bgzf_close(cf.fh);
  }
  return 0;
}
