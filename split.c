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

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: yame split [options] <in.cx> out_prefix\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -v        verbose\n");
  fprintf(stderr, "    -s        sample name list\n");
  fprintf(stderr, "    -h        This help\n");
  fprintf(stderr, "\n");

  return 1;
}

int main_split(int argc, char *argv[]) {

  int c, verbose = 0; char *fname_snames = NULL;
  while ((c = getopt(argc, argv, "s:vh"))>=0) {
    switch (c) {
    case 'v': verbose = 1; break;
    case 's': fname_snames = strdup(optarg); break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 2 > argc) { 
    usage(); 
    wzfatal("Please supply input file.\n");
  }

  cfile_t cf = open_cfile(argv[optind++]);
  char *prefix = argv[optind];

  char **snames = NULL; int snames_n = 0;
  if (fname_snames) {
    gzFile fh = wzopen(fname_snames, 1);
    char *line = NULL;
    char **fields; int nfields;
    while (gzFile_read_line(fh, &line)>0) {
      line_get_fields(line, "\t", &fields, &nfields);
      snames = realloc(snames, (snames_n+1)*sizeof(char*));
      snames[snames_n] = strdup(fields[0]);
      ++snames_n;
      free_fields(fields, nfields);
    }
    free(line);
    wzclose(fh);
    free(fname_snames);
  }
  
  int i = 0;
  for (i=0; ; ++i) {
    cdata_t c = read_cdata1(&cf);
    if (c.n == 0) break;
    char *tmp = NULL;
    if (snames_n) {
      tmp = malloc(strlen(prefix)+strlen(snames[i])+1000);
      sprintf(tmp, "%s%s.cx", prefix, snames[i]);
    } else {
      tmp = malloc(strlen(prefix) + 1000);
      sprintf(tmp, "%s_split_%i.cx", prefix, i+1);
    }
    cdata_write(tmp, &c, "wb", verbose);
    free(tmp);
    free(c.s);
  }
  
  return 0;
}
