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

#ifndef _SNAMES_H
#define _SNAMES_H

#include <zlib.h>
#include "wzio.h"

typedef struct {
  int n;
  char **s;
} snames_t;

/**
 * Loads the sample names from a given file.
 * If the filename is "-", it will read from stdin. 
 * If the file cannot be opened and the fatal parameter is non-zero, the program will exit with an error message.
 *
 * @param fname_snames The name of the file containing the sample names.
 * @param fatal A flag indicating whether to treat inability to open the file as a fatal error.
 *              If non-zero and the file cannot be opened, the program will exit with an error message.
 * @return A snames_t structure containing the sample names read from the file.
 *         If the file cannot be opened and the fatal parameter is zero, NULL is returned.
 *         If there is an error allocating memory, NULL is returned.
 */
static inline snames_t loadSampleNames(char* fname_snames, int fatal) {
  snames_t snames = {0};
  if (fname_snames == NULL) return snames;
  gzFile fp;
  if (strcmp(fname_snames, "-") == 0) {
    fp = gzdopen(fileno(stdin), "r");
  } else {
    fp = gzopen(fname_snames, "r");
    if (!fp && fatal) {
      fprintf(stderr, "[%s:%d] Fatal, cannot open file: %s\n",
              __func__, __LINE__, fname_snames);
      fflush(stderr);
      exit(1);
    }
  }
  
  if (fp == NULL) return snames;

  char *line = NULL;
  while (gzFile_read_line(fp, &line) > 0) {
    char *sname;
    if (line_get_field(line, 0, "\t", &sname)) {
      snames.s = realloc(snames.s, sizeof(*(snames.s)) * (snames.n + 1));
      if (snames.s == NULL) {
        fprintf(stderr, "Failed to allocate memory\n");
        fflush(stderr);
        exit(1);
      }
      snames.s[snames.n] = sname;
      snames.n++;
    }
  }
  free(line);
  line = NULL;
  gzclose(fp);

  return snames;
}

static inline void cleanSampleNames2(snames_t snames) {
  if (snames.n)
    for (int i=0; i< snames.n; ++i) {
      free(snames.s[i]);
    }
  free(snames.s);
}

static inline void cleanSampleNames(snames_t *snames) {
  if (snames) {
    if (snames->n)
      for (int i=0; i< snames->n; ++i) {
        free(snames->s[i]);
      }
    free(snames->s);
  }
}

#endif  /* SNAMES_H */
