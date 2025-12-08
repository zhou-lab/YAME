// SPDX-License-Identifier: AGPL-3.0-or-later
/**
 * This file is part of YAME.
 *
 * Copyright (C) 2025 Wanding Zhou
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


#ifndef _WZBED_H
#define _WZBED_H

#include <inttypes.h>
#include "wztarget.h"
#include "wzio.h"

/* This file defines several bed file parsers */

/*************************************
 ** General-purpose bed file parser **
 *************************************
 bed1_v *beds = bed_read(file_path);
**/

/****************
 ** Bed Record **
 ****************/
typedef struct bed1_t {
  /* int tid; */
  char *seqname;                /* points to bed_file_t.seqname, does not own! */
  int64_t beg;
  int64_t end;
  void *data;
} bed1_t;

DEFINE_VECTOR(bed1_v, bed1_t)

typedef void (*init_data_f)(bed1_t *b, void *aux_data);
typedef void (*parse_data_f)(bed1_t *b, char **fields, int nfields);
typedef void (*free_data_f)(void *data);

static inline bed1_t *init_bed1(init_data_f init_data, void *aux_data) {
  bed1_t *b = calloc(1, sizeof(bed1_t));
  if (init_data != NULL)
    init_data(b, aux_data);
  return b;
}

static inline void free_bed1(bed1_t *b, free_data_f free_data) {
  if (free_data != NULL)
    free_data(b->data);
  free(b);
}

/**************
 ** Bed File **
 **************/

typedef struct bed_file_t {
  char *file_path;
  gzFile fh;
  char *line;
  char *seqname;
  /* target_v *targets; */
  targets_t *targets;
} bed_file_t;

static inline bed_file_t *init_bed_file(char *file_path) {
  bed_file_t *bed = calloc(1, sizeof(bed_file_t));
  bed->file_path = strdup(file_path);
  bed->fh = wzopen(bed->file_path, 1);
  bed->targets = 0;
  bed->line = NULL;
  bed->seqname = NULL;
  return bed;
}

static inline void free_bed_file(bed_file_t *bed) {
  wzclose(bed->fh);
  /* destroy_target_v(bed->targets); */
  if (bed->targets) destroy_targets(bed->targets);
  free(bed->file_path);
  free(bed->line);
  free(bed->seqname);
  free(bed);
}

static inline int bed_read1(bed_file_t *bed, bed1_t *b, parse_data_f parse_data) {
  if (bed->fh == NULL) return 0;
  if (gzFile_read_line(bed->fh, &bed->line) == 0) return 0;

  char **fields; int nfields;
  line_get_fields(bed->line, "\t", &fields, &nfields);
  if (nfields < 3)
    wzfatal("[%s:%d] Bed file has fewer than 3 columns.\n", __func__, __LINE__);

  if (bed->seqname==NULL || strcmp(fields[0], bed->seqname) != 0) {
    free(bed->seqname); bed->seqname = strdup(fields[0]);
  }
  b->seqname = bed->seqname;
  /* b->tid = get_tid(bed->targets, fields[0], 1); */
  /* locate_or_insert_target_v(bed->targets, fields[0]); */

  ensure_number(fields[1]);
  b->beg = atoi(fields[1]);

  ensure_number(fields[2]);
  b->end = atoi(fields[2]);

  if (parse_data != NULL) parse_data(b, fields, nfields);
  else b->data = NULL;
  free_fields(fields, nfields);

  return 1;
}

#endif /* _WZBED_H */
