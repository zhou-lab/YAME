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

void decompress(cdata_t *c, cdata_t *expanded) {
  switch (c->fmt) {
  case '0': { fmt0_decompress(c, expanded); break; }
  case '1': { fmt1_decompress(c, expanded); break; }
  case '2': { fmt2_decompress(c, expanded); break; }
  case '3': { fmt3_decompress(c, expanded); break; }
  case '4': { fmt4_decompress(c, expanded); break; }
  case '5': { fmt5_decompress(c, expanded); break; }
  case '6': { fmt6_decompress(c, expanded); break; }
  default: wzfatal("Unsupported format for inflation: %c.\n", c->fmt);
  }
  /* shouldn't reach here */
}

/* decompress in situ */
void decompress2(cdata_t *c) {
  if (!c->compressed) {
    fprintf(stderr, "[%s:%d] Already decompressed.\n", __func__, __LINE__);
    fflush(stderr);
    exit(1);
  }
  cdata_t expanded = *c;
  expanded.s = NULL;
  decompress(c, &expanded);
  free_cdata(c);
  *c = expanded;
}

void cdata_compress(cdata_t *c) {
  if (c->compressed) wzfatal("Already compressed");
  switch(c->fmt) {
  case '0': { break; }
  case '1': { fmt1_compress(c); break; }
  case '2': { fmt2_compress(c); break; }
  case '3': { fmt3_compress(c); break; }
  case '4': { fmt4_compress(c); break; }
  case '5': { fmt5_compress(c); break; }
  case '6': { fmt6_compress(c); break; }
  default: wzfatal("Unrecognized format: %c.\n", c->fmt);
  }
}
