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

cdata_t fmt0_decompress(const cdata_t c);
cdata_t fmt1_decompress(const cdata_t c);
cdata_t fmt2_decompress(const cdata_t c);
cdata_t fmt3_decompress(const cdata_t c);
cdata_t fmt4_decompress(const cdata_t c);
cdata_t fmt5_decompress(const cdata_t c);
cdata_t fmt6_decompress(const cdata_t c);
cdata_t fmt7_decompress(const cdata_t c);

cdata_t decompress(cdata_t c) {
  switch (c.fmt) {
  case '0': { return fmt0_decompress(c); }
  case '1': { return fmt1_decompress(c); }
  case '2': { return fmt2_decompress(c); }
  case '3': { return fmt3_decompress(c); }
  case '4': { return fmt4_decompress(c); }
  case '5': { return fmt5_decompress(c); }
  case '6': { return fmt6_decompress(c); }
  case '7': { return fmt7_decompress(c); }
  default: wzfatal("Unsupported format for inflation: %c.\n", c.fmt);
  }
  return c; /* shouldn't reach here */
}

void decompress_in_situ(cdata_t *c) {
  if (!c->compressed) {
    fprintf(stderr, "[%s:%d] Already decompressed.\n", __func__, __LINE__);
    fflush(stderr);
    exit(1);
  }
  cdata_t expanded = decompress(*c);
  free_cdata(c);
  *c = expanded;
}

