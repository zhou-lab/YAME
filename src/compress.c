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

void fmt1_compress(cdata_t *c);
void fmt2_compress(cdata_t *c);
void fmt3_compress(cdata_t *c);
void fmt4_compress(cdata_t *c);
void fmt5_compress(cdata_t *c);
void fmt6_compress(cdata_t *c);

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
