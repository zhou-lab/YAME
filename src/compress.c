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
