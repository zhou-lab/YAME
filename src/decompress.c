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

