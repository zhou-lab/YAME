#include "cgfile.h"

void decompress(cgdata_t *cg, cgdata_t *expanded) {
  switch (cg->fmt) {
  case '0': { fmt0_decompress(cg, expanded); break; }
  case '1': { fmt1_decompress(cg, expanded); break; }
  case '2': { fmt2_decompress(cg, expanded); break; }
  case '3': { fmt3_decompress(cg, expanded); break; }
  case '4': { fmt4_decompress(cg, expanded); break; }
  case '5': { fmt5_decompress(cg, expanded); break; }
  case '6': { fmt6_decompress(cg, expanded); break; }
  default: wzfatal("Unrecognized format: %c.\n", cg->fmt);
  }
  /* shouldn't reach here */
}

void decompress2(cgdata_t *cg) {
  cgdata_t expanded = *cg;
  expanded.s = NULL;
  decompress(cg, &expanded);
  *cg = expanded;
}

void cdata_compress(cgdata_t *cg) {
  if (cg->compressed) wzfatal("Already compressed");
  switch(cg->fmt) {
  case '0': { break; }
  case '1': { fmt1_compress(cg); break; }
  case '2': { fmt2_compress(cg); break; }
  case '3': { fmt3_compress(cg); break; }
  case '4': { fmt4_compress(cg); break; }
  case '5': { fmt5_compress(cg); break; }
  case '6': { fmt6_compress(cg); break; }
  default: wzfatal("Unrecognized format: %c.\n", cg->fmt);
  }
}
