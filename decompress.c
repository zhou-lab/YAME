#include "kycg.h"

void decompress(cgdata_t *cg, cgdata_t *expanded) {
  switch (cg->fmt) {
  case '0': { fmt0_decompress(cg, expanded); break; }
  case '1': { fmt1_decompress(cg, expanded); break; }
  case '3': { fmt3_decompress(cg, expanded); break; }
  case '4': { fmt4_decompress(cg, expanded); break; }
  case '5': { fmt5_decompress(cg, expanded); break; }
  default: wzfatal("Unrecognized format: %c.\n", cg->fmt);
  }
  /* shouldn't reach here */
}

void recompress(cgdata_t *cg) {
  if (cg->compressed) wzfatal("Already compressed");
  switch(cg->fmt) {
  case '0': { break; }
  case '1': { fmt1_compress(cg); break; }
  case '3': { fmt3_compress(cg); break; }
  case '4': { fmt4_compress(cg); break; }
  case '5': { fmt5_compress(cg); break; }
  default: wzfatal("Unrecognized format: %c.\n", cg->fmt);
  }
}
