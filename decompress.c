#include "kycg.h"

cgdata_t decompress(cgdata_t *cg) {
  switch (cg->fmt) {
  case '0': { return fmt0_decompress(cg); break; }
  case '1': { return fmt1_decompress(cg); break; }
  case '3': { return fmt3_decompress(cg); break; }
  case '4': { return fmt4_decompress(cg); break; }
  case '5': { return fmt5_decompress(cg); break; }
  default: wzfatal("Unrecognized format: %c.\n", cg->fmt);
  }
  /* shouldn't reach here */
  return fmt0_decompress(cg);
}

void compress(cgdata_t *cg) {
  if (cg->compressed) wzfatal("Already compressed");
  switch(cg->fmt) {
  case '0': { fmt0_compress(cg); break; }
  case '1': { fmt1_compress(cg); break; }
  case '3': { fmt3_compress(cg); break; }
  case '4': { fmt4_compress(cg); break; }
  case '5': { fmt5_compress(cg); break; }
  default: wzfatal("Unrecognized format: %c.\n", cg->fmt);
  }
}
