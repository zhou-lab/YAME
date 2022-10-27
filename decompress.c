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
