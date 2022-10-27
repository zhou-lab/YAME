#include "kycg.h"

cgdata_t slice(cgdata_t *cg, uint64_t beg, uint64_t end) {
  if (end > cg->n) end = cg->n;
  if (end < beg) wzfatal("Slicing negative span.");

  uint8_t *s = (uint8_t*) cg->s;
  cgdata_t cg2 = {0};
  cg2.s = malloc((end-beg+1)*sizeof(uint8_t));
  memcpy(cg2.s, cg->s, (end-beg+1)*sizeof(uint8_t));
  cg2.n = end - beg + 1;
  cg2.compressed = 0;
  cg2.fmt = 5;
  return cg2;
}
