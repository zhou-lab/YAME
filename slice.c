#include "kycg.h"

void slice(cgdata_t *cg, uint64_t beg, uint64_t end, cgdata_t *cg_sliced) {
  if (end > cg->n) end = cg->n;
  if (end < beg) wzfatal("Slicing negative span.");

  cg_sliced->s = realloc(cg_sliced->s, (end-beg+1)*sizeof(uint8_t));
  memcpy(cg_sliced->s, cg->s, (end-beg+1)*sizeof(uint8_t));
  cg_sliced->n = end - beg + 1;
  cg_sliced->compressed = 0;
  cg_sliced->fmt = 5;
}
