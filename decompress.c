#include "cfile.h"

void decompress(cdata_t *c, cdata_t *expanded) {
  switch (c->fmt) {
  case '0': { fmt0_decompress(c, expanded); break; }
  case '1': { fmt1_decompress(c, expanded); break; }
  case '2': { fmt2_decompress(c, expanded); break; }
  case '3': { fmt3_decompress(c, expanded); break; }
  case '4': { fmt4_decompress(c, expanded); break; }
  case '5': { fmt5_decompress(c, expanded); break; }
  case '6': { fmt6_decompress(c, expanded); break; }
  default: wzfatal("Unrecognized format: %c.\n", c->fmt);
  }
  /* shouldn't reach here */
}

void decompress2(cdata_t *c) {
  cdata_t expanded = *c;
  expanded.s = NULL;
  decompress(c, &expanded);
  *c = expanded;
}

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
