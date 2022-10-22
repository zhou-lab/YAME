#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdio.h>
#include "wzmisc.h"
#include "wzbed.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: kycg pack [options] <in.bed> <out.cg>\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -h        This help\n");
  fprintf(stderr, "\n");

  return 1;
}

/* The header design, 17 bytes
   uint64_t x 1: num columns
   uint8_t x 1: format (0=vec; 1=rle)
   uint64_t x 1: length (n_cgs or n_bytes for rle)
*/

int main_pack(int argc, char *argv[]) {

  int c;
  while ((c = getopt(argc, argv, "h"))>=0) {
    switch (c) {
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 2 > argc) { 
    usage(); 
    wzfatal("Please supply input and output file.\n"); 
  }
  
  gzFile fh = wzopen(argv[optind++]);
  char *line = NULL;
  int64_t n=0, m=1<<22; int64_t n1=0;
  char *s=calloc(m, 1);
  while (gzFile_read_line(fh, &line) > 0) {
    if (line[0] != '0') {
      s[n>>3] |= (1<<(n&0x7));
      n1++;
    }
    n++;
    if (n>m-2) { m<<=1; s=realloc(s,m); }
  }
  free(line);
  wzclose(fh);
  fprintf(stderr, "Vector of length %"PRId64" loaded\n", n);

  /* see if rle gives shorter storage */
  int64_t n_rle=0;
  uint8_t *s_rle = NULL;
  int64_t i=0; uint16_t l=0; uint8_t u0=0;
  uint64_t sum_l=0;
  for (i=0, l=0; i<n; ++i) {
    /* unsigned char u = s[i>>3] & (1<<(n&0x7)); */
    uint8_t u = (s[i>>3]>>(i&0x7))&0x1;
    /* either not the same as before or reach block size max */
    if ((l != 0 && u != u0) || l+2 >= 1<<15) {
      s_rle = realloc(s_rle, n_rle+3);
      s_rle[n_rle] = u0;
      *((uint16_t*) (s_rle+n_rle+1)) = l;
      sum_l += l;
      n_rle += 3;
      l = 1;
    } else {
      ++l;
    }
    u0 = u;
  }
  /* the last rle */
  s_rle = realloc(s_rle, n_rle+3);
  s_rle[n_rle] = u0;
  /* s_rle[n_rle+1] = l; */
  *((uint16_t*) (s_rle+n_rle+1)) = l;
  sum_l += l;
  n_rle += 3;
  fprintf(stderr, "RLE sum: %"PRIu64"\n", sum_l);
  fprintf(stderr, "N_vec: %"PRId64"; N_rle: %"PRId64"\n", n>>3, n_rle);

  FILE *out = fopen(argv[optind], "wb");
  uint8_t fmt; uint64_t dat;
  if (n>>3 > n_rle) {           /* rle */
    dat = 1; fwrite(&dat, sizeof(uint64_t), 1, out); /* number of columns */
    fmt = 1; fwrite(&fmt, sizeof(uint8_t), 1, out);
    fwrite(&n_rle, sizeof(int64_t), 1, out);
    fwrite(s_rle, 1, n_rle, out);
    fputs("Stored as Run-length-encoded vector\n", stderr);
  } else {                      /* plain vector */
    dat = 1; fwrite(&dat, sizeof(uint64_t), 1, out); /* number of columns */
    fmt = 0; fwrite(&fmt, 1, 1, out);
    fwrite(&n, sizeof(int64_t), 1, out);
    fwrite(s, 1, (n>>3)+1, out);
    fputs("Stored as bit-vector.\n", stderr);
  }
  fclose(out);

  free(s); free(s_rle);
  return 0;
}
