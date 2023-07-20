#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdio.h>
#include "wzmisc.h"
#include "wzbed.h"
#include "cgfile.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: kycg overlap [options] <feature.cg> <query.cg>\n");
  fprintf(stderr, "Query can be a multi-sample set. only the first feature will be used.\n");
  fprintf(stderr, "If the input is format 3: test overlap of M+U.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -u        Universe .cg file (optional).\n");
  fprintf(stderr, "    -h        This help.\n");
  fprintf(stderr, "Returns four numbers, nu, nf, nq, nfq;\n");
  fprintf(stderr, "\n");

  return 1;
}

static void bit_mask(uint8_t *s, uint8_t *mask, size_t n) {
  size_t i;
  for (i=0; i<(n>>3)+1; ++i) s[i] &= mask[i];
}

static size_t bit_count(cgdata_t cg) {

  /* create a look-up table */
  int byte2cnt[256]; int p;
  for (p=0; p<256; ++p) {
    unsigned char q = p; int ii, cnt = 0;
    for (ii=0; ii<8; ++ii) { if (q&1) cnt++; q>>=1; }
    byte2cnt[p] = cnt;
  }
  
  size_t i,k,m = 0;
  for (i=0; i<(cg.n>>3); ++i) m += byte2cnt[cg.s[i]];
  for (k=0; k<(cg.n&0x7); ++k) m += (cg.s[i]>>k) & 0x1;
  return m;
}

/**
 * Reads data from a given cgfile, converts it into format 0 if needed, and returns it.
 * Format 0 refers to a bit-packed format where 1 byte represents an 8-bit binary vector.
 *
 * @param cgf: A pointer to a cgfile_t structure from which data is to be read.
 * @return: A cgdata_t structure which contains the read and converted data.
 * 
 * The function performs the following operations:
 * 
 * 1. It begins by reading data from the given cgfile into a cgdata_t structure.
 * 2. Input format is '0': no further operation is performed.
 * 3. Input format is '1': if the value is 0, return 0 else 1.
 * 4. Input format is '3': if the M+U is 0, return 0 else 1
 * 5. Other input formats are not allowed.
 */
static void convertToFmt0(cgdata_t *cg) {
  cgdata_t cg_out = {0};
  switch (cg->fmt) {
  case '0': return;
  case '1': {
    cg_out.fmt = '0';
    cg_out.compressed = 1;
    cg_out.n=0;
    uint64_t i;
    for (i=0; i<cg->n/3; ++i) {
      cg_out.n += *((uint16_t*) (cg->s+i*3+1));
    }
    cg_out.s = calloc((cg_out.n>>3)+1, 1);
    size_t sum; uint16_t l=0;
    for (i=0, sum=0; i<cg->n/3; ++i, sum+=l) {
      l = *((uint16_t*) (cg->s+i*3+1));
      if (cg->s[i*3] > '0') {
        for(size_t j=sum; j<sum+l; ++j) {
          cg_out.s[j>>3] |= (1<<(j&0x7));
        }
      }
    }
    break;
  }
  case '3': {
    cgdata_t expanded = {0};
    fmt3_decompress(cg, &expanded);

    cg_out.fmt = '0';
    cg_out.compressed = 1;
    cg_out.n = expanded.n;
    cg_out.s = calloc((cg_out.n>>3)+1,1);
    uint64_t *s = (uint64_t*) expanded.s;
    for (uint64_t i=0; i<expanded.n; ++i) {
      if (s[i] > 0) { /* equivalent to: 1 if M+U > 0 else 0 */
        cg_out.s[i>>3] |= (1<<(i&0x7));
      }
    }
    free(expanded.s);
    break;
  }
  default: wzfatal("Format %c unsupported.\n", cg->fmt);
  }
  free(cg->s);
  *cg = cg_out;
}

/* The design, first 10 bytes are uint64_t (length) + uint16_t (0=vec; 1=rle) */
int main_overlap(int argc, char *argv[]) {
  int c; char *upath = NULL;
  while ((c = getopt(argc, argv, "u:h"))>=0) {
    switch (c) {
    case 'u': upath = strdup(optarg); break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 2 > argc) { 
    usage(); 
    wzfatal("Please supply input and output file.\n"); 
  }


  /* import universe if provided */
  cgdata_t cg_uni = {0};
  size_t m_uni = 0;
  if (upath != NULL) {
    cgfile_t cgf_uni = open_cgfile(upath);
    cg_uni = read_cg(&cgf_uni);
    m_uni = bit_count(cg_uni);
    free(upath);
    bgzf_close(cgf_uni.fh);
  }
  
  cgfile_t cgf_fea = open_cgfile(argv[optind++]);
  for (uint64_t kf=0;;++kf) {
    cgdata_t cg_fea = read_cg(&cgf_fea);
    if (cg_fea.n == 0) break;
    convertToFmt0(&cg_fea);

    if (cg_uni.n > 0) {
      assert(cg_uni.n == cg_fea.n);
      bit_mask(cg_fea.s, cg_uni.s, cg_uni.n);
    } else { m_uni = cg_fea.n; }
    size_t nf = bit_count(cg_fea);
    
    cgfile_t cgf_qry = open_cgfile(argv[optind]);
    for (uint64_t kq=0;;++kq) {
      cgdata_t cg_qry = read_cg(&cgf_qry);
      if (cg_qry.n == 0) break;
      convertToFmt0(&cg_qry);

      if (cg_fea.n != cg_qry.n) {
        wzfatal("Query (%"PRId64") and feature (%"PRId64") has different lengths.",
                cg_qry.n, cg_fea.n);
      }

      if (cg_uni.n > 0) {
        assert(cg_uni.n == cg_qry.n);
        bit_mask(cg_qry.s, cg_uni.s, cg_uni.n);
      }

      size_t nq = bit_count(cg_qry);
      bit_mask(cg_qry.s, cg_fea.s, cg_fea.n);
      size_t nfq = bit_count(cg_qry);
      fprintf(stdout, "%"PRIu64"\t%"PRIu64"\t%zu\t%zu\t%zu\t%zu\n",
              kf+1, kq+1, m_uni, nf, nq, nfq);
      
      free(cg_qry.s);
    }
    free(cg_fea.s);
    bgzf_close(cgf_qry.fh);
  }
  if (cg_uni.n > 0) free(cg_uni.s);
  bgzf_close(cgf_fea.fh);
  
  return 0;
}
