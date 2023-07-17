#include <sys/stat.h>
#include <sys/types.h>
#include "cgfile.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: kycg rowsub [options] <in.cg> <row.index.list>\n");
  fprintf(stderr, "This function outputs to stdout.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -v        verbose\n");
  fprintf(stderr, "    -h        This help\n");
  fprintf(stderr, "\n");

  return 1;
}

static int64_t *load_row_indices(char *fname, int64_t *n) {

  int64_t *indices = NULL;
  /* snames_t snames = {0}; */
  if (fname == NULL) { *n = 0; return indices; }
  gzFile fp;
  if (strcmp(fname, "-") == 0) {
    fp = gzdopen(fileno(stdin), "r");
  } else {
    fp = gzopen(fname, "r");
    if (!fp) {
      fprintf(stderr, "[%s:%d] Fatal, cannot open file: %s\n",
              __func__, __LINE__, fname);
      fflush(stderr);
      exit(1);
    }
  }
  
  if (fp == NULL) return indices;

  char *line = NULL; *n=0;
  while (gzFile_read_line(fp, &line) > 0) {
    char *field = NULL;
    if (line_get_field(line, 0, "\t", &field)) {
      indices = realloc(indices, sizeof(int64_t)*((*n)+1));
      if (indices == NULL) {
        fprintf(stderr, "Failed to allocate memory\n");
        fflush(stderr);
        exit(1);
      }
      indices[(*n)++] = strtoll(field, NULL, 10);
      free(field);
    }
  }
  free(line);
  gzclose(fp);
  return indices;
}

static void sliceToIndices(cgdata_t *cg, int64_t *row_indices, int64_t n, cgdata_t *cg2) {

  assert(!cg->compressed);
  cg2->s = realloc(cg2->s, n*cgdata_unit_size(cg));
  cg2->fmt = cg->fmt;
  int64_t i;
  for (i=0; i<n; ++i) {
    uint64_t *s = (uint64_t*) cg->s;
    uint64_t *s2 = (uint64_t*) cg2->s;
    s2[i] = s[row_indices[i]-1]; // input is 1-based
  }
  cg2->n = n;
  cg2->compressed = 0;
}

int main_rowsub(int argc, char *argv[]) {

  int c;
  while ((c = getopt(argc, argv, "vh"))>=0) {
    switch (c) {
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 2 > argc) {
    usage(); 
    wzfatal("Please supply input files.\n");
  }

  char *fname = argv[optind];
  int64_t n=0;
  int64_t *row_indices = load_row_indices(argv[optind+1], &n);

  cgfile_t cgf = open_cgfile(fname);
  BGZF *fp_out = bgzf_dopen(fileno(stdout), "w");
  assert(fp_out != NULL);
  for (uint64_t k=0; ; ++k) {
    cgdata_t cg = read_cg(&cgf);
    if (cg.n == 0) break;
    
    cgdata_t cg2 = {0};
    decompress(&cg, &cg2);
    cgdata_t cg3 = {0};
    cg3.s = NULL;
    sliceToIndices(&cg2, row_indices, n, &cg3);
    recompress(&cg3);
    cgdata_write1(fp_out, &cg3);
    free(cg3.s); free(cg2.s); free(cg.s);
  }
  bgzf_close(cgf.fh);
  bgzf_close(fp_out);
  
  return 0;
}
