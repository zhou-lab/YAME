#include <sys/stat.h>
#include <sys/types.h>
#include "cfile.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: yame rowsub [options] <in.cx> <row.index.list>\n");
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

static void sliceToIndices(cdata_t *c, int64_t *row_indices, int64_t n, cdata_t *c2) {

  assert(!c->compressed);
  c2->unit = c->unit;
  c2->s = realloc(c2->s, n*c2->unit);
  c2->fmt = c->fmt;
  int64_t i;
  for (i=0; i<n; ++i) {
    memcpy(c2->s+c2->unit*i, c->s+c->unit*(row_indices[i]-1), c->unit);
  }
  c2->n = n;
  c2->compressed = 0;
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

  cfile_t cf = open_cfile(fname);
  BGZF *fp_out = bgzf_dopen(fileno(stdout), "w");
  assert(fp_out != NULL);
  while (1) {
    cdata_t c = read_cdata1(&cf);
    if (c.n == 0) break;
    
    cdata_t c2 = {0};
    decompress(&c, &c2);
    cdata_t c3 = {0};
    c3.s = NULL;
    sliceToIndices(&c2, row_indices, n, &c3);
    cdata_compress(&c3);
    cdata_write1(fp_out, &c3);
    free(c3.s); free(c2.s); free(c.s);
  }
  bgzf_close(cf.fh);
  bgzf_close(fp_out);
  
  return 0;
}
