#include <sys/stat.h>
#include <sys/types.h>
#include "cfile.h"

typedef struct config_t {
  char *fname_index;
  uint64_t beg;
  uint64_t end;
  int64_t index;
  int64_t isize;
} config_t;

static int usage(config_t *config) {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: yame rowsub [options] <in.cx>\n");
  fprintf(stderr, "This function outputs to stdout.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -v        verbose\n");
  fprintf(stderr, "    -l        File name for row index list. If this is given, other index formats are ignored.\n");
  fprintf(stderr, "    -b        begin-end format: begin index, 0-base.\n");
  fprintf(stderr, "    -e        begin-end format: end index 1-base.\n");
  fprintf(stderr, "    -i        index-block-size format: index (0-base).\n");
  fprintf(stderr, "    -s        index-block-size format: size (N=%"PRIu64").\n", config->isize);
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

static void sliceToBlock(cdata_t *c, uint64_t beg, uint64_t end, cdata_t *c2) {
  if (c->compressed) {
    fprintf(stderr, "[%s:%d] Input is compressed.\n", __func__, __LINE__);
    fflush(stderr);
    exit(1);
  }
  if (end > c->n) end = c->n-1; // 0-base
  if (beg > c->n) {
    fprintf(stderr, "[%s:%d] Begin (%"PRIu64") is bigger than the data vector size (%"PRIu64").\n", __func__, __LINE__, beg, c->n);
    fflush(stderr);
    exit(1);
  }
  
  uint64_t n = end - beg + 1;
  c2->unit = c->unit;
  c2->s = realloc(c2->s, n*c2->unit);
  c2->fmt = c->fmt;
  memcpy(c2->s, c->s+c->unit*beg, c->unit*n);
  c2->n = n;
  c2->compressed = 0;
}

int main_rowsub(int argc, char *argv[]) {

  config_t config = {
    .fname_index = NULL,
    .index = -1, .isize = 1000000, .beg = 0, .end = 1};
  int c;
  while ((c = getopt(argc, argv, "b:e:i:s:vh"))>=0) {
    switch (c) {
    case 'h': return usage(&config); break;
    case 'l': config.fname_index = optarg; break;
    case 'b': config.beg = atoi(optarg); break;   /* assume 0-index */
    case 'e': config.end = atoi(optarg)-1; break; /* convert to 0-index */
    case 'i': config.index = atoi(optarg); break; /* 0-index */
    case 's': config.isize = atoi(optarg); break;
    default: usage(&config); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (config.index >= 0) {
    config.beg = config.index*config.isize;
    config.end = (config.index+1)*config.isize-1;
  }

  if (optind + 1 > argc) {
    usage(&config); 
    wzfatal("Please supply input files.\n");
  }
  char *fname = argv[optind];

  int64_t *row_indices = NULL;
  int64_t n_indices=0;
  if (config.fname_index) {
    row_indices = load_row_indices(config.fname_index, &n_indices);
  }

  cfile_t cf = open_cfile(fname);
  BGZF *fp_out = bgzf_dopen(fileno(stdout), "w");
  if (fp_out == NULL) {
    fprintf(stderr, "[%s:%d] Cannot open output stream.\n", __func__, __LINE__);
    fflush(stderr);
    exit(1);
  }
  while (1) {
    cdata_t c = read_cdata1(&cf);
    if (c.n == 0) break;
    
    cdata_t c2 = {0};
    decompress(&c, &c2);
    cdata_t c3 = {0};
    c3.s = NULL;
    if (row_indices) {
      sliceToIndices(&c2, row_indices, n_indices, &c3);
    } else {
      sliceToBlock(&c2, config.beg, config.end, &c3);
    }
    cdata_compress(&c3);
    cdata_write1(fp_out, &c3);
    free(c3.s); free(c2.s); free(c.s);
  }
  bgzf_close(cf.fh);
  bgzf_close(fp_out);

  if (n_indices) free(row_indices);
  
  return 0;
}
