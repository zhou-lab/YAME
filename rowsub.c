#include <sys/stat.h>
#include <sys/types.h>
#include "cfile.h"

typedef struct config_t {
  char *fname_rindex;
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
  fprintf(stderr, "    -l [PATH] rows in 1-based indices. If given, other index formats are ignored.\n");
  fprintf(stderr, "    -L [PATH] rows in the format of chrm_beg1. This option requires -R.\n");
  fprintf(stderr, "    -R [PATH] co-subset row coordinates.\n");
  fprintf(stderr, "              The subset row coordinate is appended in output as the first file.\n");
  fprintf(stderr, "    -m [PATH] rows specified as a mask file (format 1 or 2).\n");
  fprintf(stderr, "    -b        rows as a range (begin-end format): begin index, 0-base.\n");
  fprintf(stderr, "    -e        rows as a range (begin-end format): end index 1-base.\n");
  fprintf(stderr, "    -i        rows as a range (index-block-size format): index (0-base).\n");
  fprintf(stderr, "    -s        rows as a range (index-block-size format): size (N=%"PRIu64").\n", config->isize);
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

static int split_string_and_number(const char* input, char** out_string, uint64_t* out_number) {
  // Find the position of the underscore ('_')
  const char *underscore = strchr(input, '_');
  if (!underscore) return -1; // underscore not found

  // Allocate memory for the string part (including the null-terminator)
  *out_string = malloc(underscore - input + 1);
  if (!*out_string) return -2; // memory allocation failed

  // Copy the string part to the output
  strncpy(*out_string, input, underscore - input);
  (*out_string)[underscore - input] = '\0'; // Null-terminate

  // Convert the number part to uint64_t
  char *endptr;
  *out_number = strtoull(underscore + 1, &endptr, 10);

  // Check if conversion was successful
  if (*endptr != '\0') { free(*out_string); return -3; }
  return 0; // success
}

static int64_t *load_row_indices_by_names(char *fname_rnindex, cdata_t *cr, int64_t *n_indices) {
  int64_t *indices = NULL;
  /* snames_t snames = {0}; */
  if (fname_rnindex == NULL) { *n_indices = 0; return indices; }
  gzFile fp;
  if (strcmp(fname_rnindex, "-") == 0) {
    fp = gzdopen(fileno(stdin), "r");
  } else {
    fp = gzopen(fname_rnindex, "r");
    if (!fp) {
      fprintf(stderr, "[%s:%d] Fatal, cannot open file: %s\n",
              __func__, __LINE__, fname_rnindex);
      fflush(stderr);
      exit(1);
    }
  }
  
  if (fp == NULL) return indices;
  
  row_finder_t fdr = init_finder(cr);
  *n_indices = 0;
  char *line = NULL;
  while (gzFile_read_line(fp, &line) > 0) {
    char *chrm = NULL;
    uint64_t beg1 = 0;
    if (split_string_and_number(line, &chrm, &beg1) < 0) {
      fprintf(stderr, "Failed to extract coordinate: %s\n", line);
      fflush(stderr);
      exit(1);
    }
    indices = realloc(indices, ((*n_indices)+1)*sizeof(int64_t));
    indices[(*n_indices)] = row_finder_search(chrm, beg1, &fdr, cr);
    if (!indices[(*n_indices)]) {
      fprintf(stderr, "[%s:%d] Cannot find coordinate: %s\n", __func__, __LINE__, line);
      fflush(stderr);
      exit(1);
    }
    (*n_indices)++;
    free(chrm);
  }
  free_row_finder(&fdr);
  free(line);
  gzclose(fp);
  return indices;
}

static cdata_t sliceToIndices(cdata_t *c, int64_t *row_indices, int64_t n) {
  if (c->compressed) {
    fprintf(stderr, "[%s:%d] Slicing compressed data.\n", __func__, __LINE__);
    fflush(stderr);
    exit(1);
  }
  cdata_t c2 = {0};
  c2.unit = c->unit;
  c2.s = realloc(c2.s, n*c2.unit);
  c2.fmt = c->fmt;
  int64_t i;
  for (i=0; i<n; ++i) {
    memcpy(c2.s+c2.unit*i, c->s+c->unit*(row_indices[i]-1), c->unit);
  }
  c2.n = n;
  c2.compressed = 0;
  return c2;
}

// beg and end are both 0-based
static cdata_t sliceToBlock(cdata_t *c, uint64_t beg, uint64_t end) {
  if (c->compressed) {
    fprintf(stderr, "[%s:%d] Cannot slice compressed data.\n", __func__, __LINE__);
    fflush(stderr);
    exit(1);
  }
  if (end > c->n-1) end = c->n-1; // 0-base
  if (beg > c->n-1) {
    fprintf(stderr, "[%s:%d] Begin (%"PRIu64") is bigger than the data vector size (%"PRIu64").\n", __func__, __LINE__, beg, c->n);
    fflush(stderr);
    exit(1);
  }

  cdata_t c2 = {0};
  c2.unit = c->unit;
  c2.s = realloc(c2.s, (end-beg+1)*c2.unit);
  c2.fmt = c->fmt;
  memcpy(c2.s, c->s+c->unit*beg, c->unit*(end-beg+1));
  c2.n = end-beg+1;
  c2.compressed = 0;
  return c2;
}

static cdata_t sliceToMask(cdata_t *c, cdata_t *c_mask) {
  if (c->compressed) {
    fprintf(stderr, "[%s:%d] Slicing compressed data.\n", __func__, __LINE__);
    fflush(stderr);
    exit(1);
  }
  if (c->n != c_mask->n) {
    fprintf(stderr, "[%s:%d] Mask (N=%"PRIu64") and data (N=%"PRIu64") are of different lengths.\n",
            __func__, __LINE__, c_mask->n, c->n);
    fflush(stderr);
    exit(1);
  }

  uint64_t n = 0;
  for (uint64_t i=0; i<c->n; ++i)
    if (c_mask->s[i>>3]&(1<<(i&0x7))) n++;

  cdata_t c2 = {0};
  c2.unit = c->unit;
  c2.s = realloc(c2.s, n*c2.unit);
  c2.fmt = c->fmt;
  for (uint64_t i=0, k=0; i<c->n; ++i) {
    if (c_mask->s[i>>3]&(1<<(i&0x7)))
      memcpy(c2.s+(k++)*c->unit, c->s+i*c->unit, c->unit);
  }
  c2.n = n;
  c2.compressed = 0;
  return c2;

}

int main_rowsub(int argc, char *argv[]) {

  config_t config = {
    .fname_rindex = NULL,
    .index = -1, .isize = 1000000, .beg = 0, .end = 1};
  int c; char *fname_row = NULL; char *fname_mask = NULL;
  char *fname_rnindex = NULL;
  while ((c = getopt(argc, argv, "R:m:l:L:b:e:i:s:vh"))>=0) {
    switch (c) {
    case 'R': fname_row = strdup(optarg); break;
    case 'm': fname_mask = strdup(optarg); break;
    case 'l': config.fname_rindex = strdup(optarg); break;
    case 'L': fname_rnindex = strdup(optarg); break;
    case 'b': config.beg = atoi(optarg); break;   /* assume 0-index */
    case 'e': config.end = atoi(optarg)-1; break; /* convert to 0-index */
    case 'i': config.index = atoi(optarg); break; /* 0-index */
    case 's': config.isize = atoi(optarg); break;
    case 'h': return usage(&config); break;
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
  if (config.fname_rindex) {
    row_indices = load_row_indices(config.fname_rindex, &n_indices);
  }

  cdata_t c_mask = {0};
  if (fname_mask) {
    cfile_t cf_mask = open_cfile(fname_mask);
    c_mask = read_cdata1(&cf_mask);
    if (c_mask.fmt >= '2') {
      fprintf(stderr, "[%s:%d] Mask is not binary.\n", __func__, __LINE__);
      fflush(stderr);
      exit(1);
    }
    convertToFmt0(&c_mask);
    bgzf_close(cf_mask.fh);
  }

  cfile_t cf = open_cfile(fname);
  BGZF *fp_out = bgzf_dopen(fileno(stdout), "w");
  if (fp_out == NULL) {
    fprintf(stderr, "[%s:%d] Cannot open output stream.\n", __func__, __LINE__);
    fflush(stderr);
    exit(1);
  }

  if (fname_rnindex && !fname_row) {
    fprintf(stderr, "[%s:%d] Missing -R for BED coordinates.\n", __func__, __LINE__);
    fflush(stderr);
    exit(1);
  }

  if (fname_row) {
    cfile_t cf_row = open_cfile(fname_row);
    cdata_t cr = read_cdata1(&cf_row);
    if (!row_indices && fname_rnindex) {
      row_indices = load_row_indices_by_names(fname_rnindex, &cr, &n_indices);
    }
    
    cdata_t cr2;
    if (row_indices) cr2 = fmt7_sliceToIndices(&cr, row_indices, n_indices);
    else if (c_mask.n) cr2 = fmt7_sliceToMask(&cr, &c_mask);
    else cr2 = fmt7_sliceToBlock(&cr, config.beg, config.end);
    cdata_write1(fp_out, &cr2);
    free_cdata(&cr);
    free_cdata(&cr2);
    bgzf_close(cf_row.fh);
  }
  
  while (1) {
    cdata_t c = read_cdata1(&cf);
    if (c.n == 0) break;
    if (c.fmt == '7') {
      cdata_t c2;
      if (row_indices) c2 = fmt7_sliceToIndices(&c, row_indices, n_indices);
      else if (c_mask.n) c2 = fmt7_sliceToMask(&c, &c_mask);
      else c2 = fmt7_sliceToBlock(&c, config.beg, config.end);
      cdata_write1(fp_out, &c2);
      free_cdata(&c2);
    } else {
      cdata_t c2 = {0};
      decompress(&c, &c2);
      cdata_t c3 = {0};
      if (row_indices) c3 = sliceToIndices(&c2, row_indices, n_indices);
      else if (c_mask.n) c3 = sliceToMask(&c2, &c_mask);
      else c3 = sliceToBlock(&c2, config.beg, config.end);
      cdata_compress(&c3);
      cdata_write1(fp_out, &c3);
      free(c3.s); free(c2.s);
    }
    free_cdata(&c);
  }
  bgzf_close(cf.fh);
  bgzf_close(fp_out);

  if (n_indices) free(row_indices);
  if (fname_row) free(fname_row);
  if (fname_rnindex) free(fname_rnindex);
  
  return 0;
}
