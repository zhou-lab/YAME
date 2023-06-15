#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "kycg.h"

index_t* loadIndex(char* fname_cg) {
  char* fname_index = malloc(strlen(fname_cg) + strlen(".idx") + 1);
  if (fname_index == NULL) {
    printf("Failed to allocate memory for index file name\n");
    return NULL;
  }
  strcpy(fname_index, fname_cg);
  strcat(fname_index, ".idx");

  gzFile file = wzopen(fname_index, 0);
  if (file == NULL) {
    printf("Failed to open index file\n");
    free(fname_index);
    return NULL;
  }

  index_t* idx = kh_init(index);
  if (idx == NULL) {
    printf("Failed to create hash table\n");
    wzclose(file);
    free(fname_index);
    return NULL;
  }

  char* line = NULL;
  while (gzFile_read_line(file, &line) > 0) {
    char* sample_name;
    if (line_get_field(line, 0, "\t", &sample_name)) {
      char* index_str;
      if (line_get_field(line, 1, "\t", &index_str)) {
        int ret;
        khiter_t k = kh_put(index, idx, sample_name, &ret);
        if (ret == -1) {
          printf("Failed to insert value into hash table\n");
          wzclose(file);
          free(fname_index);
          free(line);
          kh_destroy(index, idx);
          free(idx);
          return NULL;
        }
        kh_value(idx, k) = strtoll(index_str, NULL, 10);
        free(index_str);
      }
      /* free(sample_name); */
    }
  }
  free(line);
  wzclose(file);
  free(fname_index);

  return idx;
}
  
static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: kycg index [options] <in.cg>\n");
  fprintf(stderr, "The index file name default to <in.cg>.idx\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -s        tab-delimited sample name list (use first column) \n");
  fprintf(stderr, "    -c        output index to console\n");
  fprintf(stderr, "    -h        This help\n");
  fprintf(stderr, "\n");

  return 1;
}

typedef struct {
  int n;
  char **array;
} snames_t;

static snames_t* loadSampleNames(char* fname_snames) {
  if (fname_snames == NULL) return NULL;
  gzFile file = wzopen(fname_snames, 0);
  if (file == NULL) return NULL;

  snames_t *snames = malloc(sizeof(snames_t));
  if (snames == NULL) {
    printf("Failed to allocate memory for sample name vector\n");
    return NULL;
  }
  snames->array = NULL;
  snames->n = 0;

  char *line = NULL;
  while (gzFile_read_line(file, &line) > 0) {
    char *sname;
    if (line_get_field(line, 0, DELIMITER, &sname)) {
      snames->array = realloc(snames->array, sizeof(*(snames->array)) * (snames->n + 1));
      if (snames->array == NULL) {
        printf("Failed to allocate memory\n");
        return NULL;
      }
      snames->array[snames->n] = sname;
      snames->n++;
    }
  }
  free(line);
  line = NULL;
  wzclose(file);

  return snames;
}

static void cleanUpSnames(snames_t *sname) {
  for (int i = 0; i < sname->n; i++) {
    free(sname->array[i]);
  }
  free(sname->array);
  free(sname);
}

int main_index(int argc, char *argv[]) {

  int c, console = 0;
  char *fname_snames = NULL;
  while ((c = getopt(argc, argv, "cs:h"))>=0) {
    switch (c) {
    case 'c': console = 1; break;
    case 's': fname_snames = strdup(optarg); break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 1 > argc) { 
    usage(); 
    wzfatal("Please supply input file.\n"); 
  }

  snames_t *snames = loadSampleNames(fname_snames);
  char *fname_index = NULL;
  fname_index = malloc(strlen(argv[optind]) + strlen(".idx") + 1);
  if (fname_index == NULL) {
    printf("Failed to allocate memory for index file name\n");
    return -1;
  }
  strcpy(fname_index, argv[optind]);
  strcat(fname_index, ".idx");

  FILE *out;
  if (console) out = stdout;
  else out = fopen(fname_index, "w");

  cgfile_t cgf = open_cgfile(argv[optind]);
  cgdata_t cg = {0}; int i;
  int64_t sample_idx = bgzf_tell(cgf.fh);
  for (i=0; ; ++i) {
    if (!read_cg_(&cgf, &cg)) break;
    if (snames && i < snames->n) {
      fputs(snames->array[i], out);
    } else {
      fprintf(out, "Unknown_%d", i+1);
    }
    fprintf(out, "\t%"PRId64"\n", sample_idx);
    sample_idx = bgzf_tell(cgf.fh);
  }
  fflush(out);

  if (snames && i < snames->n) {
    fprintf(stderr, "[Warning] There are more samples from the index file than detected from data.\n");
    fflush(stderr);
  }
  
  if (console) fclose(out);
  if (snames) cleanUpSnames(snames);
  free(fname_index);
  free(fname_snames);
  
  bgzf_close(cgf.fh);
  free(cg.s);
  
  
  return 0;
}

int64_t getIndex(index_t* idx, char* sample_name) {
  khiter_t k = kh_get(index, idx, sample_name);
  if (k == kh_end(idx)) {
    // Sample name not found in the hash table
    return -1;
  } else {
    // Sample name found, return the index value
    return kh_value(idx, k);
  }
}
