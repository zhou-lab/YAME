#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "kstring.h"
#include "cfile.h"

char *get_fname_index(const char *fname_cx) {
  char *fname_index = NULL;
  fname_index = malloc(strlen(fname_cx) + strlen(".idx") + 1);
  if (fname_index == NULL) {
    printf("Failed to allocate memory for index file name\n");
    return NULL;
  }
  strcpy(fname_index, fname_cx);
  strcat(fname_index, ".idx");
  return fname_index;
}

index_t* loadIndex(char* fname_index) {
  gzFile file = wzopen(fname_index, 0);
  if (file == NULL) {
    /* fprintf(stderr, "Failed to open index file: %s\n", fname_index); */
    return NULL;
  }

  index_t* idx = kh_init(index);
  if (idx == NULL) {
    printf("Failed to create hash table\n");
    wzclose(file);
    return NULL;
  }

  char* line = NULL;
  while (gzFile_read_line(file, &line) > 0) {
    char* sname;
    if (line_get_field(line, 0, "\t", &sname)) {
      char* index_str;
      if (line_get_field(line, 1, "\t", &index_str)) {
        int ret;
        khiter_t k = kh_put(index, idx, sname, &ret);
        if (ret == -1) {
          printf("Failed to insert value into hash table\n");
          wzclose(file);
          free(line);
          kh_destroy(index, idx);
          free(idx);
          return NULL;
        }
        kh_value(idx, k) = strtoll(index_str, NULL, 10);
        free(index_str);
      }
    }
  }
  free(line);
  wzclose(file);
  return idx;
}

static int64_t last_address(index_t *idx) {
  int64_t max_addr = 0; int64_t addr;
  kh_foreach_value(idx, addr, {
      if (addr > max_addr) max_addr = addr;
    });
  return max_addr;
}

/* static int64_t first_address(index_t *idx) { */
/*   int64_t min_addr = 0; int64_t addr; */
/*   kh_foreach_value(idx, addr, { */
/*       if (addr < min_addr) min_addr = addr; */
/*     }); */
/*   return min_addr; */
/* } */

index_t *insert_index(index_t *idx, char *sname, int64_t addr) {
  if (getIndex(idx, sname) >= 0) {
    fprintf(stderr, "[Error] Sample name %s already exists in index.\n", sname);
    exit(1);
  }
  int ret;
  khiter_t k = kh_put(index, idx, sname, &ret);
  if (ret == -1) {
    printf("Failed to insert value into hash table\n");
    exit(1);
  }
  kh_value(idx, k) = addr;
  return idx;
}

static index_t* append_index(index_t *idx, cfile_t *cf, char* sname_to_append) {

  int64_t addr;
  cdata_t c = {0};
  if (kh_size(idx) == 0) {      /* first item in index */
    addr = bgzf_tell(cf->fh);
  } else {
    assert(bgzf_seek(cf->fh, last_address(idx), SEEK_SET) == 0);
    read_cdata2(cf, &c);         /* read past the last c data block */
    addr = bgzf_tell(cf->fh);
  }
  read_cdata2(cf, &c);      /* make sure we do have one additional data block */
  if (c.n > 0) {
    idx = insert_index(idx, sname_to_append, addr);
  } else {
    fprintf(stderr, "Failed to detect additional data.\n");
  }
  return idx;
}

static int comparePairs(const void* a, const void* b) {
    const index_pair_t* pairA = (const index_pair_t*)a;
    const index_pair_t* pairB = (const index_pair_t*)b;
    
    if (pairA->value < pairB->value) return -1;
    else if (pairA->value > pairB->value) return 1;
    else return 0;
}

/* return sorted key-value pairs */
index_pair_t *index_pairs(index_t *idx, int *n) {

  index_pair_t* pairs = (index_pair_t*)malloc(kh_size(idx) * sizeof(index_pair_t));
  *n = 0;

  // Iterate over key-value pairs and store them in the array
  char *key;
  int64_t value;
  kh_foreach(idx, key, value, {
      pairs[*n].key = strdup(key);
      pairs[*n].value = value;
      (*n)++;
    });

  // Sort the array in ascending order of values
  qsort(pairs, *n, sizeof(index_pair_t), comparePairs);
  return pairs;
}

void clean_index_pairs(index_pair_t *idx_pairs, int n) {
  for (int i=0; i<n; ++i) {
    free(idx_pairs[i].key);
  }
  free(idx_pairs);
}

/* index_pair_t* load_index_pairs(char *fname_cx, int *n) { */
/*   char *fname_index = get_fname_index(fname_cx); */
/*   index_t *idx = loadIndex(fname_index); */
/*   if (!idx) { */
/*     *n = 0; */
/*     return NULL; */
/*   } */
/*   index_pair_t *idx_pairs = index_pairs(idx, n); */
/*   cleanIndex(idx); free(fname_index); */
/*   return idx_pairs; */
/* } */

snames_t loadSampleNamesFromIndex(char *fname) {
  char *fname_index = get_fname_index(fname);
  index_t *idx = loadIndex(fname_index);
  snames_t snames = {0};
  free(fname_index);
  if (!idx) return snames;

  index_pair_t *idx_pairs = index_pairs(idx, &(snames.n));
  snames.s = calloc(snames.n, sizeof(char*));
  for (int i=0; i<snames.n; ++i) snames.s[i] = idx_pairs[i].key;
  free(idx_pairs);          // ownership of keys are transfered to snames.s
  cleanIndex(idx);
  return snames;
}

void writeIndex(FILE *fp, index_t *idx) {

  int n;
  index_pair_t *pairs = index_pairs(idx, &n);
  for (int i = 0; i < n; i++) {
    const char* key = pairs[i].key;
    int64_t value = pairs[i].value;
    fprintf(fp, "%s\t%"PRId64"\n", key, value);
  }
  
  // Clean up the hash map and array
  free(pairs);
}
  
static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: yame index [options] <in.cx>\n");
  fprintf(stderr, "The index file name default to <in.cx>.idx\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -s [file path]   tab-delimited sample name list (use first column) \n");
  fprintf(stderr, "    -1 [sample name] add one sample to the end of the index\n");
  fprintf(stderr, "    -c               output index to console\n");
  fprintf(stderr, "    -h               This help\n");
  fprintf(stderr, "\n");

  return 1;
}

int main_index(int argc, char *argv[]) {

  int c0, console = 0;
  char *fname_snames = NULL;
  char *sname_to_append = NULL;
  while ((c0 = getopt(argc, argv, "cs:1:h"))>=0) {
    switch (c0) {
    case 'c': console = 1; break;
    case 's': fname_snames = strdup(optarg); break;
    case '1': sname_to_append = strdup(optarg); break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c0);
    }
  }

  if (optind + 1 > argc) { 
    usage(); 
    wzfatal("Please supply input file.\n"); 
  }

  /* load index */
  char *fname_index = get_fname_index(argv[optind]);
  cfile_t cf = open_cfile(argv[optind]);
  cdata_t c = {0};
  if (sname_to_append) {        /* append new sample to existing index */
    
    index_t *idx = loadIndex(fname_index);
    if (idx) {
      idx = append_index(idx, &cf, sname_to_append);
    } else {
      idx = kh_init(index);
      idx = append_index(idx, &cf, sname_to_append);
    }

    FILE *out;
    if (console) out = stdout;
    else out = fopen(fname_index, "w");
    writeIndex(out, idx);
    if (!console) fclose(out);
    
  } else {                      /* index all samples */
    
    snames_t snames = loadSampleNames(fname_snames, 1);
    int n=0; kstring_t *sname_v = NULL; int64_t *addr_v = NULL;
    index_t* idx = kh_init(index);

    if (snames.n >0) {               /* sample names is given */
      int64_t addr = bgzf_tell(cf.fh);
      for (int i=0; i< snames.n; ++i) {
        if (!read_cdata2(&cf, &c)) {
          fprintf(stderr, "[Error] Data is shorter than the sample name list.\n");
          fflush(stderr);
          exit(1);
        }
        insert_index(idx, snames.s[i], addr);
        addr = bgzf_tell(cf.fh);
      }
      
    } else {                    /* sample names are unknown */

      for (n=0; ; ++n) {
        int64_t addr = bgzf_tell(cf.fh);
        if (!read_cdata2(&cf, &c)) break;
        sname_v = realloc(sname_v, sizeof(kstring_t)*(n+1));
        addr_v = realloc(addr_v, sizeof(int64_t)*(n+1));
        memset(&sname_v[n], 0, sizeof(kstring_t));
        ksprintf(&sname_v[n], "Unknown_%d", n+1);
        addr_v[n] = addr;
      }
      
      for (int i=0; i<n; ++i) {
        insert_index(idx, sname_v[i].s, addr_v[i]);
      }
    }
    
    FILE *out;
    if (console) out = stdout;
    else out = fopen(fname_index, "w");
    writeIndex(out, idx);
    if (!console) fclose(out);

    if (snames.n >0) {
      cleanSampleNames(&snames);
    } else if (n) {
      for (int i=0; i<n; ++i) {
        free(sname_v[i].s);
      }
    }
    if (n) { free(addr_v); free(sname_v); }
  }

  free(fname_snames);
  free(fname_index);
  bgzf_close(cf.fh);
  free(c.s);
  return 0;
}

int64_t getIndex(index_t* idx, char* sname) {
  khiter_t k = kh_get(index, idx, sname);
  if (k == kh_end(idx)) {
    // Sample name not found in the hash table
    return -1;
  } else {
    // Sample name found, return the index value
    return kh_value(idx, k);
  }
}
