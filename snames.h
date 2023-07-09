#ifndef _SNAMES_H
#define _SNAMES_H

#include <zlib.h>
#include "wzio.h"

typedef struct {
  int n;
  char **array;
} snames_t;

/**
 * Loads the sample names from a given file.
 * If the filename is "-", it will read from stdin. 
 * If the file cannot be opened and the fatal parameter is non-zero, the program will exit with an error message.
 *
 * @param fname_snames The name of the file containing the sample names.
 * @param fatal A flag indicating whether to treat inability to open the file as a fatal error.
 *              If non-zero and the file cannot be opened, the program will exit with an error message.
 * @return A snames_t structure containing the sample names read from the file.
 *         If the file cannot be opened and the fatal parameter is zero, NULL is returned.
 *         If there is an error allocating memory, NULL is returned.
 */
static inline snames_t loadSampleNames(char* fname_snames, int fatal) {
  snames_t snames = {0};
  if (fname_snames == NULL) return snames;
  gzFile fp;
  if (strcmp(fname_snames, "-") == 0) {
    fp = gzdopen(fileno(stdin), "r");
  } else {
    fp = gzopen(fname_snames, "r");
    if (!fp && fatal) {
      fprintf(stderr, "[%s:%d] Fatal, cannot open file: %s\n",
              __func__, __LINE__, fname_snames);
      fflush(stderr);
      exit(1);
    }
  }
  
  if (fp == NULL) return snames;

  char *line = NULL;
  while (gzFile_read_line(fp, &line) > 0) {
    char *sname;
    if (line_get_field(line, 0, "\t", &sname)) {
      snames.array = realloc(snames.array, sizeof(*(snames.array)) * (snames.n + 1));
      if (snames.array == NULL) {
        fprintf(stderr, "Failed to allocate memory\n");
        fflush(stderr);
        exit(1);
      }
      snames.array[snames.n] = sname;
      snames.n++;
    }
  }
  free(line);
  line = NULL;
  gzclose(fp);

  return snames;
}

static inline void cleanSampleNames(snames_t *snames) {
  if (snames) {
    if (snames->n)
      for (int i=0; i< snames->n; ++i) {
        free(snames->array[i]);
      }
    free(snames->array);
  }
}

#endif  /* SNAMES_H */
