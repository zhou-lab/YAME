#ifndef _WZIO_H
#define _WZIO_H

#include <zlib.h>
#include "wzmisc.h"
#include "wvec.h"

#define wzclose gzclose

/**
 * Opens a file for reading, whether it is a plain text file or a gzipped file.
 *
 * @param path The path to the file to be opened.
 *             Use "-" to read from standard input.
 * @fatal if non-zero, the function will exit upon failure
 * @return A gzFile handle for the opened file.
 * @throws An error message and exits the program if the file cannot be opened.
 */
static inline gzFile wzopen(char *path, int fatal) {
  gzFile fh;
  if (strcmp(path, "-") == 0) {
    fh = gzdopen(fileno(stdin), "r");
  } else {
    fh = gzopen(path, "r");
    if (!fh && fatal) {
      fprintf(stderr, "[%s:%d] Fatal, cannot open file: %s\n",
              __func__, __LINE__, path);
      fflush(stderr);
      exit(1);
    }
  }
  return fh;
}

/**
 * Opens a file for writing.
 *
 * @param path The path to the file to be opened.
 *             Use NULL to write to standard output.
 * @return A FILE pointer for the opened file.
 * @throws An error message and exits the program if the file cannot be opened.
 */
static inline FILE *wzopen_out(char *path) {
   FILE *out;
   if (path) {
      out = fopen(path, "w");
      if (!out) {
         fprintf(stderr, "[%s:%d] Fatal, cannot open file: %s\n",
                 __func__, __LINE__, path);
         fflush(stderr);
         exit(1);
      }
   } else {
      out = stdout;
   }
   return out;
}

/**
 * Reads a line from a gzFile handle.
 *
 * @usage char *line; gzFile_read_line(fh, &line); free(line);
 * @param fh The gzFile handle to read from.
 * @param s  A pointer to a character pointer that will hold the read line.
 *           Either NULL or previously allocated c-string
 *           The memory for the line is allocated by the function and should be freed by the caller.
 * @return 1 if a line is successfully read, 0 if end-of-file is reached.
 * @throws An error message and exits the program if the string pointer is NULL.
 */
static inline int gzFile_read_line(gzFile fh, char **s) {

  if (s == NULL) {
    fprintf(stderr, "[%s:%d] Fatal, empty string construct.\n", __func__, __LINE__);
    fflush(stderr);
    exit(1);
  }
  
  /* reset s */
  int m = 10, l = 0;            /* memory and string length */
  *s = realloc(*s, m);

  /* read until '\n' or EOF */
  while (1) {
    int c = gzgetc(fh);
    if (l > m-2) { m <<= 1; *s = realloc(*s, m); }
    if (c == '\n') {(*s)[l] = '\0'; return 1;}
    if (c == EOF) {(*s)[l] = '\0'; return 0;}
    (*s)[l++] = c;
  }
  return 0;                     /* should not come here */
}

/**
 * Counts the number of lines in a gzFile handle.
 *
 * @param fh The gzFile handle to count lines in.
 * @return The number of lines in the file.
 */
static inline int gzFile_count_lines(gzFile fh) {

  int n = 0;
  while(1) {
    int c = gzgetc(fh);
    if (c == '\n') n++;
    if (c == EOF) return n+1;
  }
  return 0;                     /* should not come here */
}


/**
 * Retrieves a field from a line based on its index and separator.
 *
 * @param line         The input line to extract the field from.
 * @param field_index  The index of the desired field (0-based).
 * @param sep          The separator character(s) used to split the line into fields.
 * @param field        A pointer to a character pointer that will hold the extracted field.
 *                     The memory for the field is allocated by the function and should be freed by the caller.
 * @return 1 if the field is successfully retrieved, 0 if there are not enough fields.
 */
static inline int line_get_field(const char *line, int field_index, const char *sep, char **field) {

  char *working = calloc(strlen(line) + 1, sizeof(char));
  strcpy(working, line);
  char *tok;

  tok = strtok(working, sep);
  int i;
  for (i=0; i<field_index; ++i)
    tok = strtok(NULL, sep);

  if (tok == NULL) {            /* not enough fields */
    free(working);
    return 0;
  }

  *field = strdup(tok);
  free(working);
  return 1;
}

/********************************
 ** Get all fields of one line **
 ********************************
Usage:
   char **fields; int nfields;
   line_get_fields("my line", " ", &fields, &nfields);
   free_fields(fields, nfields);

   Note: separators/delimiters are not merged - the most likely use-case. */

/**
 * Retrieves all fields from a line based on a separator.
 *
 * @param line     The input line to extract the fields from.
 * @param sep      The separator character(s) used to split the line into fields.
 * @param fields   A pointer to a character pointer array that will hold the extracted fields.
 *                 The memory for the fields and array is allocated by the function and should be freed by the caller.
 * @param nfields  A pointer to an integer that will store the number of extracted fields.
 *
 * Example Usage:
 * ```c
 * char **fields;
 * int nfields;
 * const char *line = "John,Doe,42";
 * line_get_fields(line, ",", &fields, &nfields);
 *
 * // Access the extracted fields
 * for (int i = 0; i < nfields; i++) {
 *     printf("Field %d: %s\n", i, fields[i]);
 * }
 *
 * // Free the memory allocated by line_get_fields
 * free_fields(fields, nfields);
 * ```
 */
#define free_fields(fds, nfds) free_char_array(fds, nfds)
static inline void line_get_fields(const char *line, const char *sep, char ***fields, int *nfields) {

  *nfields = 1;
  const char *s = line;
  while ((s = strpbrk(s, sep)) != NULL) { (*nfields)++; s++; }

  *fields = calloc(*nfields, sizeof(char *));
  char *working = calloc(strlen(line) + 1, sizeof(char));
  strcpy(working, line);
  char *tok; int i;

  tok = strtok(working, sep);
  for (i=0; tok != NULL; ++i) {
    (*fields)[i] = strdup(tok);
    tok = strtok(NULL, sep);
  }
  free(working);
}

/**
 * Retrieves a specific number of fields from a line based on a separator.
 *
 * @param line     The input line to extract the fields from.
 * @param sep      The separator character(s) used to split the line into fields.
 * @param fields   A pointer to a character pointer array that will hold the extracted fields.
 *                 The memory for the fields and array is allocated by the function and should be freed by the caller.
 * @param nfields  A pointer to an integer that specifies the expected number of extracted fields.
 *                 If nfields < 0, it will be set to the actual number of extracted fields.
 * @param aux      A pointer to a character pointer that will store auxiliary memory for parsing.
 *                 The memory is allocated by the function and should be freed by the caller.
 * @throws An error message and exits the program if the number of fields extracted does not match the expected number.
 *
 * Example Usage:
 * ```c
 * char **fields;
 * int nfields = 3;
 * const char *line = "John,Doe,42";
 * char *aux = NULL;
 * line_get_fields2(line, ",", &fields, &nfields, &aux);
 *
 * // Access the extracted fields
 * for (int i = 0; i < nfields; i++) {
 *     printf("Field %d: %s\n", i, fields[i]);
 * }
 *
 * // Free the memory allocated by line_get_fields2
 * free_fields(fields, nfields);
 * free(aux);
 * ```
 */
static inline void line_get_fields2(
  const char *line, const char *sep, char ***fields, int *nfields, char **aux) {

  int nfields_ = 1;
  const char *s = line;
  while ((s = strpbrk(s, sep)) != NULL) { nfields_++; s++; }

  if (*nfields < 0) {
    *nfields = nfields_;
    *fields = calloc(*nfields, sizeof(char*));
  } else if (*nfields != nfields_) {
    fprintf(stderr, "Wrong field number %d (expecting %d).\n", nfields_, *nfields);
    fflush(stderr);
    exit(1);
  }

  *aux = realloc(*aux, (strlen(line) + 1) * sizeof(char));
  strcpy(*aux, line);
  char *tok; int i;

  tok = strtok(*aux, sep);
  for (i=0; tok != NULL; ++i) {
    (*fields)[i] = realloc((*fields)[i], strlen(tok)+1);
    strcpy((*fields)[i], tok);
    tok = strtok(NULL, sep);
  }
}

#endif /* _WZIO_H */
