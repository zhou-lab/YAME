#ifndef _CFILE_H
#define _CFILE_H

#include "cdata.h"
#include "snames.h"
#include "index.h"

/* cfile for reading, see cdata_write for writing */
typedef struct cfile_t {
  BGZF *fh;
  int n;                        /* number of samples read */
} cfile_t;

/**
 * Opens a file and returns a cfile_t instance.
 * If the filename is "-", it will open stdin for reading. Otherwise, it opens the named file.
 * On any error opening the file, the program will exit.
 *
 * @param fname The name of the file to open.
 * @return A cfile_t instance that represents the opened file.
 */
cfile_t open_cfile(char *fname);

/**
 * Reads cdata from a cfile_t instance.
 *
 * @param cf The cfile_t instance to read from.
 * @return A cdata_t instance with the data read from the file.
 */
cdata_t read_cdata1(cfile_t *cf);

/**
 * Reads a cdata_t instance from a cfile_t instance. 
 * This function is a lower-level utility for reading compressed data from a file.
 * c memory will be reallocated.
 *
 * @param cf The cfile_t instance to read from.
 * @param c The cdata_t instance to store the read data into.
 * @return 1 if the read operation was successful and there was data to read, 0 if there was no data to read.
 */
int read_cdata2(cfile_t *cf, cdata_t *c);

DEFINE_VECTOR(cdata_v, cdata_t)

/**
 * Reads cdata from a specified range in a cfile_t instance.
 * If "end" is smaller than "beg", the program will exit with an error.
 *
 * @param cf The cfile_t instance to read from.
 * @param beg The beginning of the range to read from.
 * @param end The end of the range to read from.
 * @return A cdata_v instance with the data read from the file.
 */
cdata_v* read_cdata(cfile_t *cf, int64_t beg, int64_t end);

/**
 * Reads all cdata from a cfile_t instance. This function can be memory intensive if there are many samples.
 *
 * @param cf The cfile_t instance to read from.
 * @return A cdata_v instance with all the data read from the file.
 */
cdata_v* read_cdata_all(cfile_t *cf);

/**
 * Reads the first n cdata from a cfile_t instance.
 *
 * @param cf The cfile_t instance to read from.
 * @param n The number of cdata to read from the head of the file.
 * @return A cdata_v instance with the data read from the file.
 */
cdata_v* read_cdata_from_head(cfile_t *cf, int64_t n);

/**
 * Reads the last n cdata from a cfile_t instance.
 *
 * @param cf The cfile_t instance to read from.
 * @param idx The index of the data to read.
 * @param n The number of cdata to read from the tail of the file.
 * @return A cdata_v instance with the data read from the file.
 */
cdata_v* read_cdata_from_tail(cfile_t *cf, index_t *idx, int64_t n);

/**
 * Reads cdata from a cfile_t instance at specified indices.
 *
 * @param cf The cfile_t instance to read from.
 * @param indices The indices of the data to read.
 * @param n The number of cdata to read.
 * @return A cdata_v instance with the data read from the file.
 */
cdata_v* read_cdata_with_indices(cfile_t *cf, const int64_t* indices, int n);

/**
 * Reads cdata from a cfile_t instance with specified sample names.
 * If any sample name is not found in the index, the program will exit with an error.
 *
 * @param cf The cfile_t instance to read from.
 * @param idx The index of the data to read.
 * @param snames The sample names to read.
 * @return A cdata_v instance with the data read from the file.
 */
cdata_v* read_cdata_with_snames(cfile_t *cf, index_t *idx, snames_t *snames);

/**
 * Writes a cdata_t instance to a BGZF file stream.
 * 
 * This function takes a pointer to a BGZF file stream and a cdata_t instance, 
 * and writes the data from the cdata_t instance to the BGZF file stream. 
 * The data is expected to be in a specific format, matching the structure of the cdata_t type.
 *
 * @param fp A pointer to the BGZF file stream to write to.
 * @param c A pointer to the cdata_t instance to be written.
 */
void cdata_write1(BGZF *fp, cdata_t *c);

/**
 * Writes the cdata to the specified file.
 *
 * @param fname_out The name of the output file. If NULL, output to stdout.
 * @param c The cdata_t instance to write to the file.
 * @param mode The mode to open the file. This should be either "w" for write mode or "a" for append mode.
 * @param verbose A flag to control verbosity. If non-zero, additional information will be printed during the write process.
 */
void cdata_write(char *fname_out, cdata_t *c, const char *mode, int verbose);

#endif
