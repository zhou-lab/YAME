#ifndef _CGFILE_H
#define _CGFILE_H

#include "cgdata.h"
#include "snames.h"
#include "index.h"

/* cg file for reading, see cgdata_write for writing */
/* TODO, is n really necessary? */
typedef struct cgfile_t {
  BGZF *fh;
  int n;                        /* number of samples read */
} cgfile_t;

/**
 * Opens a file and returns a cgfile_t instance.
 * If the filename is "-", it will open stdin for reading. Otherwise, it opens the named file.
 * On any error opening the file, the program will exit.
 *
 * @param fname The name of the file to open.
 * @return A cgfile_t instance that represents the opened file.
 */
cgfile_t open_cgfile(char *fname);

/**
 * Reads cgdata from a cgfile_t instance.
 *
 * @param cgf The cgfile_t instance to read from.
 * @return A cgdata_t instance with the data read from the file.
 */
cgdata_t read_cg(cgfile_t *cgf);

/**
 * Reads a cgdata_t instance from a cgfile_t instance. 
 * This function is a lower-level utility for reading compressed data from a file.
 * cg memory will be reallocated.
 *
 * @param cgf The cgfile_t instance to read from.
 * @param cg The cgdata_t instance to store the read data into.
 * @return 1 if the read operation was successful and there was data to read, 0 if there was no data to read.
 */
int read_cg2(cgfile_t *cgf, cgdata_t *cg);

/**
 * Reads cgdata from a specified range in a cgfile_t instance.
 * If "end" is smaller than "beg", the program will exit with an error.
 *
 * @param cgf The cgfile_t instance to read from.
 * @param beg The beginning of the range to read from.
 * @param end The end of the range to read from.
 * @return A cgdata_v instance with the data read from the file.
 */
cgdata_v* read_cgs(cgfile_t *cgf, int64_t beg, int64_t end);

/**
 * Reads all cgdata from a cgfile_t instance. This function can be memory intensive if there are many samples.
 *
 * @param cgf The cgfile_t instance to read from.
 * @return A cgdata_v instance with all the data read from the file.
 */
cgdata_v* read_cgs_all(cgfile_t *cgf);

/**
 * Reads the first n cgdata from a cgfile_t instance.
 *
 * @param cgf The cgfile_t instance to read from.
 * @param n The number of cgdata to read from the head of the file.
 * @return A cgdata_v instance with the data read from the file.
 */
cgdata_v* read_cgs_from_head(cgfile_t *cgf, int64_t n);

/**
 * Reads the last n cgdata from a cgfile_t instance.
 *
 * @param cgf The cgfile_t instance to read from.
 * @param idx The index of the data to read.
 * @param n The number of cgdata to read from the tail of the file.
 * @return A cgdata_v instance with the data read from the file.
 */
cgdata_v* read_cgs_from_tail(cgfile_t *cgf, index_t *idx, int64_t n);

/**
 * Reads cgdata from a cgfile_t instance at specified indices.
 *
 * @param cgf The cgfile_t instance to read from.
 * @param indices The indices of the data to read.
 * @param n The number of cgdata to read.
 * @return A cgdata_v instance with the data read from the file.
 */
cgdata_v* read_cgs_with_indices(cgfile_t *cgf, const int64_t* indices, int n);

/**
 * Reads cgdata from a cgfile_t instance with specified sample names.
 * If any sample name is not found in the index, the program will exit with an error.
 *
 * @param cgf The cgfile_t instance to read from.
 * @param idx The index of the data to read.
 * @param snames The sample names to read.
 * @return A cgdata_v instance with the data read from the file.
 */
cgdata_v* read_cgs_with_snames(cgfile_t *cgf, index_t *idx, snames_t *snames);

/* void cgdata_write(char *fname_out, cgdata_t *cg, const char *mode, int verbose); */

/**
 * Writes a cgdata_t instance to a BGZF file stream.
 * 
 * This function takes a pointer to a BGZF file stream and a cgdata_t instance, 
 * and writes the data from the cgdata_t instance to the BGZF file stream. 
 * The data is expected to be in a specific format, matching the structure of the cgdata_t type.
 *
 * @param fp A pointer to the BGZF file stream to write to.
 * @param cg A pointer to the cgdata_t instance to be written.
 */
void cgdata_write1(BGZF *fp, cgdata_t *cg);

/**
 * Writes the cgdata to the specified file.
 *
 * @param fname_out The name of the output file. If NULL, output to stdout.
 * @param cg The cgdata_t instance to write to the file.
 * @param mode The mode to open the file. This should be either "w" for write mode or "a" for append mode.
 * @param verbose A flag to control verbosity. If non-zero, additional information will be printed during the write process.
 */
void cgdata_write(char *fname_out, cgdata_t *cg, const char *mode, int verbose);

#endif
