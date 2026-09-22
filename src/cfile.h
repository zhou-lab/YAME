// SPDX-License-Identifier: LicenseRef-CHOP-Academic-BSD-2-Clause
// 2-Clause BSD for academic and non-profit research use; commercial use by
// inquiry to zhouw3@chop.edu (see LICENSE).
/**
 * This file is part of YAME.
 *
 * Copyright (C) 2021-present The Children's Hospital of Philadelphia
 *
 * Use of this software is available to academic and non-profit institutions
 * for research purposes under the 2-Clause BSD License; for use or transfers
 * to commercial entities, inquire with Dr. Wanding Zhou at zhouw3@chop.edu.
 * See the LICENSE file at the root of the repository for the full terms.
 */

#ifndef _CFILE_H
#define _CFILE_H

#include "cdata.h"
#include "snames.h"
#include "index.h"

/** The header design, 17 bytes
 *  uint64_t: signature, used for validation
 *  uint8_t: format (0=vec; 1=rle)
 *  uint64_t: length (n_cs or n_bytes for rle)
 */
/* cfile for reading, see cdata_write for writing */
typedef struct cfile_t {
  BGZF *fh;
  int n;                        /* number of samples read */
  const char *fname;            /* what open_cfile() was given, for messages;
                                 * points at the caller's string, not owned */
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
 * Outcome of one record-read attempt. The core reports; the caller sets
 * policy. YAME's own commands treat everything past CX_READ_END as fatal
 * (see read_cdata2 below); a reader whose CX data is a PREFIX of a larger
 * file passes a limit instead and never sees the other statuses.
 */
typedef enum {
  CX_READ_OK = 0,     /* a record was read                                  */
  CX_READ_END,        /* clean end: a zero-length read at a record boundary,
                       * or the caller's limit was reached                   */
  CX_READ_NOT_CX,     /* the next bytes are not a BGZF block header          */
  CX_READ_TRUNCATED,  /* a block or record ends short of what it promises    */
  /* NOT_CX and TRUNCATED are best-effort labels taken from the BGZF error
   * code, and the two cannot always be told apart: a truncation that leaves a
   * stub shorter than a block header is indistinguishable from a foreign tail.
   * Do not use NOT_CX as "my data ends here" -- pass a limit and read END. */
  CX_READ_BADSIG,     /* a record header without CDSIG                       */
  CX_READ_NOMEM       /* the record's bytes could not be allocated           */
} cx_read_t;

/* "read to the end of the stream", the value every YAME command passes */
#define CX_NO_LIMIT ((int64_t)-1)
/* enough for the longest message cx_read_record() formats */
#define CX_ERRBUF 256

/**
 * Read one record and say what happened, without ever exiting.
 *
 * @param cf    the open file
 * @param c     the record, overwritten
 * @param limit raw file offset at which this reader's CX data ends; no block
 *              starting at or after it is pulled, and the read reports
 *              CX_READ_END instead. CX_NO_LIMIT reads to end of stream.
 *              Meaningful only for a seekable input: stdin must pass
 *              CX_NO_LIMIT.
 * @param err   buffer for the message describing a failure status, or NULL
 * @param errn  size of `err`
 *
 * A bundle whose CX portion is a prefix (a methscope MSBNDL1 file) passes the
 * container offset as `limit`. Without it the reader walks into the container
 * and reports CX_READ_NOT_CX, which is how a documented bundle read turned
 * into a fatal in v1.40.
 */
cx_read_t cx_read_record(cfile_t *cf, cdata_t *c, int64_t limit,
                         char *err, size_t errn);

/**
 * Raw block passthrough helpers.
 *
 * A record that begins on a BGZF block boundary can be moved between files by
 * copying its compressed bytes, with no inflate and no deflate. The unit that
 * can be copied is a whole block, so this only applies when the record starts
 * a block; otherwise the first block also holds the previous record's tail.
 * Producers that write each record through its own BGZF stream, or that flush
 * before each record, leave every record aligned.
 */

/** The 28-byte empty BGZF member that terminates a file. */
extern const uint8_t CX_BGZF_EOF[28];

/** 1 if the virtual offset is block-aligned and a record signature sits there.
 *  Reading it costs one block inflate. Used to catch a stale index before a
 *  raw copy, which otherwise moves bytes without ever looking at them. */
int cx_record_at(BGZF *fh, int64_t voffset);

/** Narrow [*beg,*end) to the record's own blocks, dropping the empty BGZF
 *  members a concatenated store leaves at either end. Copying those through
 *  can make a reader stop early -- see the note on the definition. Pass
 *  at_eof non-zero only when *end is end of file, the one case where an extent
 *  can close with an empty member. */
void cx_trim_empty_members(FILE *in, int64_t *beg, int64_t *end, int at_eof);

/** Copy [beg,end) of `in` to `out`. Both are plain handles; beg and end are
 *  byte offsets, which for an aligned record are block boundaries. */
int cx_copy_bytes(FILE *in, int64_t beg, int64_t end, FILE *out);

/** 1 if an empty BGZF member sits at this byte offset. Checked at record
 *  boundaries after a raw copy: an empty member there can strand every record
 *  that follows, present in the bytes but unreachable to a reader. */
int cx_empty_member_at(const char *fname, int64_t off);

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

/**
 * stdout as a BGZF writer, refusing a terminal.
 *
 * Every command that writes a CX stream with no -o sends it here. Written to a
 * terminal it is compressed binary: it garbles the display and can leave the
 * terminal in a state the user has to reset. A reader copying
 * `yame pairwise -H 1 -c 5 -d 0.2 a.cg b.cg` off the documentation hit exactly
 * that. Redirected or piped it behaves as before, because then the bytes have
 * somewhere to go.
 */
BGZF *yame_bgzf_stdout(const char *mode, const char *cmd);

#endif
