// SPDX-License-Identifier: AGPL-3.0-or-later
/**
 * This file is part of YAME.
 *
 * Copyright (C) 2021-present Wanding Zhou
 *
 * YAME is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * YAME is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with YAME.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "cfile.h"

#include "yame_ui.h"
/*
 * yame subset
 * ===========
 *
 * Overview
 * --------
 * This subcommand has two distinct modes:
 *
 *   (A) Sample subsetting (default): subset a multi-sample .cx by sample name.
 *   (B) Format-2 state subsetting (-s): convert a single fmt2 state vector into
 *       one binary (fmt0) track per requested state term.
 *
 * The requested names can be provided as trailing command-line arguments or
 * via -l <list.txt> (one name per line). If no names are provided, the default
 * is to take head/tail samples from the input index (see below).
 *
 *
 * Mode A: subset_samples()
 * -----------------------
 * Purpose:
 *   Extract specific samples (cdata records) from a .cx stream by name.
 *
 * Requirements:
 *   The input must be indexed (.cxi). The index maps sample name -> BGZF byte offset.
 *   Without an index, random access via bgzf_seek() is not possible, so the code
 *   errors out.
 *
 * Name selection:
 *   - If snames.n > 0 (explicit list), use that list in the given order.
 *   - Else, build a list from the input index:
 *       * if tail > 0: take last  tail samples
 *       * else:        take first head samples (head defaults to 1 if < 1)
 *
 * Extraction:
 *   For each requested sample name:
 *     - lookup offset = getIndex(idx, name)
 *     - bgzf_seek(cf.fh, offset, SEEK_SET)
 *     - read_cdata2(&cf, &c)
 *     - cdata_write1(fp_out, &c)
 *
 * Output indexing:
 *   If -o is provided (writing to a file), we generate an output index by
 *   re-reading the output file sequentially and recording bgzf_tell() offsets
 *   for each emitted record, keyed by the corresponding sample name.
 *   If output is stdout, no index is written.
 *
 *
 * Mode B: subset_fmt2_states()  [enabled by -s]
 * ------------------------------------------------
 * Purpose:
 *   Given a single format-2 (state) dataset, produce a separate binary (fmt0)
 *   mask for each requested term/state name.
 *
 * Input expectations:
 *   - Reads ONE cdata record from the input stream and requires c.fmt == '2'.
 *   - The fmt2 record is decompressed (decompress_in_situ).
 *   - fmt2_set_aux() is called if needed to populate aux->keys[] (the list of
 *     state labels/terms).
 *
 * Term resolution:
 *   For each requested term name:
 *     - scan aux->keys[] to find a unique matching key index (i_term)
 *     - error if missing or if multiple matches occur
 *
 * Mask creation:
 *   Creates a fmt0 output vector c0 with length c.n (same number of rows as input).
 *   For each row ii:
 *     - if f2_get_uint64(&c, ii) == i_term, set bit ii in c0
 *   Writes c0 as one output record per requested term.
 *
 * Output indexing:
 *   If -o is used, an output index is written mapping term name -> record offset,
 *   using the same “re-read output file and bgzf_tell() offsets” approach.
 *
 *
 * Practical notes / gotchas
 * -------------------------
 * - In default sample-subset mode, the tool is fundamentally index-driven.
 * - In -s mode, the term list is taken from the fmt2 key dictionary (aux->keys),
 *   and the output contains *multiple* fmt0 records (one per requested term).
 * - If you pass both explicit names and -H/-T, the explicit list wins because
 *   snames.n is non-zero and head/tail logic is skipped.
 */

/* BGZF write mode carrying a zlib level. bgzf_open2()/bgzf_dopen() pass the
 * mode through mode2level(), which picks up the first digit, so "w1" is a
 * level-1 writer and "w0" stores. level < 0 leaves the zlib default. */
static void write_mode(char *buf, size_t n, int level) {
  if (level >= 0 && level <= 9) snprintf(buf, n, "w%c", '0' + level);
  else snprintf(buf, n, "w");
}

void subset_fmt2_states(cfile_t cf, snames_t snames, char *fname_out,
                        int level) {
  cdata_t c = read_cdata1(&cf);
  decompress_in_situ(&c);
  if (c.fmt != '2') {
    wzfatal("To subset states, please provide a format 2 input. Give %c", c.fmt);
  }

  // output
  BGZF *fp;
  char wmode[8]; write_mode(wmode, sizeof(wmode), level);
  if (fname_out) fp = bgzf_open2(fname_out, wmode);
  else fp = bgzf_dopen(fileno(stdout), wmode);
  if (fp == NULL) {
    fprintf(stderr, "Error opening file for writing: %s\n", fname_out);
    exit(1);
  }

  /* The index is built from the writer: bgzf_tell() before a record is the
   * same virtual offset a reader recovers by reading up to it, so there is no
   * need to inflate the finished file a second time to find the offsets. */
  index_t *idx2 = fname_out ? kh_init(index) : NULL;

  if (!c.aux) fmt2_set_aux(&c);
  f2_aux_t *aux = (f2_aux_t*) c.aux;
  cdata_t c0 = {.n = c.n, .fmt = '0', .compressed=1}; // output data
  c0.s = calloc(cdata_nbytes(&c0), 1);
  for (int64_t i = 0; i<snames.n; ++i) {
    uint64_t i_term = 0; int found = 0;
    for (uint64_t j = 0; j<aux->nk; ++j) {
      if (strcmp(snames.s[i], aux->keys[j]) == 0) {
        if (found) wzfatal("Multiple match for %s.", snames.s[i]);
        else i_term = j;
        found = 1;
      }
    }
    if (!found) wzfatal("Cannot find term %s.", snames.s[i]);

    memset(c0.s, 0, cdata_nbytes(&c0));
    for (uint64_t ii = 0; ii < c.n; ++ii) {
      if (f2_get_uint64(&c, ii) == i_term) FMT0_SET(c0, ii);
    }
    bgzf_flush(fp);             /* start this record on a block boundary */
    if (idx2) insert_index(idx2, snames.s[i], bgzf_tell(fp));
    cdata_write1(fp, &c0);
  }
  free_cdata(&c);
  free_cdata(&c0);
  bgzf_close(fp);

  if (fname_out) {              // output index
    char *fname_index2 = get_fname_index(fname_out);
    FILE *out = fopen(fname_index2, "w");
    writeIndex(out, idx2);
    fclose(out);
    free(fname_index2);
    free(fname_out);
    freeIndex(idx2);
  }
}

static int cmp_int64(const void *a, const void *b) {
  int64_t x = *(const int64_t*)a, y = *(const int64_t*)b;
  return (x > y) - (x < y);
}

/**
 * Raw block passthrough.
 *
 * subset copies whole records that come out byte-identical, so when a record
 * starts on a BGZF block boundary its compressed bytes can be moved with
 * fread/fwrite -- no inflate, no deflate, no CRC. A record is aligned when its
 * index offset has a zero in-block part, which is what happens when it was
 * written to its own file and the files were concatenated (mergeCG2). Records
 * appended to a shared writer land mid-block, so this is a fast path: it
 * returns 0 without writing anything if any requested record is unaligned, and
 * the caller falls back to the decode/encode path.
 *
 * A record's compressed extent runs to the start of the record that follows it
 * *in the file*, hence the sorted copy of all index offsets; the last record
 * runs to end of file. Any embedded end-of-file marker inside that extent is
 * copied along with it, which is what concatenated stores already contain.
 */
static int subset_raw(char *fname_in, index_t *idx, snames_t snames,
                      char *fname_out, int verbose) {

  int npairs = 0;
  index_pair_t *pairs = index_pairs(idx, &npairs);
  if (npairs <= 0) { clean_index_pairs(pairs, npairs); return 0; }

  /* EVERY record must begin at a block boundary, not just the requested ones.
   * A record's extent ends where the next record begins, so a neighbour that
   * starts mid-block puts the tail of an aligned record inside the neighbour's
   * first block; copying up to that block address would silently truncate it.
   * This is arithmetic on the index alone, no I/O. */
  for (int i = 0; i < npairs; ++i)
    if (pairs[i].value & 0xFFFF) { clean_index_pairs(pairs, npairs); return 0; }
  for (int i = 0; i < snames.n; ++i)
    if (getIndex(idx, snames.s[i]) < 0) {
      clean_index_pairs(pairs, npairs); return 0;
    }

  /* record boundaries, in file order */
  int64_t *addr = malloc(npairs * sizeof(int64_t));
  for (int i = 0; i < npairs; ++i) addr[i] = pairs[i].value >> 16;
  qsort(addr, npairs, sizeof(int64_t), cmp_int64);

  FILE *in = fopen(fname_in, "rb");
  if (!in) { free(addr); clean_index_pairs(pairs, npairs); return 0; }
  if (fseeko(in, 0, SEEK_END) != 0) {
    fclose(in); free(addr); clean_index_pairs(pairs, npairs); return 0;
  }
  int64_t fsize = ftello(in);

  /* A raw copy moves bytes without reading them, so it would not notice an
   * index pointing at the wrong place -- the decode path catches that on the
   * record signature. Check both ends of every extent about to be copied,
   * before anything is written; a failure hands the job back to the decode
   * path, which reports it. Costs one block inflate per boundary. */
  {
    BGZF *fh = bgzf_open(fname_in, "r");
    if (!fh) {
      fclose(in); free(addr); clean_index_pairs(pairs, npairs); return 0;
    }
    for (int i = 0; i < snames.n; ++i) {
      int64_t beg = getIndex(idx, snames.s[i]) >> 16, end = fsize;
      for (int j = 0; j < npairs; ++j)
        if (addr[j] > beg) { end = addr[j]; break; }
      if (!cx_record_at(fh, beg << 16) ||
          (end != fsize && !cx_record_at(fh, end << 16))) {
        bgzf_close(fh); fclose(in);
        free(addr); clean_index_pairs(pairs, npairs);
        return 0;
      }
    }
    bgzf_close(fh);
  }

  FILE *out = fname_out ? fopen(fname_out, "wb") : stdout;
  if (!out) {
    fclose(in); free(addr); clean_index_pairs(pairs, npairs);
    fprintf(stderr, "Error opening file for writing: %s\n", fname_out);
    exit(1);
  }

  index_t *idx2 = fname_out ? kh_init(index) : NULL;
  int64_t written = 0;

  for (int i = 0; i < snames.n; ++i) {
    int64_t beg = getIndex(idx, snames.s[i]) >> 16, end = fsize;
    for (int j = 0; j < npairs; ++j)     /* first boundary past this record */
      if (addr[j] > beg) { end = addr[j]; break; }

    /* the output offset is a block address too, since only whole blocks are
     * ever appended, so the in-block part of the virtual offset stays 0 */
    if (idx2) insert_index(idx2, snames.s[i], written << 16);

    /* drop the concatenation's empty members, or two can end up adjacent in
     * the output and a reader stops there */
    cx_trim_empty_members(in, &beg, &end);

    if (!cx_copy_bytes(in, beg, end, out))
      wzfatal("Failed copying %s out of %s.\n", snames.s[i], fname_in);
    written += end - beg;
  }
  if (fwrite(CX_BGZF_EOF, 1, sizeof(CX_BGZF_EOF), out) != sizeof(CX_BGZF_EOF))
    wzfatal("Short write.\n");

  free(addr); clean_index_pairs(pairs, npairs);
  fclose(in);
  if (fname_out) fclose(out); else fflush(out);
  if (verbose)
    fprintf(stderr, "[subset] raw block passthrough: %"PRId64" bytes copied.\n",
            written);

  if (fname_out) {
    /* The records asked for are all in the bytes by construction, so counting
     * what was written proves nothing. What can go wrong is that a reader
     * cannot walk to them -- that is how an earlier version of this path lost
     * records silently -- so check reachability before claiming success. */
    int64_t bad = -1, nmem = cx_walk_members(fname_out, &bad);
    if (bad >= 0)
      wzfatal("[%s:%d] %s: member %"PRId64" of %"PRId64" is empty and not "
              "last, so a reader stops there and the records after it are "
              "unreachable. Refusing to leave a short file.\n",
              __func__, __LINE__, fname_out, bad, nmem);

    char *fname_index2 = get_fname_index(fname_out);
    FILE *fidx = fopen(fname_index2, "w");
    writeIndex(fidx, idx2);
    fclose(fidx);
    free(fname_index2);
    freeIndex(idx2);
  }
  return 1;
}

void subset_samples(cfile_t cf, index_t *idx, snames_t snames, char *fname_in,
                    char *fname_out, int head, int tail, int level,
                    int verbose) {

  // check if we have index
  if (!idx) {
    fprintf(stderr, "Error, the cx file needs indexing for subsetting.\n");
    fflush(stderr);
    exit(1);
  }

  // if sample names are not explicitly given, use head and tails
  if (snames.n == 0) {
    int npairs = 0;
    index_pair_t *pairs = index_pairs(idx, &npairs);
    if (tail > 0) {
      if (tail > npairs) tail = npairs;
      for (int i=0; i<tail; ++i) {
        snames.s = realloc(snames.s, (snames.n+1)*sizeof(const char*));
        snames.s[snames.n++] = strdup(pairs[npairs-tail+i].key);
      }
    } else {
      if (head < 1) head = 1;   // default to head 1
      if (head > npairs) head = npairs;
      for (int i=0; i<head; ++i) {
        snames.s = realloc(snames.s, (snames.n+1)*sizeof(const char*));
        snames.s[snames.n++] = strdup(pairs[i].key);
      }
    }
    clean_index_pairs(pairs, npairs);
  }

  /* Fast path: move the compressed bytes without touching the codec. Only
   * possible when every requested record starts on a block boundary, and only
   * when no output level was asked for, since a raw copy keeps whatever
   * compression the input has. */
  if (level < 0 && fname_in &&
      subset_raw(fname_in, idx, snames, fname_out, verbose)) return;
  if (verbose)
    fprintf(stderr, "[subset] re-encoding (records are not block-aligned, "
            "or -z was given).\n");

  // output
  BGZF *fp;
  char wmode[8]; write_mode(wmode, sizeof(wmode), level);
  if (fname_out) fp = bgzf_open2(fname_out, wmode);
  else fp = bgzf_dopen(fileno(stdout), wmode);
  if (fp == NULL) {
    fprintf(stderr, "Error opening file for writing: %s\n", fname_out);
    exit(1);
  }

  /* Index from the writer, not from a second pass over the finished file:
   * bgzf_tell() before a record gives the block address already flushed plus
   * the offset inside the block being filled, which is exactly where a reader
   * has to seek to find that record. */
  index_t *idx2 = fname_out ? kh_init(index) : NULL;

  cdata_t c = {0};              // output data
  int emitted = 0;
  for (int i=0; i<snames.n; ++i) {
    int64_t index = getIndex(idx, snames.s[i]);
    /* a name absent from the index used to abort here with a bare assertion */
    if (index < 0)
      wzfatal("[%s:%d] %s is not in the index of the input.\n",
              __func__, __LINE__, snames.s[i]);
    assert(bgzf_seek(cf.fh, index, SEEK_SET) == 0);
    read_cdata2(&cf, &c);
    if (c.n <= 0) {
      fprintf(stderr, "[%s:%d] Error, cannot find %s.\n", __func__, __LINE__, snames.s[i]);
      fflush(stderr);
      exit(1);
    }
    bgzf_flush(fp);             /* start this record on a block boundary */
    if (idx2) insert_index(idx2, snames.s[i], bgzf_tell(fp));
    cdata_write1(fp, &c);
    emitted++;
  }
  bgzf_close(fp);
  free(c.s);

  if (emitted != snames.n)
    wzfatal("[%s:%d] emitted %d records for %"PRId64" requested names.\n",
            __func__, __LINE__, emitted, (int64_t) snames.n);

  if (fname_out) {              // output index
    char *fname_index2 = get_fname_index(fname_out);
    FILE *out = fopen(fname_index2, "w");
    writeIndex(out, idx2);
    fclose(out);
    free(fname_index2);
    free(fname_out);
    freeIndex(idx2);
  }
}

static int usage(void) {
  yame_usage_head("yame subset [options] <in.cx> [sample1 sample2 ...] > out.cx");
  yame_usage_sec("Purpose:");
  yame_usage_text("Subset a multi-sample .cx by sample names (requires an index), or");
  yame_usage_text("(with -s) convert a format-2 state track into one binary track per state.");
  yame_usage_sec("Modes:");
  yame_usage_text("(A) Sample subsetting (default):");
  yame_usage_cont("Select named samples from <in.cx> and emit them in the given order.");
  yame_usage_cont("Requires <in.cx>.cxi index.");
  yame_usage_text("(B) Subset format-2 states (-s):");
  yame_usage_cont("Interpret <in.cx> as a single format-2 dataset (must be fmt2).");
  yame_usage_cont("For each requested state name, emit one format-0 bitset where");
  yame_usage_cont("bit=1 iff row state == that term.");
  yame_usage_sec("Input sample list:");
  yame_usage_text("Provide sample names either:");
  yame_usage_text("  * as trailing arguments on the command line, OR");
  yame_usage_text("  * via -l <list.txt> (one name per line).");
  yame_usage_sec("Options:");
  yame_usage_opt("-o <out.cx>", "Write output to a file. If provided, an output index (.cxi)");
  yame_usage_cont("is also generated. If omitted, writes to stdout (no index).");
  yame_usage_opt("-l <list>", "Path to sample/state list. Ignored if names are provided as");
  yame_usage_cont("trailing command-line arguments.");
  yame_usage_opt("-s", "Format-2 state filtering mode (output format 0; one record per term).");
  yame_usage_opt("-H <N>", "If no names are provided, take the first N samples from the input index.");
  yame_usage_opt("-T <N>", "If no names are provided, take the last  N samples from the input index.");
  yame_usage_opt("-z <0-9>", "Re-encode the output at this zlib level, turning off the raw");
  yame_usage_cont("copy below. Only useful when the copy cannot run, or when you");
  yame_usage_cont("want a different compression than the input has: -z1 is much");
  yame_usage_cont("faster than the default 6, -z0 stores. Any level reads back");
  yame_usage_cont("identically.");
  yame_usage_opt("-v", "Say which path ran (raw copy or re-encode).");
  yame_usage_opt("-h", "Show this help message.");
  yame_usage_sec("Notes:");
  yame_usage_text("* A subset copies whole records unchanged, so when every requested record");
  yame_usage_text("  starts on a BGZF block boundary its compressed bytes are moved verbatim");
  yame_usage_text("  -- no decompression, no re-compression. Records written to their own");
  yame_usage_text("  file and concatenated are aligned; records appended to a shared writer");
  yame_usage_text("  are not, and those fall back to decode/encode. -v says which ran.");
  yame_usage_text("* -H/-T only apply when you did NOT provide an explicit name list.");
  yame_usage_text("* -T requires an index (same as default sample subsetting).");
  yame_usage_text("* In -s mode, the input is expected to be a single fmt2 record; the output");
  yame_usage_text("  contains one fmt0 record per requested term/state.");
  fprintf(stderr, "\n");

  return 1;
}

int main_subset(int argc, char *argv[]) {

  int c0; char *fname_snames = NULL;
  int head = -1, tail = -1;
  char *fname_out = NULL;
  int filter_fmt2_states = 0;
  int level = -1;               // -1 = zlib default (6), and lets raw copy run
  int verbose = 0;
  while ((c0 = getopt(argc, argv, "o:l:sH:T:z:vh"))>=0) {
    switch (c0) {
    case 'o': fname_out = strdup(optarg); break;
    case 'l': fname_snames = strdup(optarg); break;
    case 's': filter_fmt2_states = 1; break;
    case 'H': head = atoi(optarg); break;
    case 'T': tail = atoi(optarg); break;
    case 'z': {                 /* not atoi: "abc" would pass as level 0,
                                 * i.e. silently store instead of compress */
      char *endp = NULL;
      long v = strtol(optarg, &endp, 10);
      if (endp == optarg || *endp || v < 0 || v > 9)
        wzfatal("-z takes a zlib level 0-9, given \"%s\".\n", optarg);
      level = (int)v;
      break;
    }
    case 'v': verbose = 1; break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c0);
    }
  }

  if (argc < optind + 1) { 
    usage(); 
    wzfatal("Please supply input file and output file.\n"); 
  }

  // input
  char *fname_in = argv[optind];
  cfile_t cf = open_cfile(fname_in);
  char *fname_index = get_fname_index(fname_in);
  index_t *idx = loadIndex(fname_index);
  free(fname_index);
  optind++;

  // get sample names
  snames_t snames = {0};
  if (optind < argc) {          // sample names from command line
    for(int i = optind; i < argc; ++i) {
      snames.s = realloc(snames.s, (snames.n+1)*sizeof(const char*));
      snames.s[snames.n++] = strdup(argv[i]);
    }
  } else {                      // from a file list
    snames = loadSampleNames(fname_snames, 1);
  }

  if (filter_fmt2_states) {
    subset_fmt2_states(cf, snames, fname_out, level);
  } else {
    subset_samples(cf, idx, snames, fname_in, fname_out, head, tail, level,
                   verbose);
  }
  
  // clean up
  bgzf_close(cf.fh);
  cleanIndex(idx);
  cleanSampleNames(&snames);
  
  return 0;
}

