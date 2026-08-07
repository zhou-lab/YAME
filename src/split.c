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
static int usage() {
  yame_usage_head("yame split [options] <in.cx> out_prefix");
  yame_usage_sec("Options:");
  yame_usage_opt("-v", "verbose; also says which path ran");
  yame_usage_opt("-s", "sample name list");
  yame_usage_opt("-h", "This help");
  yame_usage_sec("Notes:");
  yame_usage_text("* Each output holds one record unchanged, so when the input is indexed");
  yame_usage_text("  and its records start on BGZF block boundaries, the compressed bytes");
  yame_usage_text("  are copied out verbatim -- no decompression, no re-compression. An");
  yame_usage_text("  unindexed or unaligned input falls back to decode/encode.");
  fprintf(stderr, "\n");

  return 1;
}

/* Build the output name for record i, matching the streaming path below. */
static char *out_name(char *prefix, char **snames, int snames_n, int i) {
  char *tmp;
  if (snames_n) {
    if (prefix) {
      tmp = malloc(strlen(prefix)+strlen(snames[i])+1000);
      sprintf(tmp, "%s%s.cx", prefix, snames[i]);
    } else {
      tmp = malloc(strlen(snames[i])+1000);
      sprintf(tmp, "%s.cx", snames[i]);
    }
  } else {
    if (prefix) {
      tmp = malloc(strlen(prefix) + 1000);
      sprintf(tmp, "%s_split_%i.cx", prefix, i+1);
    } else {
      tmp = malloc(1000);
      sprintf(tmp, "split_%i.cx", i+1);
    }
  }
  return tmp;
}

static int cmp_int64(const void *a, const void *b) {
  int64_t x = *(const int64_t*)a, y = *(const int64_t*)b;
  return (x > y) - (x < y);
}

/**
 * Raw block passthrough for split.
 *
 * Every output file holds exactly one record and starts a fresh BGZF stream,
 * so there is nothing to align on the way out. All that is needed is the
 * record's compressed extent in the input, which the index supplies: from its
 * block address to the next record's, or to end of file for the last one.
 *
 * Returns 0 without writing anything when the input has no index, when any
 * record starts mid-block, or when a -s list does not cover every record; the
 * caller then streams as before.
 */
static int split_raw(char *fname_in, char *prefix, char **snames, int snames_n,
                     int verbose) {

  char *fname_index = get_fname_index(fname_in);
  index_t *idx = loadIndex(fname_index);
  free(fname_index);
  if (!idx) return 0;

  int npairs = 0;
  index_pair_t *pairs = index_pairs(idx, &npairs);
  if (npairs <= 0 || (snames_n && snames_n < npairs)) {
    clean_index_pairs(pairs, npairs); cleanIndex(idx); return 0;
  }

  int64_t *addr = malloc(npairs * sizeof(int64_t));
  for (int i = 0; i < npairs; ++i) addr[i] = pairs[i].value;
  qsort(addr, npairs, sizeof(int64_t), cmp_int64);   /* file order */

  BGZF *fh = bgzf_open(fname_in, "r");
  if (!fh) { free(addr); clean_index_pairs(pairs, npairs); cleanIndex(idx); return 0; }
  for (int i = 0; i < npairs; ++i) {
    if (!cx_record_at(fh, addr[i])) {   /* unaligned, or a stale index */
      bgzf_close(fh); free(addr);
      clean_index_pairs(pairs, npairs); cleanIndex(idx);
      return 0;
    }
  }
  bgzf_close(fh);

  FILE *in = fopen(fname_in, "rb");
  if (!in || fseeko(in, 0, SEEK_END) != 0) {
    if (in) fclose(in);
    free(addr); clean_index_pairs(pairs, npairs); cleanIndex(idx);
    return 0;
  }
  int64_t fsize = ftello(in);

  for (int i = 0; i < npairs; ++i) {
    int64_t beg = addr[i] >> 16;
    int64_t end = (i+1 < npairs) ? (addr[i+1] >> 16) : fsize;
    cx_trim_empty_members(in, &beg, &end);
    char *tmp = out_name(prefix, snames, snames_n, i);
    FILE *out = fopen(tmp, "wb");
    if (!out) wzfatal("Cannot write %s.\n", tmp);
    if (!cx_copy_bytes(in, beg, end, out)) wzfatal("Failed copying %s.\n", tmp);
    if (fwrite(CX_BGZF_EOF, 1, sizeof(CX_BGZF_EOF), out) != sizeof(CX_BGZF_EOF))
      wzfatal("Short write to %s.\n", tmp);
    fclose(out);
    if (verbose) fprintf(stderr, "%s\n", tmp);
    free(tmp);
  }
  fclose(in);
  free(addr); clean_index_pairs(pairs, npairs); cleanIndex(idx);
  if (verbose)
    fprintf(stderr, "[split] raw block passthrough: %d records copied.\n", npairs);
  return 1;
}

int main_split(int argc, char *argv[]) {

  int c, verbose = 0; char *fname_snames = NULL;
  while ((c = getopt(argc, argv, "s:vh"))>=0) {
    switch (c) {
    case 'v': verbose = 1; break;
    case 's': fname_snames = strdup(optarg); break;
    case 'h': return usage(); break;
    default: usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }

  if (optind + 1 > argc) { 
    usage(); 
    wzfatal("Please supply input file.\n");
  }

  char *fname_in = strcmp(argv[optind], "-") == 0 ? NULL : argv[optind];
  cfile_t cf = open_cfile(argv[optind++]);
  char *prefix = NULL;
  if (optind + 1 == argc) {
    prefix = argv[optind];
  }
  
  char **snames = NULL; int snames_n = 0;
  if (fname_snames) {
    gzFile fh = wzopen(fname_snames, 1);
    char *line = NULL;
    char **fields; int nfields;
    while (gzFile_read_line(fh, &line)>0) {
      line_get_fields(line, "\t", &fields, &nfields);
      snames = realloc(snames, (snames_n+1)*sizeof(char*));
      snames[snames_n] = strdup(fields[0]);
      ++snames_n;
      free_fields(fields, nfields);
    }
    free(line);
    wzclose(fh);
    free(fname_snames);
  }
  
  /* Fast path: copy each record's compressed bytes straight out. */
  if (fname_in && split_raw(fname_in, prefix, snames, snames_n, verbose)) {
    bgzf_close(cf.fh);
    return 0;
  }
  if (verbose)
    fprintf(stderr, "[split] re-encoding (no index, or records are not "
            "block-aligned).\n");

  int i = 0;
  for (i=0; ; ++i) {
    cdata_t c = read_cdata1(&cf);
    if (c.n == 0) break;
    char *tmp = out_name(prefix, snames, snames_n, i);
    cdata_write(tmp, &c, "wb", verbose);
    free(tmp);
    free(c.s);
  }
  
  return 0;
}
