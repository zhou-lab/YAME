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

int read_cdata2(cfile_t *cf, cdata_t *c) {
  c->n = 0;
  uint64_t sig;
  int64_t size;
  if (cf->fh->block_length == 0) bgzf_read_block(cf->fh); /* somehow this is needed for concat'ed bgzipped files */
  size = bgzf_read(cf->fh, &sig, sizeof(uint64_t));
  if(size != sizeof(uint64_t)) return 0;
  if (sig != CDSIG) wzfatal("Unmatched signature. File corrupted.\n");
  bgzf_read(cf->fh, &(c->fmt), sizeof(char));
  bgzf_read(cf->fh, &(c->n), sizeof(uint64_t));
  c->compressed = 1;
  c->s = realloc(c->s, cdata_nbytes(c));
  bgzf_read(cf->fh, c->s, cdata_nbytes(c));
  cf->n++;
  return 1;
}

const uint8_t CX_BGZF_EOF[28] =
  "\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\033\0\3\0\0\0\0\0\0\0\0";

int cx_record_at(BGZF *fh, int64_t voffset) {
  if (voffset < 0 || (voffset & 0xFFFF) != 0) return 0;
  if (bgzf_seek(fh, voffset, SEEK_SET) < 0) return 0;
  /* after a seek into a concatenated file the block has to be pulled in
   * explicitly before the first read, same as read_cdata2() above */
  if (fh->block_length == 0) bgzf_read_block(fh);
  uint64_t sig = 0;
  if (bgzf_read(fh, &sig, sizeof(sig)) != (ssize_t)sizeof(sig)) return 0;
  return sig == CDSIG;
}

/**
 * Pull [*beg,*end) in to the record's own blocks, dropping the empty BGZF
 * members that sit at either end.
 *
 * A concatenated store is files joined end to end, and each file finishes with
 * the 28-byte empty member every BGZF writer emits. `yame index` records a
 * record's offset as wherever the reader stood after the previous record,
 * which is that trailing empty member -- so a record's byte range normally
 * OPENS with one, and the store's last record also CLOSES with the store's own.
 *
 * Copying those through looks harmless and is not. read_cdata2() tolerates a
 * single leading empty member (its bgzf_read_block() pre-call is exactly that
 * allowance), but two in a row exhaust it: bgzf_read() comes up short,
 * read_cdata2() reports end of stream, and every record after that point is
 * silently invisible. Two in a row is what you get by placing the store's last
 * record anywhere but last in the output -- its trailing empty member meets
 * the next record's leading one.
 */
void cx_trim_empty_members(FILE *in, int64_t *beg, int64_t *end, int at_eof) {
  uint8_t m[sizeof(CX_BGZF_EOF)];
  while (*end - *beg >= (int64_t) sizeof(m)) {   /* leading */
    if (fseeko(in, *beg, SEEK_SET) != 0) return;
    if (fread(m, 1, sizeof(m), in) != sizeof(m)) return;
    if (memcmp(m, CX_BGZF_EOF, sizeof(m)) != 0) break;
    *beg += sizeof(m);
  }
  /* Only the extent that runs to end of file can close with an empty member.
   * Every other extent stops at the next record's offset, and that offset
   * points AT the empty member, so it lies outside the range. Checking anyway
   * would seek to the far end of a multi-megabyte record and back for nothing,
   * throwing away the readahead the copy is about to use. */
  if (!at_eof) return;
  while (*end - *beg >= (int64_t) sizeof(m)) {   /* trailing */
    if (fseeko(in, *end - (int64_t) sizeof(m), SEEK_SET) != 0) return;
    if (fread(m, 1, sizeof(m), in) != sizeof(m)) return;
    if (memcmp(m, CX_BGZF_EOF, sizeof(m)) != 0) break;
    *end -= sizeof(m);
  }
}

/**
 * Is there an empty BGZF member sitting at this offset?
 *
 * A count of what was *written* cannot catch the failure this guards against:
 * the raw path always emits one record per requested name, and the records
 * that went missing were present in the bytes -- a reader just could not walk
 * to them, because two adjacent empty members stopped it early. So the check
 * has to be about reachability.
 *
 * It only has to look at record boundaries, though. An empty member never
 * occurs inside a record's own blocks -- a writer emits one exactly once, at
 * the end of a file -- so the only place a stray one can land is where two
 * copied extents meet, and those offsets are already known: they are the ones
 * being recorded for the output index. That makes the check O(records) rather
 * than O(blocks), which is the reason to prefer it; on NFS the two were not
 * far enough apart to separate by timing, so the argument is the cost model,
 * not a measurement.
 */
int cx_empty_member_at(const char *fname, int64_t off) {
  FILE *f = fopen(fname, "rb");
  if (!f) return 0;
  uint8_t m[sizeof(CX_BGZF_EOF)];
  int hit = 0;
  if (fseeko(f, off, SEEK_SET) == 0 &&
      fread(m, 1, sizeof(m), f) == sizeof(m) &&
      memcmp(m, CX_BGZF_EOF, sizeof(m)) == 0)
    hit = 1;
  fclose(f);
  return hit;
}

int cx_copy_bytes(FILE *in, int64_t beg, int64_t end, FILE *out) {
  if (fseeko(in, beg, SEEK_SET) != 0) return 0;
  size_t bufsz = 1<<20;
  uint8_t *buf = malloc(bufsz);
  if (!buf) return 0;
  int ok = 1;
  for (int64_t left = end - beg; left > 0; ) {
    size_t want = left < (int64_t)bufsz ? (size_t)left : bufsz;
    size_t got = fread(buf, 1, want, in);
    if (got != want || fwrite(buf, 1, got, out) != got) { ok = 0; break; }
    left -= got;
  }
  free(buf);
  return ok;
}

cfile_t open_cfile(char *fname) { /* for read */
  cfile_t cf = {0};
  if (strcmp(fname, "-")==0) {
    cf.fh = bgzf_dopen(fileno(stdin), "r");
  } else {
    cf.fh = bgzf_open(fname, "r");
  }
  if (cf.fh == NULL) {
    fprintf(stderr, "Error opening file %s\n", fname);
    exit(1);
  }
  cf.n = 0;
  return cf;
}

cdata_t read_cdata1(cfile_t *cf) {
  cdata_t c = {0};
  if (!read_cdata2(cf, &c)) return c;
  return c;
}

cdata_v* read_cdata(cfile_t *cf, int64_t beg, int64_t end) {

  if (beg < 0) beg = 0;
  if (end >= 0 && end < beg) wzfatal("End is smaller than beg");

  cdata_v *cs = init_cdata_v(10);
  cdata_t c = {0};
  int64_t i=0;
  for (i=0; end<0 || i<=end; ++i) {
    read_cdata2(cf, &c);
    if (i<beg) continue;
    if (c.n>0) {
      (*next_ref_cdata_v(cs)) = c;
      c.s = NULL;
    } else {
      break;
    }
  }
  return cs;
}

/* this function is memory intensive if there are many samples */
cdata_v* read_cdata_all(cfile_t *cf) {

  cdata_v *cs = init_cdata_v(10);
  while (1) {
    cdata_t c = read_cdata1(cf);
    if (c.n >0) push_cdata_v(cs, c);
    else break;
  }
  return cs;
}

cdata_v* read_cdata_from_head(cfile_t *cf, int64_t n) {
  cdata_v *cs = init_cdata_v(10);
  cdata_t c = {0};
  for (int64_t i=0; i<n; ++i) {
    read_cdata2(cf, &c);
    if (c.n <= 0) break;
    (*next_ref_cdata_v(cs)) = c;
    c.s = NULL;
  }
  return cs;
}

cdata_v* read_cdata_from_tail(cfile_t *cf, index_t *idx, int64_t n) {
  int npairs = 0;
  index_pair_t *pairs = index_pairs(idx, &npairs);
  if (n > npairs) n = npairs;
  int64_t *indices = malloc(n*sizeof(int64_t));
  for (int64_t i=npairs-n; i<npairs; ++i) {
    indices[i-npairs+n] = pairs[i].value;
  }
  cdata_v *cs = read_cdata_with_indices(cf, indices, n);
  free(indices);
  return cs;
}

cdata_v* read_cdata_with_indices(cfile_t *cf, const int64_t* indices, int n) {
  cdata_v *cs = init_cdata_v(10);
  cdata_t c = {0};

  for (int i = 0; i < n; i++) {
    int64_t index = indices[i];
    if (index < 0) {
      fprintf(stderr, "\n");
      fprintf(stderr, "[%s:%d] Index is negative.\n", __func__, __LINE__);
      fflush(stderr);
      exit(1);
    }

    // Reposition the file pointer using bgzf_seek
    if (bgzf_seek(cf->fh, index, SEEK_SET) != 0) {
      fprintf(stderr, "[%s:%d] Cannot seek input.\n", __func__, __LINE__);
      fflush(stderr);
      exit(1);
    }
    read_cdata2(cf, &c);
    if (c.n > 0) {
      (*next_ref_cdata_v(cs)) = c;
      c.s = NULL;
    } else {
      break;
    }
  }

  return cs;
}

cdata_v* read_cdata_with_snames(cfile_t *cf, index_t *idx, snames_t *snames) {
  // check if we have all sample names in index
  int64_t* indices = malloc(snames->n * sizeof(int64_t));
  for (int i = 0; i < snames->n; i++) {
    indices[i] = getIndex(idx, snames->s[i]);
    if (indices[i] == -1) {
      fprintf(stderr, "Cannot find sample %s in index.\n", snames->s[i]);
      fflush(stderr);
      exit(1);
    }
  }
  cdata_v *cs = read_cdata_with_indices(cf, indices, snames->n);
  free(indices);
  return cs;
}

void cdata_write1(BGZF *fp, cdata_t *c) {
  // Write the signature
  uint64_t sig = CDSIG;
  if (bgzf_write(fp, &sig, sizeof(uint64_t)) < 0) {
    fprintf(stderr, "Error writing signature to file\n");
    return;
  }

  // Write the format
  if (bgzf_write(fp, &(c->fmt), sizeof(uint8_t)) < 0) {
    fprintf(stderr, "Error writing format to file\n");
    return;
  }

  // Write the count
  if (bgzf_write(fp, &(c->n), sizeof(uint64_t)) < 0) {
    fprintf(stderr, "Error writing count to file\n");
    return;
  }

  // Write the data
  if (bgzf_write(fp, c->s, cdata_nbytes(c)) < 0) {
    fprintf(stderr, "Error writing data to file\n");
    return;
  }
}

void cdata_write(char *fname_out, cdata_t *c, const char *mode, int verbose) {

  if (!c->compressed) cdata_compress(c);  
  BGZF* fp;
  if (fname_out) fp = bgzf_open2(fname_out, mode);
  else fp = bgzf_dopen(fileno(stdout), mode);
  
  if (fp == NULL) {
    fprintf(stderr, "Error opening file for writing: %s\n", fname_out);
    return;
  }
  cdata_write1(fp, c);
  bgzf_close(fp);

  if (verbose) {
    fprintf(stderr, "[%s:%d] Stored as Format %c\n", __func__, __LINE__, c->fmt);
    fflush(stderr);
  }
}
