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
  c->s = realloc(c->s, cdata_nbytes(c));
  bgzf_read(cf->fh, c->s, cdata_nbytes(c));
  c->compressed = 1;
  cf->n++;
  return 1;
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
