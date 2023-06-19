#include "cgfile.h"

int read_cg2(cgfile_t *cgf, cgdata_t *cg) {
  cg->n = 0;
  uint64_t sig;
  int64_t size;
  if (cgf->fh->block_length == 0) bgzf_read_block(cgf->fh); /* somehow this is needed for concat'ed bgzipped files */
  size = bgzf_read(cgf->fh, &sig, sizeof(uint64_t));
  if(size != sizeof(uint64_t)) return 0;
  if (sig != CGSIG) wzfatal("Unmatched signature. File corrupted.\n");
  bgzf_read(cgf->fh, &(cg->fmt), sizeof(char));
  bgzf_read(cgf->fh, &(cg->n), sizeof(uint64_t));
  cg->s = realloc(cg->s, cgdata_nbytes(cg));
  bgzf_read(cgf->fh, cg->s, cgdata_nbytes(cg));
  cg->compressed = 1;
  cgf->n++;
  return 1;
}

cgfile_t open_cgfile(char *fname) { /* for read */
  cgfile_t cgf = {0};
  if (strcmp(fname, "-")==0) {
    cgf.fh = bgzf_dopen(fileno(stdin), "r");
  } else {
    cgf.fh = bgzf_open(fname, "r");
  }
  if (cgf.fh == NULL) {
    fprintf(stderr, "Error opening file %s\n", fname);
    exit(1);
  }
  cgf.n = 0;
  return cgf;
}

cgdata_t read_cg(cgfile_t *cgf) {
  cgdata_t cg = {0};
  if (!read_cg2(cgf, &cg)) return cg;
  return cg;
}

cgdata_v* read_cgs(cgfile_t *cgf, int64_t beg, int64_t end) {

  if (beg < 0) beg = 0;
  if (end >= 0 && end < beg) wzfatal("End is smaller than beg");

  cgdata_v *cgs = init_cgdata_v(10);
  cgdata_t cg = {0};
  int64_t i=0;
  for (i=0; end<0 || i<=end; ++i) {
    read_cg2(cgf, &cg);
    if (i<beg) continue;
    if (cg.n>0) {
      (*next_ref_cgdata_v(cgs)) = cg;
      cg.s = NULL;
    } else {
      break;
    }
  }
  return cgs;
}

/* this function is memory intensive if there are many samples */
cgdata_v* read_cgs_all(cgfile_t *cgf) {

  cgdata_v *cgs = init_cgdata_v(10);
  while (1) {
    cgdata_t cg = read_cg(cgf);
    if (cg.n >0) push_cgdata_v(cgs, cg);
    else break;
  }
  return cgs;
}

cgdata_v* read_cgs_from_head(cgfile_t *cgf, int64_t n) {
  cgdata_v *cgs = init_cgdata_v(10);
  cgdata_t cg = {0};
  for (int64_t i=0; i<n; ++i) {
    read_cg2(cgf, &cg);
    if (cg.n <= 0) break;
    (*next_ref_cgdata_v(cgs)) = cg;
    cg.s = NULL;
  }
  return cgs;
}

cgdata_v* read_cgs_from_tail(cgfile_t *cgf, index_t *idx, int64_t n) {
  int npairs = 0;
  index_pair_t *pairs = index_pairs(idx, &npairs);
  if (n > npairs) n = npairs;
  int64_t *indices = malloc(n*sizeof(int64_t));
  for (int64_t i=npairs-n; i<npairs; ++i) {
    indices[i-npairs+n] = pairs[i].value;
  }
  cgdata_v *cgs = read_cgs_with_indices(cgf, indices, n);
  free(indices);
  return cgs;
}

cgdata_v* read_cgs_with_indices(cgfile_t *cgf, const int64_t* indices, int n) {
  cgdata_v *cgs = init_cgdata_v(10);
  cgdata_t cg = {0};

  for (int i = 0; i < n; i++) {
    int64_t index = indices[i];
    assert(index >= 0);

    // Reposition the file pointer using bgzf_seek
    assert(bgzf_seek(cgf->fh, index, SEEK_SET) == 0);
    read_cg2(cgf, &cg);
    if (cg.n > 0) {
      (*next_ref_cgdata_v(cgs)) = cg;
      cg.s = NULL;
    } else {
      break;
    }
  }

  return cgs;
}

cgdata_v* read_cgs_with_snames(cgfile_t *cgf, index_t *idx, snames_t *snames) {
  // check if we have all sample names in index
  int64_t* indices = malloc(snames->n * sizeof(int64_t));
  for (int i = 0; i < snames->n; i++) {
    indices[i] = getIndex(idx, snames->array[i]);
    if (indices[i] == -1) {
      fprintf(stderr, "Cannot find sample %s in index.\n", snames->array[i]);
      fflush(stderr);
      exit(1);
    }
  }
  cgdata_v *cgs = read_cgs_with_indices(cgf, indices, snames->n);
  free(indices);
  return cgs;
}

void cgdata_write1(BGZF *fp, cgdata_t *cg) {
  // Write the signature
  uint64_t sig = CGSIG;
  if (bgzf_write(fp, &sig, sizeof(uint64_t)) < 0) {
    fprintf(stderr, "Error writing signature to file\n");
    bgzf_close(fp);
    return;
  }

  // Write the format
  if (bgzf_write(fp, &(cg->fmt), sizeof(uint8_t)) < 0) {
    fprintf(stderr, "Error writing format to file\n");
    bgzf_close(fp);
    return;
  }

  // Write the count
  if (bgzf_write(fp, &(cg->n), sizeof(uint64_t)) < 0) {
    fprintf(stderr, "Error writing count to file\n");
    bgzf_close(fp);
    return;
  }

  // Write the data
  if (bgzf_write(fp, cg->s, cgdata_nbytes(cg)) < 0) {
    fprintf(stderr, "Error writing data to file\n");
    bgzf_close(fp);
    return;
  }
}

void cgdata_write(char *fname_out, cgdata_t *cg, const char *mode, int verbose) {

  if (!cg->compressed) recompress(cg);
  
  BGZF* fp;
  if (fname_out) fp = bgzf_open2(fname_out, mode);
  else fp = bgzf_dopen(fileno(stdout), mode);
    
  if (fp == NULL) {
    fprintf(stderr, "Error opening file for writing: %s\n", fname_out);
    return;
  }
  cgdata_write1(fp, cg);
  bgzf_close(fp);

  if (verbose) {
    fprintf(stderr, "[%s:%d] Stored as Format %c\n", __func__, __LINE__, cg->fmt);
    fflush(stderr);
  }
}
