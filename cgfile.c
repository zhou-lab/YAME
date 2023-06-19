#include "cgfile.h"

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
  if (!read_cg_(cgf, &cg)) return cg;
  return cg;
}

cgdata_v* read_cgs(cgfile_t *cgf, int64_t beg, int64_t end) {

  if (beg < 0) beg = 0;
  if (end >= 0 && end < beg) wzfatal("End is smaller than beg");

  cgdata_v *cgs = init_cgdata_v(10);
  cgdata_t cg = {0};
  int64_t i=0;
  for (i=0; end<0 || i<=end; ++i) {
    read_cg_(cgf, &cg);
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
  for (int64_t i=0; i<=n; ++i) {
    read_cg_(cgf, &cg);
    if (cg.n <= 0) break;
    (*next_ref_cgdata_v(cgs)) = cg;
    cg.s = NULL;
  }
  return cgs;
}

cgdata_v* read_cgs_from_tail(cgfile_t *cgf, index_t *idx, int64_t n) {
  int64_t last = kh_size(idx);
  if (n > last) n = last;
  int64_t *indices = malloc(n*sizeof(int64_t));
  index_pair_t *pairs = index_pairs(idx);
  for (int64_t i=last-n+1; i<=last; ++i) {
    indices[i] = pairs[i].value;
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
    read_cg_(cgf, &cg);
    if (cg.n > 0) {
      (*next_ref_cgdata_v(cgs)) = cg;
      cg.s = NULL;
    } else {
      break;
    }
  }

  return cgs;
}

cgdata_v* read_cgs_with_snames(cgfile_t *cgf, snames_t *snames) {
  // check if we have all sample names in index
  int64_t* indices = malloc(snames->n * sizeof(int64_t));
  for (unsigned i = 0; i < snames->n; i++) {
    indices[i] = getIndex(idx, snames->array[i]);
    if (indices[i] == -1) {
      fprintf(stderr, "Cannot find sample %s in index.\n", snames->array[i]);
      fflush(stderr);
      exit(1);
    }
  }
  cgdata_v *cgs = read_cgs_with_indices(&cgf, indices, sample_names->size);
  free(indices);
  return cgs;
}
