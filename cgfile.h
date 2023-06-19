typedef struct cgfile_t {
  BGZF *fh;
  int n;                        /* number of samples read */
} cgfile_t;

static inline int read_cg_(cgfile_t *cgf, cgdata_t *cg) {
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
