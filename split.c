#include "cgfile.h"

static int usage() {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: kycg split [options] <in.cg> out_prefix\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -v        verbose\n");
  fprintf(stderr, "    -s        sample name list\n");
  fprintf(stderr, "    -h        This help\n");
  fprintf(stderr, "\n");

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

  if (optind + 2 > argc) { 
    usage(); 
    wzfatal("Please supply input file.\n");
  }

  cgfile_t cgf = open_cgfile(argv[optind++]);
  char *prefix = argv[optind];

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
  
  int i = 0;
  for (i=0; ; ++i) {
    cgdata_t cg = read_cg(&cgf);
    if (cg.n == 0) break;
    char *tmp = NULL;
    if (snames_n) {
      tmp = malloc(strlen(prefix)+strlen(snames[i])+1000);
      sprintf(tmp, "%s%s.cg", prefix, snames[i]);
    } else {
      tmp = malloc(strlen(prefix) + 1000);
      sprintf(tmp, "%s_split_%i.cg", prefix, i+1);
    }
    cgdata_write(tmp, &cg, "wb", verbose);
    free(tmp);
    free(cg.s);
  }
  
  return 0;
}
