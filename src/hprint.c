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

#include <string.h>
#include "cfile.h"

/* ANSI color codes */
#define ANSI_RESET   "\x1b[0m"
#define ANSI_METH    "\x1b[33m"   /* yellow - methylated   */
#define ANSI_UNMETH  "\x1b[34m"   /* blue   - unmethylated */
#define ANSI_NA      "\x1b[90m"   /* grey   - NA           */
#define ANSI_HIGH    "\x1b[31m"   /* red    - high beta    */
#define ANSI_MED     "\x1b[33m"   /* yellow - medium beta  */
#define ANSI_LOW     "\x1b[34m"   /* blue   - low beta     */

/* UTF-8 block characters for fmt6 */
#define CH_METH   "\xe2\x96\x88"  /* █ U+2588 FULL BLOCK  */
#define CH_UNMETH "\xe2\x96\x91"  /* ░ U+2591 LIGHT SHADE */
#define CH_NA     "."

typedef struct {
  char    *chrm;
  uint64_t first_row; /* 1-based */
  uint64_t last_row;  /* 1-based */
  uint64_t n_cpgs;
  uint64_t n_cols;    /* output columns after windowing */
} chrom_info_t;

static int usage(void) {
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: yame hprint [options] <in.cx>\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Modes:\n");
  fprintf(stderr, "  -R <ref>          Whole-genome view: one column per CpG window across all chroms.\n");
  fprintf(stderr, "  -R <ref> -r <reg> Region view: rows=samples, columns=CpG sites in region.\n");
  fprintf(stderr, "  (neither)         Legacy full-dataset dump (fmt6 only).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "  -c            Disable ANSI color output (default: color on)\n");
  fprintf(stderr, "  -g            Granular output: 0-9 deciles instead of H/M/L\n");
  fprintf(stderr, "  -R <ref.cr>   Reference coordinate file (format 7)\n");
  fprintf(stderr, "  -r <region>   Genomic region: chr16  or  chr16:10000000-10100000\n");
  fprintf(stderr, "  -l <int>      Sample label column width (default: 20)\n");
  fprintf(stderr, "  -t <int>      Ruler tick interval in columns (default: 10)\n");
  fprintf(stderr, "  -w <int>      Max data columns; wider views are window-averaged (default: 80)\n");
  fprintf(stderr, "  -h            This help\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Symbols (per-site):  fmt6: " CH_METH " meth  " CH_UNMETH " unmeth  . NA   fmt3/4: H/M/L/.   fmt0: 1/0\n");
  fprintf(stderr, "Symbols (windowed):  H >0.67   M 0.33-0.67   L <0.33   . no coverage\n");
  fprintf(stderr, "\n");
  return 1;
}

/* ------------------------------------------------------------------ */
/* Helpers shared by region and whole-genome modes                     */
/* ------------------------------------------------------------------ */

/* Parse "chr16:1000000-2000000" or bare "chr16" → chrm (caller frees), beg1, end1. */
static int parse_region(const char *reg, char **chrm_out,
                         uint64_t *beg1_out, uint64_t *end1_out) {
  const char *colon = strchr(reg, ':');
  if (!colon) {
    *chrm_out = strdup(reg);
    *beg1_out = 0;
    *end1_out = UINT64_MAX;
    return 0;
  }
  const char *dash = strchr(colon + 1, '-');
  if (!dash) return -1;
  *chrm_out  = strndup(reg, (size_t)(colon - reg));
  *beg1_out  = strtoull(colon + 1, NULL, 10);
  *end1_out  = strtoull(dash  + 1, NULL, 10);
  return (*beg1_out > 0 && *end1_out >= *beg1_out) ? 0 : -1;
}

/*
 * Scan format 7 to find row range and CpG count for a region.
 */
static void get_region_info(cdata_t *cr, const char *chrm,
                             uint64_t beg1, uint64_t end1,
                             uint64_t *n_out,
                             uint64_t *first_row, uint64_t *last_row,
                             uint64_t *last_cpg) {
  row_reader_t rdr = {0};
  uint64_t n = 0;
  *first_row = *last_row = 0;

  while (row_reader_next_loc(&rdr, cr)) {
    if (strcmp(rdr.chrm, chrm) != 0) {
      if (n > 0) break;
      continue;
    }
    if (rdr.value < beg1) continue;
    if (rdr.value > end1) break;
    if (n == 0) *first_row = rdr.index;
    n++;
    *last_row = rdr.index;
    if (last_cpg) *last_cpg = rdr.value;
  }
  *n_out = n;
}

/*
 * Second pass to grab window anchor coordinates.
 */
static void get_win_pos(cdata_t *cr, const char *chrm,
                         uint64_t beg1, uint64_t end1,
                         uint64_t win_size, uint64_t n_cols,
                         uint64_t *win_pos) {
  row_reader_t rdr = {0};
  uint64_t n = 0;
  uint64_t col = 0;

  while (row_reader_next_loc(&rdr, cr)) {
    if (strcmp(rdr.chrm, chrm) != 0) {
      if (n > 0) break;
      continue;
    }
    if (rdr.value < beg1) continue;
    if (rdr.value > end1) break;

    if (n == col * win_size) {
      win_pos[col++] = rdr.value;
      if (col == n_cols) break;
    }
    n++;
  }
}

/*
 * Walk the entire compressed fmt7 reference and record, per chromosome,
 * the 1-based first/last row index and CpG count.
 */
static chrom_info_t *collect_genome_chroms(cdata_t *cr, int *n_out) {
  row_reader_t  rdr = {0};
  chrom_info_t *ch  = NULL;
  int           n   = 0;

  while (row_reader_next_loc(&rdr, cr)) {
    if (n == 0 || strcmp(rdr.chrm, ch[n - 1].chrm) != 0) {
      ch = realloc(ch, (size_t)(n + 1) * sizeof(chrom_info_t));
      ch[n].chrm      = strdup(rdr.chrm);
      ch[n].first_row = rdr.index;
      ch[n].last_row  = rdr.index;
      ch[n].n_cpgs    = 1;
      ch[n].n_cols    = 0;
      n++;
    } else {
      ch[n - 1].last_row = rdr.index;
      ch[n - 1].n_cpgs++;
    }
  }
  *n_out = n;
  return ch;
}

/*
 * Single-line genome ruler: '|' at each chromosome boundary followed by
 * as many characters of the chromosome name as fit, '.' filling the rest.
 */
static void print_genome_ruler(chrom_info_t *ch, int n_chroms,
                                uint64_t total_cols, int label_w) {
  char *buf = malloc(total_cols + 1);
  memset(buf, '.', total_cols);
  buf[total_cols] = '\0';

  uint64_t col = 0;
  for (int i = 0; i < n_chroms; i++) {
    buf[col] = '|';
    if (ch[i].n_cols > 1) {
      const char *name = ch[i].chrm;
      if (strncmp(name, "chr", 3) == 0) name += 3;
      size_t avail = ch[i].n_cols - 1;
      size_t nlen  = strlen(name);
      memcpy(buf + col + 1, name, nlen < avail ? nlen : avail);
    }
    col += ch[i].n_cols;
  }
  fprintf(stdout, "%*s  %s\n", label_w, "", buf);
  free(buf);
}

/*
 * Three-line ruler.
 * pos[n_cols] holds the genomic coordinate anchoring each output column
 * (in windowed mode, that is the first CpG of each window).
 * n_cpgs is the true CpG count; win_size > 1 triggers a "(win=N)" annotation.
 */
static void print_ruler(const char *chrm, uint64_t *pos, uint64_t n_cols,
                         uint64_t n_cpgs, uint64_t win_size, uint64_t last_cpg,
                         int label_w, int tick_every) {
  if (win_size > 1)
    fprintf(stdout, "%*s  %s:%"PRIu64"-%"PRIu64" (%"PRIu64" CpGs, win=%"PRIu64")\n",
            label_w, "", chrm, pos[0], last_cpg, n_cpgs, win_size);
  else
    fprintf(stdout, "%*s  %s:%"PRIu64"-%"PRIu64" (%"PRIu64" CpGs)\n",
            label_w, "", chrm, pos[0], last_cpg, n_cpgs);

  fprintf(stdout, "%*s  ", label_w, "");
  for (uint64_t i = 0; i < n_cols; ++i)
    fputc(i % (uint64_t)tick_every == 0 ? '|' : '.', stdout);
  fputc('\n', stdout);

  char *buf = malloc(n_cols + 1);
  memset(buf, ' ', n_cols);
  buf[n_cols] = '\0';
  for (uint64_t i = 0; i < n_cols; i += (uint64_t)tick_every) {
    char num[24];
    int len = snprintf(num, sizeof(num), "%"PRIu64, pos[i]);
    size_t avail = n_cols - i;
    memcpy(buf + i, num, (size_t)len < avail ? (size_t)len : avail);
  }
  fprintf(stdout, "%*s  %s\n", label_w, "", buf);
  free(buf);
}

/* Emit a single H/M/L character or 0-9 decile (with optional color) for a beta in [0,1]. */
static void print_beta(double b, int color, int granular) {
  if (granular) {
    int decile = (int)(b * 10.0);
    if (decile > 9) decile = 9;
    if (decile < 0) decile = 0;

    if (color) {
      /* Blue-White-Red spectrum using ANSI 256 colors */
      static const char *colors[] = {
        "\x1b[38;5;21m",  /* 0: deep blue */
        "\x1b[38;5;27m",  /* 1 */
        "\x1b[38;5;33m",  /* 2 */
        "\x1b[38;5;39m",  /* 3 */
        "\x1b[38;5;45m",  /* 4: light blue */
        "\x1b[38;5;231m", /* 5: white */
        "\x1b[38;5;224m", /* 6: light red/pink */
        "\x1b[38;5;217m", /* 7 */
        "\x1b[38;5;210m", /* 8 */
        "\x1b[38;5;196m"  /* 9: deep red */
      };
      fputs(colors[decile], stdout);
    }
    fputc('0' + decile, stdout);
    if (color) fputs(ANSI_RESET, stdout);
  } else {
    char ch         = b > 0.67 ? 'H' : b > 0.33 ? 'M' : 'L';
    const char *col = b > 0.67 ? ANSI_HIGH : b > 0.33 ? ANSI_MED : ANSI_LOW;
    if (color) fputs(col, stdout);
    fputc(ch, stdout);
    if (color) fputs(ANSI_RESET, stdout);
  }
}

/* Decode one fmt3 record at byte *ip in the compressed stream s[].
 * Advances *ip. Returns number of rows covered; sets *M and *U.
 * Zero-runs have M=U=0 and nrows > 1; non-zero records have nrows = 1. */
static inline uint64_t fmt3_next(const uint8_t *s, uint64_t *ip,
                                  uint64_t *M, uint64_t *U) {
  uint8_t type = s[*ip] & 0x3;
  if (type == 0) {
    uint16_t v; memcpy(&v, s + *ip, 2);
    *ip += 2; *M = 0; *U = 0;
    return (uint64_t)(v >> 2);
  }
  if (type == 1) {
    *M = s[*ip] >> 5; *U = (s[*ip] >> 2) & 0x7; *ip += 1;
  } else if (type == 2) {
    uint16_t v; memcpy(&v, s + *ip, 2);
    uint64_t mv = (uint64_t)(v >> 2); *U = mv & 0x7f; *M = mv >> 7; *ip += 2;
  } else {
    uint64_t v; memcpy(&v, s + *ip, 8);
    uint64_t mv = v >> 2; *U = mv & 0x7fffffff; *M = mv >> 31; *ip += 8;
  }
  return 1;
}

/* Single-pass fmt3 region printer.
 * Walks the compressed stream once; accumulates windowed averages and emits
 * H/M/L/. directly — zero intermediate allocation. */
static void stream_fmt3_region(const cdata_t *c,
                                uint64_t first_row, uint64_t last_row,
                                uint64_t win_size, uint64_t n_cols, int color, int granular) {
  uint64_t i = 0, row = 0, pM = 0, pU = 0, prem = 0;
  double   wsum = 0.0;
  uint64_t wvalid = 0, wpos = 0, wemit = 0;

  while (wemit < n_cols && (prem > 0 || i < c->n)) {
    if (prem == 0) { if (i >= c->n) break; prem = fmt3_next(c->s, &i, &pM, &pU); }
    if (row < first_row) {
      uint64_t skip = first_row - row;
      if (skip >= prem) { row += prem; prem = 0; continue; }
      row += skip; prem -= skip;
    }
    if (row > last_row) break;

    uint64_t avail = last_row - row + 1;
    uint64_t take  = prem < avail ? prem : avail;
    int      cov   = (pM + pU > 0);
    double   beta  = cov ? (double)pM / (pM + pU) : 0.0;

    uint64_t cnt = take;
    while (cnt > 0) {
      uint64_t wleft = win_size - wpos, t = cnt < wleft ? cnt : wleft;
      if (cov) { wsum += (double)t * beta; wvalid += t; }
      wpos += t; cnt -= t;
      if (wpos == win_size) {
        if (wvalid == 0) { if (color) fputs(ANSI_NA, stdout); fputc('.', stdout); if (color) fputs(ANSI_RESET, stdout); }
        else             print_beta(wsum / wvalid, color, granular);
        wemit++; wsum = 0.0; wvalid = 0; wpos = 0;
      }
    }
    row += take; prem -= take;
  }
  if (wpos > 0 && wemit < n_cols) {
    if (wvalid == 0) { if (color) fputs(ANSI_NA, stdout); fputc('.', stdout); if (color) fputs(ANSI_RESET, stdout); }
    else             print_beta(wsum / wvalid, color, granular);
  }
}

/* Single-pass fmt3 genome printer.
 * One forward scan covers all chromosomes in order; emits windowed H/M/L/.
 * per chromosome — zero intermediate allocation. */
static void stream_fmt3_genome(const cdata_t *c,
                                const chrom_info_t *ch, int n_chroms,
                                uint64_t win_size, int color, int granular) {
  uint64_t i = 0, row = 0, pM = 0, pU = 0, prem = 0;

  for (int ci = 0; ci < n_chroms; ci++) {
    uint64_t cstart = ch[ci].first_row - 1, cend = ch[ci].last_row - 1;
    double   wsum = 0.0;
    uint64_t wvalid = 0, wpos = 0, wemit = 0;

    while (row <= cend && (prem > 0 || i < c->n)) {
      if (prem == 0) { if (i >= c->n) break; prem = fmt3_next(c->s, &i, &pM, &pU); }
      if (row < cstart) {
        uint64_t skip = cstart - row;
        if (skip >= prem) { row += prem; prem = 0; continue; }
        row += skip; prem -= skip;
      }
      if (row > cend) break;

      uint64_t avail = cend - row + 1;
      uint64_t take  = prem < avail ? prem : avail;
      int      cov   = (pM + pU > 0);
      double   beta  = cov ? (double)pM / (pM + pU) : 0.0;

      uint64_t cnt = take;
      while (cnt > 0) {
        uint64_t wleft = win_size - wpos, t = cnt < wleft ? cnt : wleft;
        if (cov) { wsum += (double)t * beta; wvalid += t; }
        wpos += t; cnt -= t;
        if (wpos == win_size) {
          if (wemit < ch[ci].n_cols) {
            if (wvalid == 0) { if (color) fputs(ANSI_NA, stdout); fputc('.', stdout); if (color) fputs(ANSI_RESET, stdout); }
            else             print_beta(wsum / wvalid, color, granular);
            wemit++;
          }
          wsum = 0.0; wvalid = 0; wpos = 0;
        }
      }
      row += take; prem -= take;
    }
    if (wpos > 0 && wemit < ch[ci].n_cols) {
      if (wvalid == 0) { if (color) fputs(ANSI_NA, stdout); fputc('.', stdout); if (color) fputs(ANSI_RESET, stdout); }
      else             print_beta(wsum / wvalid, color, granular);
    }
  }
}


/*
 * Print one sample row in windowed mode.
 * Each output column is the average beta over win_size consecutive sites.
 */
static void print_region_sample_windowed(cdata_t *c, uint64_t start_idx, uint64_t n_sites,
                                          uint64_t win_size, uint64_t n_cols, int color, int granular) {
  for (uint64_t j = 0; j < n_cols; ++j) {
    uint64_t i0 = start_idx + j * win_size;
    uint64_t i1 = i0 + win_size;
    if (i1 > start_idx + n_sites) i1 = start_idx + n_sites;

    double sum = 0.0;
    int    valid = 0;

    if (c->fmt == '6') {
      for (uint64_t i = i0; i < i1; ++i) {
        if (FMT6_IN_UNI(*c, i)) {
          sum += FMT6_IN_SET(*c, i) ? 1.0 : 0.0;
          valid++;
        }
      }
    } else if (c->fmt == '0') {
      for (uint64_t i = i0; i < i1; ++i) {
        sum += FMT0_IN_SET(*c, i) ? 1.0 : 0.0;
        valid++;
      }
    } else if (c->fmt == '4') {
      for (uint64_t i = i0; i < i1; ++i) {
        float v; memcpy(&v, c->s + i * c->unit, sizeof(float));
        if (v >= 0.0f) { sum += (double)v; valid++; }
      }
    } else {
      fputc('?', stdout);
      continue;
    }

    if (valid == 0) {
      if (color) fputs(ANSI_NA, stdout);
      fputc('.', stdout);
      if (color) fputs(ANSI_RESET, stdout);
    } else {
      print_beta(sum / valid, color, granular);
    }
  }
}

/* Print one sample row in region mode, dispatching by format. */
static void print_region_sample(cdata_t *c, uint64_t start_idx, uint64_t n_sites, int color, int granular) {
  for (uint64_t i = start_idx; i < start_idx + n_sites; ++i) {

    if (c->fmt == '6') {
      if (!FMT6_IN_UNI(*c, i)) {
        if (color) fputs(ANSI_NA,     stdout);
        fputs(CH_NA,     stdout);
        if (color) fputs(ANSI_RESET, stdout);
      } else if (granular) {
        print_beta(FMT6_IN_SET(*c, i) ? 1.0 : 0.0, color, 1);
      } else if (FMT6_IN_SET(*c, i)) {
        if (color) fputs(ANSI_METH,   stdout);
        fputs(CH_METH,   stdout);
        if (color) fputs(ANSI_RESET, stdout);
      } else {
        if (color) fputs(ANSI_UNMETH, stdout);
        fputs(CH_UNMETH, stdout);
        if (color) fputs(ANSI_RESET, stdout);
      }

    } else if (c->fmt == '0') {
      double b = FMT0_IN_SET(*c, i) ? 1.0 : 0.0;
      if (granular) {
        print_beta(b, color, 1);
      } else {
        if (color) fputs(b > 0.5 ? ANSI_METH : ANSI_UNMETH, stdout);
        fputc(b > 0.5 ? '1' : '0', stdout);
        if (color) fputs(ANSI_RESET, stdout);
      }

    } else if (c->fmt == '4') {
      float v;
      memcpy(&v, c->s + i * c->unit, sizeof(float));
      if (v < 0.0f) {
        if (color) fputs(ANSI_NA, stdout);
        fputs(CH_NA, stdout);
        if (color) fputs(ANSI_RESET, stdout);
      } else {
        print_beta((double)v, color, granular);
      }

    } else {
      fputc('?', stdout);
    }
  }
}



/* ------------------------------------------------------------------ */
/* main                                                                 */
/* ------------------------------------------------------------------ */

int main_hprint(int argc, char *argv[]) {
  int c;
  int color = 1, label_w = 20, tick_every = 10, max_cols = 80, granular = 0;
  char *fname_cr = NULL, *region = NULL;

  while ((c = getopt(argc, argv, "cgR:r:l:t:w:h")) >= 0) {
    switch (c) {
    case 'c': color      = 0;              break;
    case 'g': granular   = 1;              break;
    case 'R': fname_cr   = strdup(optarg); break;
    case 'r': region     = strdup(optarg); break;
    case 'l': label_w    = atoi(optarg);   break;
    case 't': tick_every = atoi(optarg);   break;
    case 'w': max_cols   = atoi(optarg);   break;
    case 'h': return usage();
    default:  usage(); wzfatal("Unrecognized option: %c.\n", c);
    }
  }
  if (max_cols < 1) max_cols = 1;

  if (optind + 1 > argc) { usage(); wzfatal("Please supply input file.\n"); }
  char *fname = argv[optind];

  /* ---- region mode ---- */
  if (region) {
    if (!fname_cr) wzfatal("-R <ref.cr> is required when -r is used.\n");

    char *chrm = NULL;
    uint64_t beg1 = 0, end1 = 0;
    if (parse_region(region, &chrm, &beg1, &end1) < 0)
      wzfatal("Cannot parse region '%s'. Expected 'chrN' or 'chrN:beg-end' (1-based).\n", region);

    cfile_t cf_cr = open_cfile(fname_cr);
    cdata_t cr = read_cdata1(&cf_cr);
    bgzf_close(cf_cr.fh);
    if (cr.fmt != '7')
      wzfatal("Reference must be format 7 (.cr), got '%c'.\n", cr.fmt);

    uint64_t cr_n = fmt7_data_length(&cr);

    uint64_t n_pos = 0, first_row = 0, last_row = 0, last_cpg = 0;
    get_region_info(&cr, chrm, beg1, end1, &n_pos, &first_row, &last_row, &last_cpg);

    if (n_pos == 0) {
      fprintf(stderr, "[hprint] No CpGs found in %s:%"PRIu64"-%"PRIu64"\n",
              chrm, beg1, end1);
      free_cdata(&cr); free(chrm); free(fname_cr); free(region);
      return 1;
    }

    /* preflight: check dimensions against the first sample before printing anything */
    {
      cfile_t cf_check = open_cfile(fname);
      cdata_t c_check = read_cdata1(&cf_check);
      uint64_t data_n = cdata_n(&c_check);
      free_cdata(&c_check);
      bgzf_close(cf_check.fh);
      if (data_n != cr_n)
        wzfatal("[hprint] Dimension mismatch: reference has %"PRIu64" CpGs "
                "but data in %s has %"PRIu64". "
                "Ensure the same reference was used to create both files.\n",
                cr_n, fname, data_n);
    }

    /* Compute window size so output fits within max_cols columns. */
    uint64_t win_size = 1;
    uint64_t n_cols   = n_pos;
    if (n_pos > (uint64_t)max_cols) {
      win_size = (n_pos + (uint64_t)max_cols - 1) / (uint64_t)max_cols;
      n_cols   = (n_pos + win_size - 1) / win_size;
    }

    /* Build per-column position anchors (first CpG of each window). */
    uint64_t *win_pos = malloc(n_cols * sizeof(uint64_t));
    get_win_pos(&cr, chrm, beg1, end1, win_size, n_cols, win_pos);
    free_cdata(&cr);

    snames_t snames = loadSampleNamesFromIndex(fname);

    print_ruler(chrm, win_pos, n_cols, n_pos, win_size, last_cpg, label_w, tick_every);
    free(win_pos);

    cfile_t cf = open_cfile(fname);
    int si = 0;
    for (;;) {
      cdata_t cin = read_cdata1(&cf);
      if (cin.n == 0) { free_cdata(&cin); break; }

      const char *label = (si < snames.n) ? snames.s[si] : "";
      fprintf(stdout, "%-*.*s  ", label_w, label_w, label);

      if (cin.fmt == '3') {
        stream_fmt3_region(&cin, first_row - 1, last_row - 1, win_size, n_cols, color, granular);
      } else {
        decompress_in_situ(&cin);
        if (last_row > cin.n)
          wzfatal("[hprint] Region rows %"PRIu64"-%"PRIu64" exceed data size "
                  "(%"PRIu64"). This should not happen if dimensions match.\n",
                  first_row, last_row, cin.n);
        if (win_size > 1)
          print_region_sample_windowed(&cin, first_row - 1, n_pos, win_size, n_cols, color, granular);
        else
          print_region_sample(&cin, first_row - 1, n_pos, color, granular);
      }
      fputc('\n', stdout);
      free_cdata(&cin);
      si++;
    }

    bgzf_close(cf.fh);
    cleanSampleNames2(snames);
    free(chrm); free(fname_cr); free(region);
    return 0;
  }

  /* ---- whole-genome mode: -R without -r ---- */
  if (fname_cr) {
    cfile_t cf_cr = open_cfile(fname_cr);
    cdata_t cr    = read_cdata1(&cf_cr);
    bgzf_close(cf_cr.fh);
    if (cr.fmt != '7')
      wzfatal("Reference must be format 7 (.cr), got '%c'.\n", cr.fmt);

    uint64_t cr_n = fmt7_data_length(&cr);

    int           n_chroms = 0;
    chrom_info_t *ch       = collect_genome_chroms(&cr, &n_chroms);
    free_cdata(&cr);

    if (n_chroms == 0)
      wzfatal("[hprint] Reference contains no entries.\n");

    /* preflight dimension check */
    {
      cfile_t  cf_check = open_cfile(fname);
      cdata_t  c_check  = read_cdata1(&cf_check);
      uint64_t data_n   = cdata_n(&c_check);
      free_cdata(&c_check);
      bgzf_close(cf_check.fh);
      if (data_n != cr_n)
        wzfatal("[hprint] Dimension mismatch: reference has %"PRIu64" CpGs "
                "but data in %s has %"PRIu64". "
                "Ensure the same reference was used to create both files.\n",
                cr_n, fname, data_n);
    }

    /* global win_size keeps chromosome boundaries intact */
    uint64_t total_cpgs = 0;
    for (int i = 0; i < n_chroms; i++) total_cpgs += ch[i].n_cpgs;
    uint64_t win_size = (total_cpgs + (uint64_t)max_cols - 1) / (uint64_t)max_cols;
    if (win_size < 1) win_size = 1;
    uint64_t total_cols = 0;
    for (int i = 0; i < n_chroms; i++) {
      ch[i].n_cols = (ch[i].n_cpgs + win_size - 1) / win_size;
      total_cols  += ch[i].n_cols;
    }

    snames_t snames = loadSampleNamesFromIndex(fname);
    print_genome_ruler(ch, n_chroms, total_cols, label_w);

    cfile_t cf = open_cfile(fname);
    int si = 0;
    for (;;) {
      cdata_t cin = read_cdata1(&cf);
      if (cin.n == 0) { free_cdata(&cin); break; }

      const char *label = (si < snames.n) ? snames.s[si] : "";
      fprintf(stdout, "%-*.*s  ", label_w, label_w, label);

      if (cin.fmt == '3') {
        stream_fmt3_genome(&cin, ch, n_chroms, win_size, color, granular);
      } else {
        decompress_in_situ(&cin);
        for (int ci = 0; ci < n_chroms; ci++) {
          if (ch[ci].last_row > cin.n)
            wzfatal("[hprint] Chromosome %s rows exceed data size. "
                    "Dimensions should have matched.\n", ch[ci].chrm);
          if (win_size > 1)
            print_region_sample_windowed(&cin, ch[ci].first_row - 1, ch[ci].n_cpgs, win_size, ch[ci].n_cols, color, granular);
          else
            print_region_sample(&cin, ch[ci].first_row - 1, ch[ci].n_cpgs, color, granular);
        }
      }
      fputc('\n', stdout);
      free_cdata(&cin);
      si++;
    }


    bgzf_close(cf.fh);
    cleanSampleNames2(snames);
    for (int i = 0; i < n_chroms; i++) free(ch[i].chrm);
    free(ch);
    free(fname_cr);
    return 0;
  }

  /* ---- original full-dataset mode (fmt6 only) ---- */
  cfile_t cf = open_cfile(fname);
  for (;;) {
    cdata_t c2 = read_cdata1(&cf);
    if (c2.n == 0) { free_cdata(&c2); break; }
    decompress_in_situ(&c2);

    if (c2.fmt != '6') {
      fprintf(stderr, "[hprint] Only format 6 supported in full-dataset mode "
              "(got '%c'). Use -r/-R for other formats.\n", c2.fmt);
      free_cdata(&c2);
      bgzf_close(cf.fh);
      return 1;
    }

    for (uint64_t i = 0; i < c2.n; ++i) {
      if (FMT6_IN_UNI(c2, i)) {
        if (FMT6_IN_SET(c2, i)) {
          if (color) fprintf(stdout, ANSI_METH   "1" ANSI_RESET);
          else fputc('1', stdout);
        } else {
          if (color) fprintf(stdout, ANSI_UNMETH "0" ANSI_RESET);
          else fputc('0', stdout);
        }
      } else {
        if (color) fprintf(stdout, ANSI_NA "2" ANSI_RESET);
        else fputc('2', stdout);
      }
    }
    fputc('\n', stdout);
    free_cdata(&c2);
  }

  bgzf_close(cf.fh);
  if (fname_cr) free(fname_cr);
  if (region)   free(region);
  return 0;
}
