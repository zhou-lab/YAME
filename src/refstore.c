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

/**
 * Which reference a file belongs to, from the one thing it always carries.
 *
 * A CX record knows how many rows it has and nothing about what they are, so
 * a command that needs coordinates has had to be told -- every time, for a
 * file whose row count already implies the answer. No two row spaces in the
 * catalogue share a count, so the count identifies the reference outright:
 * 29,401,795 rows is hg38 and can be nothing else.
 *
 * This only ever RESOLVES a reference; it never fetches one. A command that
 * finds the reference missing should say which `yame fetch` would supply it
 * and stop, rather than start a multi-megabyte download nobody asked for in
 * the middle of a print.
 */

#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <dirent.h>

#include "assets.h"
#include "cfile.h"
#include "registry.h"

int yame_ref_for_rows(uint64_t rows, const char *store_override,
                      const char *want_kind, char *path, size_t n,
                      const char **name, const char **fetch) {
  if (name) *name = NULL;
  if (fetch) *fetch = NULL;
  if (path && n) path[0] = '\0';
  if (!rows) return YAME_REF_UNKNOWN;

  const yame_ref_rows_t *e = NULL;
  for (size_t i = 0; i < YAME_REF_ROWS_N; ++i)
    if (YAME_REF_ROWS[i].rows == rows) { e = &YAME_REF_ROWS[i]; break; }

  if (!e) return YAME_REF_UNKNOWN;
  if (name) *name = e->name;
  if (fetch) *fetch = e->fetch;

  /* Identified, but not the kind of reference the caller can use: an array's
   * row space is a probe ordering, which is not a coordinate stream. Saying
   * so by name beats "not found", because the file is not the problem. */
  if (want_kind && strcmp(e->kind, want_kind) != 0) return YAME_REF_WRONG_KIND;

  char root[4096];
  yame_assets_root(store_override, NULL, root, sizeof(root));
  if (yame_assets_join(path, n, root, e->store_path) != 0)
    return YAME_REF_MISSING;

  return yame_assets_is_file(path) ? YAME_REF_OK : YAME_REF_MISSING;
}

void yame_ref_explain(FILE *out, uint64_t rows, int status, const char *name,
                      const char *fetch, const char *flag) {
  switch (status) {
  case YAME_REF_UNKNOWN:
    fprintf(out, "%s is required: %" PRIu64 " rows matches no reference this "
                 "build knows, so it cannot be inferred.\n", flag, rows);
    break;
  case YAME_REF_WRONG_KIND:
    fprintf(out, "%s is required: %" PRIu64 " rows is %s, whose row space is a "
                 "probe ordering rather than a coordinate stream.\n",
            flag, rows, name ? name : "an array");
    break;
  case YAME_REF_MISSING:
    fprintf(out, "%s is required: %" PRIu64 " rows is %s, but its reference is "
                 "not in the store.\n  Run `yame fetch %s`, or pass %s "
                 "yourself.\n",
            flag, rows, name ? name : "?", fetch ? fetch : "?", flag);
    break;
  default:
    break;
  }
}

/* ------------------------------------------------------ resolving by name */

/* The part of a filename before the first dot: "ChromHMM.20220303.cm" is the
 * ChromHMM set, whatever date it carries. */
static void set_name_of(const char *file, char *out, size_t n) {
  const char *dot = strchr(file, '.');
  size_t l = dot ? (size_t)(dot - file) : strlen(file);
  if (l >= n) l = n - 1;
  memcpy(out, file, l);
  out[l] = '\0';
}

/**
 * The newest file in `dir` whose set name is `want`.
 *
 * Newest by name, which for these is newest by date: the sets are versioned
 * as <set>.<YYYYMMDD>.cm, so lexicographic order is chronological order and
 * "ChromHMM" means the most recent ChromHMM rather than an arbitrary one.
 * A full filename matches too, so being exact still works.
 */
static int newest_named(const char *dir, const char *want, char *out,
                        size_t n) {
  DIR *d = opendir(dir);
  if (!d) return 0;

  char best[512] = "";
  struct dirent *e;
  while ((e = readdir(d))) {
    if (e->d_name[0] == '.') continue;
    if (strcmp(e->d_name, "SHA256SUMS") == 0) continue;
    /* Never an index: it shares its data file's set name and sorts after it,
     * so "genes" would resolve to genes.bed.gz.tbi. */
    if (yame_assets_index_suffix(e->d_name)) continue;

    char sn[256];
    set_name_of(e->d_name, sn, sizeof(sn));
    if (strcasecmp(sn, want) != 0 && strcasecmp(e->d_name, want) != 0) continue;
    if (strcmp(e->d_name, best) > 0) snprintf(best, sizeof(best), "%s", e->d_name);
  }
  closedir(d);

  if (!best[0]) return 0;
  return yame_assets_join(out, n, dir, best) == 0;
}

uint64_t yame_ref_file_rows(const char *path) {
  if (!path || !yame_assets_is_file(path)) return 0;
  cfile_t cf = open_cfile((char *)path);
  cdata_t c = read_cdata1(&cf);
  uint64_t rows = c.n ? cdata_n(&c) : 0;
  free_cdata(&c);
  bgzf_close(cf.fh);
  return rows;
}

int yame_ref_resolve(const char *spec, uint64_t rows, const char *store_override,
                     const char *want_kind, char *path, size_t n,
                     const char **name, const char **fetch) {
  if (name) *name = NULL;
  if (fetch) *fetch = NULL;
  if (!spec || !*spec) return YAME_REF_UNKNOWN;

  /* A path is a path. Only something that is not one is looked up, so the
   * ordinary spelling costs nothing and never changes meaning. */
  if (yame_assets_is_file(spec)) {
    snprintf(path, n, "%s", spec);
    return YAME_REF_OK;
  }

  const yame_ref_rows_t *e = NULL;
  for (size_t i = 0; i < YAME_REF_ROWS_N; ++i)
    if (YAME_REF_ROWS[i].rows == rows) { e = &YAME_REF_ROWS[i]; break; }
  if (!e) return YAME_REF_UNKNOWN;
  if (name) *name = e->name;
  if (fetch) *fetch = e->fetch;
  if (want_kind && strcmp(e->kind, want_kind) != 0) return YAME_REF_WRONG_KIND;

  char root[4096], dir[4096];
  yame_assets_root(store_override, NULL, root, sizeof(root));

  /* The row space's own name resolves to its reference: -R hg38. */
  if (strcasecmp(spec, e->name) == 0) {
    if (yame_assets_join(path, n, root, e->store_path) != 0)
      return YAME_REF_MISSING;
    return yame_assets_is_file(path) ? YAME_REF_OK : YAME_REF_MISSING;
  }

  /* Otherwise any file of that name the row space owns, its own directory
   * before its knowledgebase. Searching only the knowledgebase left the files
   * that come WITH a platform or a genome -- ordering, mask, coord, seqinfo,
   * genes -- unreachable by name, though the browser lists them in the same
   * unit as the sets. */
  for (size_t k = 0; e->dirs && e->dirs[k]; ++k)
    if (yame_assets_join(dir, sizeof(dir), root, e->dirs[k]) == 0 &&
        newest_named(dir, spec, path, n))
      return YAME_REF_OK;

  return YAME_REF_NO_NAME;
}

void yame_ref_explain_name(FILE *out, const char *spec, uint64_t rows,
                           int status, const char *name, const char *fetch,
                           const char *flag) {
  /* `-R -r chr16` parses as -R with the argument "-r": the flag may be
   * omitted, but if it is written it takes an argument, and getopt will take
   * whatever follows. Saying so beats reporting a reference named "-r". */
  if (spec && spec[0] == '-' && spec[1]) {
    fprintf(out, "%s takes an argument, and '%s' looks like the next option "
                 "rather than one.\n  Omit %s entirely to have it inferred, or "
                 "give it a path or a name.\n", flag, spec, flag);
    return;
  }

  switch (status) {
  case YAME_REF_UNKNOWN:
    fprintf(out, "%s: no file named '%s', and %" PRIu64 " rows matches no row "
                 "space this build knows, so the name cannot be looked up.\n",
            flag, spec, rows);
    break;
  case YAME_REF_WRONG_KIND:
    fprintf(out, "%s: no file named '%s', and %s is not the kind of row space "
                 "this takes.\n", flag, spec, name ? name : "that row space");
    break;
  case YAME_REF_NO_NAME:
    fprintf(out, "%s: no file named '%s', and nothing called that is in the "
                 "store for %s.\n  Run `yame fetch %s` to see what is "
                 "available, or give a path.\n",
            flag, spec, name ? name : "?", fetch ? fetch : "?");
    break;
  case YAME_REF_MISSING:
    fprintf(out, "%s: %s is not in the store.\n  Run `yame fetch %s`, or give "
                 "a path.\n", flag, name ? name : spec, fetch ? fetch : "?");
    break;
  default:
    break;
  }
}
