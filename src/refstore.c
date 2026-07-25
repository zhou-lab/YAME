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

#include "assets.h"
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
