// SPDX-License-Identifier: LicenseRef-CHOP-Academic-BSD-2-Clause
/**
 * probe_ui: the parts of yame_ui.h that only a downstream caller reaches.
 *
 * yame itself never asks for a line of text inside a widget, and never takes
 * a signal while one is up in a test, so yame_ui_panel_ask and the raw-mode
 * signal handler had no coverage -- and kycg, which uses both, carries its UI
 * from this file. This is the smallest caller that reaches them: a one-row
 * tree whose `s` key opens a panel and asks, the way kycg asks for a store.
 * t_ui.sh drives it through a pty.
 *
 *   probe_ui pick
 *       the catalogue picker as kycg's test/annotate would open it, through
 *       yame_browse_pick_opt: only HM27 and MSA, a title and verb of its own
 *       (t test), HM27 open with CGI arriving checked. Prints "PICKED <n>".
 *
 *   PROBE_UI_LOG=<file> probe_ui [initial-value]
 *       s asks (the answer starts as initial-value, default /old), q quits;
 *       one "RESULT rc=<0|1> buf=[...]" line per ask goes to the log, since
 *       stdout is the pty the widget draws on.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "yame_ui.h"
#include "assets.h"
#include "registry.h"

static char initial[512] = "/old";

static int on_key(void *ctx, char key, const char *path, const char *node_key) {
  (void) ctx; (void) path; (void) node_key;
  if (key != 's') return 0;
  char buf[512];
  snprintf(buf, sizeof buf, "%s", initial);
  yame_ui_panel_open(2);
  int rc = yame_ui_panel_ask(1, "store:", buf, sizeof buf);
  yame_ui_panel_close();
  static FILE *log = NULL;
  if (!log) log = fopen(getenv("PROBE_UI_LOG"), "w");
  if (log) { fprintf(log, "RESULT rc=%d buf=[%s]\n", rc, buf); fflush(log); }
  return 0;
}

static int pick(void) {
  static const yame_fetch_cfg_t cfg = { YAME_FILES, YAME_FILES_N, "probe", NULL, 0 };
  static const char *const units[] = { "HM27", "MSA" };
  yame_pick_opt_t o;
  memset(&o, 0, sizeof o);
  o.units = units; o.n_units = 2;
  o.open_unit = "HM27";
  o.preselect = "cgi";                 /* set name, any case */
  o.title = "probe pick";
  o.verb_key = 't'; o.verb = "test";
  char **paths = NULL;
  size_t n = yame_browse_pick_opt(&cfg, &o, &paths);
  printf("PICKED %zu\n", n);
  for (size_t i = 0; i < n; ++i) free(paths[i]);
  free(paths);
  return 0;
}

int main(int argc, char *argv[]) {
  if (argc > 1 && strcmp(argv[1], "pick") == 0) return pick();
  if (argc > 1) snprintf(initial, sizeof initial, "%s", argv[1]);
  char *roots[] = { "item\tone" };
  unsigned char branch[] = { 0 };
  yame_ui_tree_t spec = {0};
  spec.title = "probe_ui";
  spec.roots = roots;
  spec.root_branch = branch;
  spec.n_roots = 1;
  spec.on_key = on_key;
  spec.hint = "s ask";
  int r = yame_ui_tree(&spec);
  printf("TREE %d\n", r);
  return r < 0 ? 1 : 0;
}
