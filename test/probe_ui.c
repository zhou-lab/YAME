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
 *   PROBE_UI_LOG=<file> probe_ui [initial-value]
 *       s asks (the answer starts as initial-value, default /old), q quits;
 *       one "RESULT rc=<0|1> buf=[...]" line per ask goes to the log, since
 *       stdout is the pty the widget draws on.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "yame_ui.h"

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

int main(int argc, char *argv[]) {
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
