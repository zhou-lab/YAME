// SPDX-License-Identifier: AGPL-3.0-or-later
/**
 * This file is part of YAME.
 *
 * Copyright (C) 2021-present Wanding Zhou
 *
 * YAME is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * YAME is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with YAME.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef _YAME_UI_H
#define _YAME_UI_H

#include <stdint.h>
#include <stddef.h>

/**
 * Terminal presentation and interaction.
 *
 * Everything here degrades to plain text off a TTY. See ui.c for why that is
 * a correctness property and not a nicety.
 */

/** Is the terminal capable of animation / color? (stderr is a TTY, TERM sane,
 *  NO_COLOR unset.) */
int yame_ui_fancy(void);

/** Can we ask the user a question? (stdin AND stderr are TTYs.) */
int yame_ui_interactive(void);

/** Does the locale advertise UTF-8? Selects glyphs vs ASCII fallbacks. */
int yame_ui_unicode(void);

/* SGR sequences, or "" when color is off. */
const char *yame_ui_dim(void);
const char *yame_ui_bold(void);
const char *yame_ui_green(void);
const char *yame_ui_red(void);
const char *yame_ui_yellow(void);
const char *yame_ui_cyan(void);
const char *yame_ui_reset(void);

/** Glyphs: check, cross, bullet. UTF-8 or ASCII depending on locale. */
const char *yame_ui_check(void);
const char *yame_ui_cross(void);
const char *yame_ui_bullet(void);

/** Human-readable byte count into `buf`; returns buf. */
const char *yame_ui_human(uint64_t bytes, char *buf, size_t n);

/* ------------------------------------------------------------- progress */

/**
 * A single file's transfer, rendered as one self-updating line while it runs
 * and one settled line when it finishes.
 */
typedef struct {
  char     label[128];
  uint64_t total;        /* 0 when the server does not say */
  uint64_t now;
  int      frame;        /* spinner position */
  int      active;       /* a line is currently on screen */
  double   started;      /* monotonic seconds, for the rate readout */
  double   last_draw;
} yame_prog_t;

void yame_prog_begin(yame_prog_t *p, const char *label, uint64_t total);
/** Safe to call at any rate; redraw is throttled internally. */
void yame_prog_update(yame_prog_t *p, uint64_t now, uint64_t total);
/** Settle the line: check mark on success, cross on failure. */
void yame_prog_done(yame_prog_t *p, const char *detail, int ok);

/** A one-off settled line for work that did not transfer (e.g. cached). */
void yame_ui_line(const char *glyph_color, const char *glyph,
                  const char *label, const char *detail);

/* -------------------------------------------------------------- prompts */

/**
 * Yes/no. Returns 1 for yes, 0 for no.
 *
 * Callers MUST check yame_ui_interactive() first and choose a non-interactive
 * default themselves; this function assumes it may block on stdin.
 */
int yame_ui_confirm(const char *question, int default_yes);

/**
 * Free-text answer with a default. Returns a malloc'd string (the default when
 * the user just hits return), or NULL on EOF.
 */
char *yame_ui_ask(const char *question, const char *def);

/**
 * Single-choice menu. Returns the chosen index, or -1 if cancelled.
 * `notes` may be NULL; when given, notes[i] is shown dimmed after items[i].
 *
 * On a capable terminal this is an in-place cursor list navigated with the
 * arrow keys; otherwise it degrades to a numbered prompt.
 */
long yame_ui_choose(const char *title, const char **items, const char **notes,
                    size_t n);

/**
 * Multi-select. Returns a malloc'd array of n flags, or NULL if cancelled.
 *
 * In-place on a capable terminal: arrows or j/k to move, space to toggle,
 * a/n for all/none, / to filter, enter to accept. Otherwise a numbered prompt
 * accepting "all", "none", and comma/range lists such as "1-5,8,12".
 */
int *yame_ui_multiselect(const char *title, const char **items,
                         const char **notes, size_t n, int default_all);

/**
 * Per-row emphasis.
 *
 * Passed alongside the rows rather than embedded in them: the widgets measure
 * and truncate row text to fit the terminal, and escape sequences inside that
 * text would be counted as visible cells and wreck the alignment.
 */
typedef enum {
  YAME_ROW_PLAIN = 0,
  YAME_ROW_HAVE,      /* present locally -- green */
  YAME_ROW_MISSING,   /* not present -- dimmed */
  YAME_ROW_REQUIRED,  /* comes with any fetch, not selectable -- red */
} yame_row_style_t;

/**
 * Scrollable in-place viewer for tabular output.
 *
 * `header` is a column header held fixed above the rows; it and `rows` are
 * tab-separated and rendered into aligned columns. `styles` may be NULL, or
 * one yame_row_style_t per row. Returns 0 when the viewer ran, -1 when the
 * terminal cannot support it and the caller should print plainly instead.
 * Never call this when stdout is redirected: piped output must stay
 * machine-readable.
 */
int yame_ui_browse(const char *title, const char *header,
                   const char **rows, const unsigned char *styles, size_t n);

/** Child rows of one expanded node, owned by the caller of the expand fn. */
typedef struct {
  char **rows;            /* preformatted lines; the tree indents, not aligns */
  char **keys;            /* optional; opaque, handed back on accept */
  unsigned char *styles;  /* optional, one per row */
  size_t n;
} yame_ui_kids_t;

/**
 * Called once per checked row when the user accepts, before the tree returns.
 * `root` is the parent's row text, `key` the child's key from yame_ui_kids_t.
 */
typedef void (*yame_ui_accept_fn)(void *ctx, const char *root,
                                  const char *key);

/**
 * Fill `out` with the children of the node whose row is `row`. Called at most
 * once per node, the first time it is opened. Leaving out->n at 0 marks the
 * node as having nothing to show.
 */
typedef void (*yame_ui_expand_fn)(void *ctx, const char *row,
                                  yame_ui_kids_t *out);

/**
 * Two-level tree viewer: a table whose rows unfold in place.
 *
 * Right arrow (or l, or enter) opens the row under the cursor and splices its
 * children in beneath it; left (or h) closes it. Children are requested lazily
 * and kept, so opening a row twice costs one call.
 *
 * Children are rendered as given, only indented — the tree cannot align them,
 * because it sees them one parent at a time and column widths that shifted
 * with each expansion would be worse than none. Pre-format them to a common
 * width if they should line up.
 *
 * Same contract as yame_ui_browse: returns -1 when the terminal cannot support
 * it, and must not be used when stdout is redirected.
 */

/* ------------------------------------------------------ inline detail pane */

/**
 * Lines describing the row under the cursor, rendered inside the tree's own
 * frame rather than over it.
 *
 * The pane is part of the redraw, not a modal overlay, which is the whole
 * point: the arrow keys keep moving the cursor while it is open and the pane
 * follows. A blocking panel meant you had to close it, move, and reopen it to
 * compare two sets -- which is the thing you most want to do with it.
 *
 * The callback owns nothing: it fills `rows` with malloc'd strings and sets
 * `n`, and the widget frees them after drawing. `cols` is the usable width, so
 * the callback can wrap; it is called on every redraw, so it must be cheap
 * (a table lookup and some formatting, not a file read).
 */
typedef struct {
  char **rows;
  size_t n;
} yame_ui_detail_t;

typedef void (*yame_ui_detail_fn)(void *ctx, const char *root,
                                  const char *child_key, int cols,
                                  yame_ui_detail_t *out);

/* ------------------------------------------------- panel inside a widget */

/**
 * A transient region at the foot of a full-screen widget, for work that
 * happens without leaving it: a plan, a confirmation, a progress bar.
 *
 * Suspending the widget to the normal screen for a download and coming back
 * would lose the user's place and make fetching feel like a separate errand.
 * The panel keeps the catalogue visible above whatever is being fetched.
 *
 * Line indices are relative to the top of the panel. Only meaningful while a
 * full-screen widget is running; harmless otherwise.
 */
void yame_ui_panel_open(int height);
void yame_ui_panel_line(int i, const char *fmt, ...);
void yame_ui_panel_close(void);
int  yame_ui_panel_active(void);

/**
 * Usable terminal width, already clamped to the range the widgets assume.
 * Callers that lay out prose in a panel need it: a panel line is truncated at
 * the right margin, so text long enough to matter has to be wrapped by
 * whoever knows where it may be broken.
 */
int yame_ui_cols(void);

/**
 * Terminal height. A panel is clamped to two rows less than this and any
 * line past its height is silently dropped -- so a caller laying out more
 * text than fits must truncate deliberately, or the key-to-dismiss prompt is
 * what disappears and the widget looks hung.
 */
int yame_ui_rows(void);

/** Yes/no on one panel line, answered with a single keypress. */
int yame_ui_panel_confirm(int line, const char *question, int default_yes);

/** Show a message on a panel line and wait for any key. */
void yame_ui_panel_pause(int line, const char *msg);

/**
 * Edit a line of text on one panel line. `buf` carries the current value in
 * and the new one out. Returns 1 if accepted, 0 if cancelled.
 */
int yame_ui_panel_ask(int line, const char *prompt, char *buf, size_t n);

/**
 * Called once, after every checked row has been reported through `accept`.
 *
 * Runs while the widget is still on screen; use the panel calls above for any
 * output. Whatever it changes in the `roots` and `root_styles` arrays is
 * picked up when the tree resumes — but it must replace the entries in place
 * rather than reallocating the arrays, whose addresses the tree holds.
 */
typedef void (*yame_ui_commit_fn)(void *ctx);

/**
 * One thing that can be done with the checked rows.
 *
 * `accept` is called per checked row, then `commit`. A NULL commit means the
 * action ends the session and the caller acts on what it collected; a non-NULL
 * one runs inside the widget (use the panel calls) and the tree carries on,
 * reloaded, so several actions can be taken in one sitting.
 */
typedef struct {
  char               key;
  const char        *verb;    /* shown in the footer */
  yame_ui_accept_fn  accept;
  yame_ui_commit_fn  commit;  /* NULL: the action returns instead */
} yame_ui_action_t;

/**
 * Offered any key the tree does not handle itself, so a caller can bind
 * actions of its own without the widget knowing what they mean. Receives the
 * row under the cursor -- its parent's text, and its own key when it is a
 * child -- so a binding can act on what is being looked at. Runs with the
 * widget still on screen; use the panel calls for output. Return nonzero if
 * the roots changed and the view should reload.
 */
typedef int (*yame_ui_key_fn)(void *ctx, char key, const char *root,
                              const char *child_key);

/**
 * Should this child start checked? Called once per child of the preselected
 * root as the tree opens, so a named target can arrive already chosen and the
 * user only has to look at it and press f.
 */
typedef int (*yame_ui_preselect_fn)(void *ctx, const char *root,
                                    const char *key);

/**
 * Everything yame_ui_tree() needs. A struct because the list outgrew a
 * readable argument list, and because most fields are optional.
 */
typedef struct {
  const char *title;
  const char *header;          /* tab-separated column names, may be NULL */

  /* Not const: `commit` is allowed to rewrite entries in place so the view
   * reflects what it just did. The arrays themselves must not move. */
  char          **roots;
  unsigned char  *root_styles; /* may be NULL */
  size_t          n_roots;

  yame_ui_expand_fn expand;    /* may be NULL for a flat tree */

  /* What can be done with a selection. Any action makes the tree a picker.
   * More than one lets a single screen serve more than one verb -- fetching
   * what is missing and then testing against it, without leaving. */
  yame_ui_action_t  actions[4];
  size_t            n_actions;

  yame_ui_key_fn    on_key;    /* may be NULL */

  /* Inline detail pane. `detail_key` toggles it; while open it redraws for
   * whatever the cursor is on. NULL `detail` disables the whole feature. */
  char              detail_key;
  const char       *detail_verb;   /* footer label, e.g. "info" */
  yame_ui_detail_fn detail;
  const char       *hint;      /* extra footer text for the caller's keys */

  /* Whether a row already present can still be checked. False for fetching,
   * where having it means there is nothing to ask for; true when the action
   * is something other than acquiring it -- testing against it, say. */
  int               have_selectable;

  /* Open this root on entry and offer its children to `preselect`. NULL for
   * the plain catalogue view. */
  const char          *open_root;
  yame_ui_preselect_fn preselect;

  /* Bound to `r` when set: opens every collection and checks the rows this
   * accepts. A curated default matters because the catalogue is large and
   * most of it is situational -- the answer to "which of these should I
   * actually use" should be one keystroke, not a reading exercise. */
  yame_ui_preselect_fn recommend;
  void *ctx;
} yame_ui_tree_t;

/**
 * Two-level tree viewer: a table whose rows unfold in place, optionally with
 * checkboxes.
 *
 * When `accept` is set, rows carrying a key and not already marked
 * YAME_ROW_HAVE get a checkbox; space toggles one, space on a parent toggles
 * all of its children, and `f` commits. Rows styled YAME_ROW_HAVE are shown as
 * already present and cannot be checked — there is nothing to ask for.
 *
 * With `commit` set, `f` does not end the session: the widget suspends, the
 * callback runs on the normal screen, and the tree resumes with its children
 * reloaded so the results are visible. Only `q` leaves.
 *
 * Returns:
 *   -1  the terminal cannot host the widget; the caller should print plainly
 *    0  the user quit without invoking an action
 *   >0  the 1-based index into `actions` of the action they invoked, so a
 *       caller with several verbs can tell which one ended the session
 *       (kycg's browser checks for 2 to tell `t test` from `f fetch`)
 *
 * An action with a `commit` does not end the session, so it never becomes the
 * return value -- only an action whose commit is NULL does.
 */
int yame_ui_tree(const yame_ui_tree_t *spec);

#endif /* _YAME_UI_H */
