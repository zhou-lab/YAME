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

/**
 * Terminal presentation: colors, spinners, progress bars, and prompts.
 *
 * WHY THIS IS A SEPARATE MODULE
 *   Every routine here has to answer the same question first -- is anyone
 *   actually watching? -- and getting that wrong is not a cosmetic bug.
 *   These tools are meant to run inside Nextflow pipelines and Docker builds,
 *   where
 *   stdin is closed and stderr is a log file. A prompt in that setting hangs
 *   the job forever with no indication of why, and an animated progress bar
 *   fills the log with thousands of carriage returns. Centralizing the
 *   capability checks means each caller cannot forget them independently.
 *
 * THE POLICY
 *   Interaction requires stdin AND stderr to both be TTYs. Animation requires
 *   stderr to be a TTY. Neither is ever assumed. Callers that want to prompt
 *   must check yame_ui_interactive() and supply their own non-interactive
 *   behavior; nothing here silently blocks on a closed stdin.
 *
 *   This is the refinement of the original "never prompt" rule. That rule
 *   existed because a prompt would hang automation; gating prompts on an
 *   interactive terminal preserves the guarantee exactly, while letting a
 *   human at a keyboard get a usable interface.
 *
 * DEGRADATION
 *   Off a TTY: no escape sequences, no spinner, one plain line per event.
 *   Without a UTF-8 locale: ASCII glyphs ([ok], [xx], -) and an ASCII spinner.
 *   With NO_COLOR set or TERM=dumb: the full-screen widgets are skipped
 *   entirely and callers print plainly instead.
 *   The information content is identical in every mode -- only the ink differs.
 *
 * TWO KINDS OF DRAWING
 *   The progress bar is inline: it owns one line of the normal buffer and
 *   rewrites it, because a download is something you watch alongside whatever
 *   else is on screen.
 *
 *   The selector, browser and tree are full-screen: raw_enter() switches to
 *   the alternate screen buffer, so each frame is absolutely positioned from
 *   the home cell and the terminal restores the user's screen untouched on
 *   exit. That is what makes a fixed-height viewport reasonable -- a
 *   full-height widget drawn inline would shove the user's scrollback away --
 *   and it removes a whole class of redraw bug, since a frame that shrinks
 *   cannot leave the tail of a taller predecessor behind it.
 */

/* strcasestr is a GNU extension: the filter box matches case-insensitively,
 * and this is the one call that needs the feature macro. */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "yame_ui.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <strings.h>
#include <sys/select.h>
#include <unistd.h>
#include <ctype.h>
#include <time.h>
#include <signal.h>
#include <termios.h>
#include <sys/time.h>
#include <sys/ioctl.h>

/* ------------------------------------------------------- capability probes */

static int cached_fancy = -1, cached_inter = -1, cached_uni = -1;

int yame_ui_fancy(void) {
  if (cached_fancy >= 0) return cached_fancy;

  const char *term = getenv("TERM");
  int ok = isatty(STDERR_FILENO)
        && !getenv("NO_COLOR")
        && !(term && strcmp(term, "dumb") == 0);
  cached_fancy = ok;
  return ok;
}

int yame_ui_interactive(void) {
  if (cached_inter >= 0) return cached_inter;
  /* Both directions matter: a question we cannot display is as useless as one
   * whose answer we cannot read. */
  cached_inter = isatty(STDIN_FILENO) && isatty(STDERR_FILENO);
  return cached_inter;
}

int yame_ui_unicode(void) {
  if (cached_uni >= 0) return cached_uni;
  const char *v = getenv("LC_ALL");
  if (!v || !*v) v = getenv("LC_CTYPE");
  if (!v || !*v) v = getenv("LANG");
  cached_uni = (v && (strcasestr(v, "utf-8") || strcasestr(v, "utf8"))) ? 1 : 0;
  return cached_uni;
}

/* --------------------------------------------------------------- palette */

#define COLOR(fn, seq) \
  const char *fn(void) { return yame_ui_fancy() ? seq : ""; }

COLOR(yame_ui_dim,    "\033[2m")
COLOR(yame_ui_bold,   "\033[1m")
COLOR(yame_ui_green,  "\033[32m")
COLOR(yame_ui_red,    "\033[31m")
COLOR(yame_ui_yellow, "\033[33m")
COLOR(yame_ui_cyan,   "\033[36m")
COLOR(yame_ui_blue,   "\033[34m")
COLOR(yame_ui_reset,  "\033[0m")

const char *yame_ui_check(void)  { return yame_ui_unicode() ? "✓" : "ok"; }
const char *yame_ui_cross(void)  { return yame_ui_unicode() ? "✗" : "XX"; }
const char *yame_ui_bullet(void) { return yame_ui_unicode() ? "•" : "-"; }

/* Braille dots read as a smooth rotation; the ASCII fallback is the classic
 * four-frame spin. */
static const char *SPIN_U[] = {"⠋","⠙","⠹","⠸","⠼",
                               "⠴","⠦","⠧","⠇","⠏"};
static const char *SPIN_A[] = {"|","/","-","\\"};

static const char *spin_frame(int i) {
  if (yame_ui_unicode()) return SPIN_U[((unsigned)i) % 10];
  return SPIN_A[((unsigned)i) % 4];
}

/* ---------------------------------------------------------------- helpers */

const char *yame_ui_human(uint64_t bytes, char *buf, size_t n) {
  const char *unit[] = {"B", "KB", "MB", "GB", "TB"};
  double v = (double)bytes;
  int u = 0;
  while (v >= 1024.0 && u < 4) { v /= 1024.0; ++u; }
  snprintf(buf, n, "%.*f %s", (u == 0 || v >= 100) ? 0 : 1, v, unit[u]);
  return buf;
}

static double now_sec(void) {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (double)tv.tv_sec + (double)tv.tv_usec / 1e6;
}

/* Terminal width, clamped to something sane. */
static int term_cols(void) {
  /* Ask the terminal first, exactly as term_rows does. $COLUMNS is a shell
   * convenience that bash does not export to children, so consulting it alone
   * meant every widget rendered at the 80-column fallback no matter how wide
   * the window actually was. */
  struct winsize ws;
  int w = 0;
  if (ioctl(STDERR_FILENO, TIOCGWINSZ, &ws) == 0 && ws.ws_col > 0)
    w = ws.ws_col;
  if (w <= 0) {
    const char *c = getenv("COLUMNS");
    w = c ? atoi(c) : 0;
  }
  if (w <= 0) w = 80;
  if (w < 40) w = 40;
  if (w > 200) w = 200;
  return w;
}

/* Erase the current line, leaving the cursor at column 0. */
static void clear_line(void) {
  fputs("\r\033[2K", stderr);
}

/* Terminal height, clamped. */
static int term_rows(void) {
  struct winsize ws;
  if (ioctl(STDERR_FILENO, TIOCGWINSZ, &ws) == 0 && ws.ws_row > 0)
    return ws.ws_row;
  const char *l = getenv("LINES");
  int r = l ? atoi(l) : 0;
  return r > 0 ? r : 24;
}

int yame_ui_cols(void) { return term_cols(); }
int yame_ui_rows(void) { return term_rows(); }

/* --------------------------------------------------------------- raw mode */

/*
 * The in-place widgets need keys as they are pressed, so the terminal goes
 * into raw mode: no line buffering, no echo. That is a global change to
 * something the user owns, and leaving it applied would hand back a shell
 * with no echo and no working Ctrl-C -- so restoration is wired to atexit and
 * to the fatal signals as well as to the normal return path.
 */
static struct termios saved_tio;
static int raw_active = 0;
static volatile sig_atomic_t interrupted = 0;

static void raw_leave(void) {
  if (!raw_active) return;
  /* Leave the alternate screen first, so the terminal restores whatever the
   * user had before we drew over it, then hand back cooked mode. */
  fputs("\033[?1049l\033[?25h", stderr);
  tcsetattr(STDIN_FILENO, TCSAFLUSH, &saved_tio);
  fflush(stderr);
  raw_active = 0;
}

static void on_signal(int sig) {
  interrupted = 1;
  raw_leave();
  signal(sig, SIG_DFL);
  raise(sig);
}

static int raw_enter(void) {
  if (!yame_ui_interactive()) return -1;
  if (tcgetattr(STDIN_FILENO, &saved_tio) != 0) return -1;

  struct termios t = saved_tio;
  t.c_lflag &= (tcflag_t)~(ICANON | ECHO);
  t.c_cc[VMIN] = 1;
  t.c_cc[VTIME] = 0;
  if (tcsetattr(STDIN_FILENO, TCSAFLUSH, &t) != 0) return -1;

  raw_active = 1;
  interrupted = 0;

  static int hooked = 0;
  if (!hooked) {
    atexit(raw_leave);
    signal(SIGINT, on_signal);
    signal(SIGTERM, on_signal);
    hooked = 1;
  }

  /* Alternate screen: the widget gets the whole terminal and gives it back
   * untouched on exit, which is what makes a fixed-height viewport reasonable
   * rather than something that shoves the user's scrollback off the top. */
  fputs("\033[?1049h\033[H\033[2J\033[?25l", stderr);
  fflush(stderr);
  return 0;
}

/* ------------------------------------------------------------ key decoding */

typedef enum {
  K_NONE = 0, K_UP, K_DOWN, K_LEFT, K_RIGHT, K_PGUP, K_PGDN,
  K_HOME, K_END, K_ENTER, K_SPACE, K_ESC, K_BACKSPACE, K_CHAR
} keycode_t;

/**
 * One keypress. Arrow and navigation keys arrive as escape sequences
 * (ESC [ A and friends); a bare ESC is distinguished from the start of a
 * sequence by a short poll rather than by blocking, so pressing Escape does
 * not hang waiting for a second byte that never comes.
 */
static keycode_t read_key(char *out) {
  unsigned char c;
  if (read(STDIN_FILENO, &c, 1) != 1) return K_NONE;

  if (c == '\r' || c == '\n') return K_ENTER;
  if (c == ' ')  return K_SPACE;
  if (c == 127 || c == 8) return K_BACKSPACE;

  if (c == 27) {
    struct timeval tv = {0, 40000};    /* 40 ms is well above key repeat */
    fd_set fds;
    FD_ZERO(&fds);
    FD_SET(STDIN_FILENO, &fds);
    if (select(STDIN_FILENO + 1, &fds, NULL, NULL, &tv) <= 0) return K_ESC;

    unsigned char a, b;
    if (read(STDIN_FILENO, &a, 1) != 1) return K_ESC;
    if (a != '[' && a != 'O') return K_ESC;
    if (read(STDIN_FILENO, &b, 1) != 1) return K_ESC;

    switch (b) {
    case 'A': return K_UP;
    case 'B': return K_DOWN;
    case 'C': return K_RIGHT;
    case 'D': return K_LEFT;
    case 'H': return K_HOME;
    case 'F': return K_END;
    default: break;
    }
    if (b >= '0' && b <= '9') {        /* ESC [ params... final */
      /* Consume the whole parameter list, up to the final byte (0x40-0x7E).
       * Reading exactly one more byte was right for ESC [ 5 ~ but wrong for a
       * modified key: Ctrl-Up is ESC [ 1 ; 5 A, so "5A" was left in the buffer
       * and delivered as two stray keypresses -- the cursor jumped and a
       * stray 'A' reached the action dispatch. */
      unsigned char params[16];
      size_t np = 0;
      params[np++] = b;
      unsigned char t = 0;
      for (;;) {
        if (read(STDIN_FILENO, &t, 1) != 1) return K_NONE;
        if (t >= 0x40 && t <= 0x7E) break;          /* final byte */
        if (np < sizeof(params)) params[np++] = t;
      }

      if (t == '~') {
        switch (params[0]) {
        case '5': return K_PGUP;
        case '6': return K_PGDN;
        case '1': case '7': return K_HOME;
        case '4': case '8': return K_END;
        default: return K_NONE;
        }
      }
      /* A modified arrow ends in the same final byte as the plain one, so
       * treat it as the unmodified key rather than discarding it. */
      switch (t) {
      case 'A': return K_UP;
      case 'B': return K_DOWN;
      case 'C': return K_RIGHT;
      case 'D': return K_LEFT;
      case 'H': return K_HOME;
      case 'F': return K_END;
      default: return K_NONE;
      }
    }
    return K_NONE;
  }

  if (out) *out = (char)c;
  return K_CHAR;
}

/* ----------------------------------------------------- in-place rendering */

/*
 * A frame is drawn as N lines and erased by moving the cursor back up N lines
 * and clearing each one. Tracking the count we actually wrote -- rather than
 * assuming it -- is what keeps the widget stable when the list is shorter than
 * the viewport or the terminal is resized mid-run.
 */
typedef struct { int lines; } frame_t;

/*
 * Widgets draw on the alternate screen buffer, so a frame is absolutely
 * positioned: home the cursor, write the frame, erase whatever is below it.
 *
 * The earlier approach redrew in the normal buffer by counting lines back up.
 * That only works while frames stay the same height, and these frames do not:
 * folding a tree node or narrowing a search makes the next frame shorter, and
 * the tail of the taller one stayed on screen underneath it. Absolute
 * positioning plus erase-to-end removes the whole class of bug rather than
 * patching the arithmetic, and it is what a full-screen viewer wants anyway.
 */
static void frame_begin(frame_t *f) {
  fputs("\033[H", stderr);
  f->lines = 0;
}

static void frame_finish(frame_t *f) {
  (void)f;
  fputs("\033[J", stderr);   /* clear from the cursor to the end of screen */
}

static void frame_line(frame_t *f, const char *fmt, ...) {
  va_list ap;
  fputs("\r\033[2K", stderr);
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
  fputc('\n', stderr);
  ++f->lines;
}

/* Display cells in a UTF-8 string. Padding a column by bytes puts every
 * column after a multi-byte glyph two cells to the left, which is what a
 * heading marker like "\u25aa " did to every row it headed. */
static int cells_of(const char *s) {
  int n = 0;
  for (size_t i = 0; s[i]; ++i)
    if ((s[i] & 0xC0) != 0x80) ++n;     /* count lead bytes only */
  return n;
}

/* Truncate to `cols` display cells, ASCII-safe and multi-byte tolerant. */
static void fit(const char *s, int cols, char *out, size_t outsz) {
  int cells = 0;
  size_t i = 0, o = 0;
  while (s[i] && cells < cols && o + 4 < outsz) {
    unsigned char c = (unsigned char)s[i];
    size_t len = 1;
    if      ((c & 0xE0) == 0xC0) len = 2;
    else if ((c & 0xF0) == 0xE0) len = 3;
    else if ((c & 0xF8) == 0xF0) len = 4;
    for (size_t k = 0; k < len && s[i]; ++k) out[o++] = s[i++];
    ++cells;
  }
  out[o] = '\0';
}

/**
 * Lay a tab-separated line into fixed-width columns.
 *
 * Widths are in display cells; a field too wide for its column is truncated
 * rather than allowed to push the ones after it.
 */
static void columns(const char *line, const int *w, int ncol,
                    char *out, size_t cap) {
  size_t o = 0;
  out[0] = '\0';

  const char *p = line;
  for (int c = 0; c < ncol; ++c) {
    if (!*p) break;
    const char *tab = strchr(p, '\t');
    size_t flen = tab ? (size_t)(tab - p) : strlen(p);

    char field[512];
    size_t n = flen < sizeof(field) - 1 ? flen : sizeof(field) - 1;
    memcpy(field, p, n);
    field[n] = '\0';

    char cut[512];
    fit(field, w[c] > 2 ? w[c] - 2 : 1, cut, sizeof(cut));
    size_t l = strlen(cut);
    int used = cells_of(cut);
    int pad = w[c] - used;
    if (pad < 0) pad = 0;

    if (o + l + (size_t)pad + 1 >= cap) break;
    memcpy(out + o, cut, l);
    o += l;
    while (pad-- > 0) out[o++] = ' ';
    out[o] = '\0';

    if (!tab) break;
    p = tab + 1;
  }
}

/* ----------------------------------------------- panel inside a widget */

/*
 * The panel is absolutely positioned at the foot of the alternate screen, so
 * a download can render over the bottom few rows while the catalogue stays
 * visible above it. Without this, fetching would mean dropping the widget,
 * printing, and coming back -- losing the user's place every time.
 */
static int panel_top = 0;   /* 1-based screen row of the panel's first line */
static int panel_h = 0;

int yame_ui_panel_active(void) { return panel_h > 0; }

void yame_ui_panel_open(int height) {
  if (!yame_ui_fancy()) return;
  int rows = term_rows();
  if (height < 1) height = 1;
  if (height > rows - 2) height = rows - 2;
  panel_h = height;
  panel_top = rows - height + 1;
  for (int i = 0; i < panel_h; ++i)
    fprintf(stderr, "\033[%d;1H\033[2K", panel_top + i);
  fflush(stderr);
}

/* Where a panel line actually lands. panel_open clamps the panel to the
 * terminal, so a caller's fixed line index can fall outside it on a short
 * window; for an *interactive* line, dropping the write leaves the user
 * staring at nothing while read_key() blocks, which reads as a hang. Prompts
 * therefore clamp to the last row rather than disappear. */
static int panel_clamp(int i) {
  if (i < 0) return 0;
  if (i >= panel_h) return panel_h - 1;
  return i;
}

void yame_ui_panel_line(int i, const char *fmt, ...) {
  if (panel_h <= 0 || i < 0 || i >= panel_h) return;
  fprintf(stderr, "\033[%d;1H\033[2K", panel_top + i);
  va_list ap;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
  fflush(stderr);
}

void yame_ui_panel_close(void) { panel_h = 0; panel_top = 0; }

int yame_ui_panel_confirm(int line, const char *question, int default_yes) {
  line = panel_clamp(line);
  for (;;) {
    yame_ui_panel_line(line, "  %s%s%s %s[%s]%s ",
                       yame_ui_bold(), question, yame_ui_reset(),
                       yame_ui_dim(), default_yes ? "Y/n" : "y/N",
                       yame_ui_reset());
    char ch = 0;
    keycode_t k = read_key(&ch);
    if (k == K_ENTER) return default_yes;
    if (k == K_ESC || k == K_NONE) return 0;
    if (k == K_CHAR) {
      if (ch == 'y' || ch == 'Y') return 1;
      if (ch == 'n' || ch == 'N' || ch == 'q') return 0;
    }
  }
}

void yame_ui_panel_pause(int line, const char *msg) {
  line = panel_clamp(line);
  if (panel_h <= 0) return;
  yame_ui_panel_line(line, "  %s%s%s", yame_ui_dim(), msg, yame_ui_reset());
  char ch = 0;
  read_key(&ch);
}

int yame_ui_panel_ask(int line, const char *prompt, char *buf, size_t n) {
  size_t len = strlen(buf);
  for (;;) {
    int cols = term_cols();
    /* Show the tail when the value is longer than the room for it, so the
     * part being edited stays visible. */
    int room = cols - (int)strlen(prompt) - 6;
    if (room < 10) room = 10;
    const char *shown = (int)len > room ? buf + len - room : buf;

    yame_ui_panel_line(line, "  %s%s%s %s%s%s_",
                       yame_ui_bold(), prompt, yame_ui_reset(),
                       yame_ui_cyan(), shown, yame_ui_reset());

    char ch = 0;
    keycode_t k = read_key(&ch);
    if (k == K_ENTER) return 1;
    if (k == K_ESC || k == K_NONE) return 0;
    if (k == K_BACKSPACE) { if (len) buf[--len] = '\0'; continue; }
    if (k == K_CHAR && (unsigned char)ch >= 32 && len + 1 < n) {
      buf[len++] = ch;
      buf[len] = '\0';
    }
  }
}

/* ------------------------------------------------------------- progress */

#define LABEL_W 34
#define BAR_W   22
#define REDRAW_INTERVAL 0.08   /* seconds; ~12 fps is smooth and cheap */

void yame_prog_begin(yame_prog_t *p, const char *label, uint64_t total) {
  memset(p, 0, sizeof(*p));
  snprintf(p->label, sizeof(p->label), "%s", label);
  p->total = total;
  p->started = now_sec();
  p->last_draw = 0.0;
  p->active = yame_ui_fancy();
}

/* Where a progress line goes: its own panel row inside a widget, or the
 * current line of the normal screen. */
#define PANEL_ROW_DONE 1
#define PANEL_ROW_BAR  2

void yame_prog_update(yame_prog_t *p, uint64_t now, uint64_t total) {
  if (!p->active) return;

  p->now = now;
  if (total) p->total = total;

  double t = now_sec();
  if (t - p->last_draw < REDRAW_INTERVAL) return;
  p->last_draw = t;
  ++p->frame;

  char hb_now[24], hb_tot[24], hb_rate[24];
  yame_ui_human(p->now, hb_now, sizeof(hb_now));

  double elapsed = t - p->started;
  double rate = elapsed > 0.05 ? (double)p->now / elapsed : 0.0;
  yame_ui_human((uint64_t)rate, hb_rate, sizeof(hb_rate));

  int cols = term_cols();
  if (yame_ui_panel_active())
    fprintf(stderr, "\033[%d;1H\033[2K", panel_top + PANEL_ROW_BAR);
  else
    clear_line();

  if (p->total) {
    double frac = (double)p->now / (double)p->total;
    if (frac > 1.0) frac = 1.0;
    int fill = (int)(frac * BAR_W + 0.5);

    yame_ui_human(p->total, hb_tot, sizeof(hb_tot));

    /* The bar is drawn with block glyphs where the locale allows and '#'/'.'
     * where it does not; both are the same width in cells. */
    const char *full = yame_ui_unicode() ? "█" : "#";
    const char *rest = yame_ui_unicode() ? "░" : ".";

    fprintf(stderr, "%s%s%s %-*.*s ", yame_ui_cyan(), spin_frame(p->frame),
            yame_ui_reset(), LABEL_W, LABEL_W, p->label);

    fputs(yame_ui_cyan(), stderr);
    for (int i = 0; i < BAR_W; ++i) fputs(i < fill ? full : rest, stderr);
    fputs(yame_ui_reset(), stderr);

    fprintf(stderr, " %3d%%  %s%s/%s", (int)(frac * 100.0 + 0.5),
            yame_ui_dim(), hb_now, hb_tot);
    if (cols > 96 && rate > 0.0) fprintf(stderr, "  %s/s", hb_rate);
    fputs(yame_ui_reset(), stderr);
  } else {
    /* No Content-Length: show what has arrived, without a fake percentage. */
    fprintf(stderr, "%s%s%s %-*.*s %s%s%s",
            yame_ui_cyan(), spin_frame(p->frame), yame_ui_reset(),
            LABEL_W, LABEL_W, p->label,
            yame_ui_dim(), hb_now, yame_ui_reset());
  }

  fflush(stderr);
}

void yame_prog_done(yame_prog_t *p, const char *detail, int ok) {
  if (yame_ui_panel_active()) {
    /* Inside a widget the settled line replaces the bar rather than scrolling
     * past it: there is no scrollback to keep it in. */
    yame_ui_panel_line(PANEL_ROW_BAR, " ");
    yame_ui_panel_line(PANEL_ROW_DONE, "  %s%s%s %-40.40s %s%s%s",
                       ok ? yame_ui_green() : yame_ui_red(),
                       ok ? yame_ui_check() : yame_ui_cross(),
                       yame_ui_reset(), p->label,
                       yame_ui_dim(), detail ? detail : "", yame_ui_reset());
    p->active = 0;
    return;
  }
  if (p->active) clear_line();
  yame_ui_line(ok ? yame_ui_green() : yame_ui_red(),
               ok ? yame_ui_check() : yame_ui_cross(),
               p->label, detail);
  p->active = 0;
}

void yame_ui_line(const char *glyph_color, const char *glyph,
                  const char *label, const char *detail) {
  fprintf(stderr, "%s%s%s %-*.*s %s%s%s\n",
          glyph_color ? glyph_color : "", glyph ? glyph : " ", yame_ui_reset(),
          LABEL_W, LABEL_W, label ? label : "",
          yame_ui_dim(), detail ? detail : "", yame_ui_reset());
  fflush(stderr);
}

/* --------------------------------------------------------------- prompts */

/** Read a line from stdin, trimmed. Returns NULL on EOF. */
static char *read_line(void) {
  char buf[4096];
  if (!fgets(buf, sizeof(buf), stdin)) return NULL;

  size_t n = strlen(buf);
  while (n && (buf[n-1] == '\n' || buf[n-1] == '\r')) buf[--n] = '\0';

  char *p = buf;
  while (*p && isspace((unsigned char)*p)) ++p;
  size_t len = strlen(p);
  while (len && isspace((unsigned char)p[len-1])) p[--len] = '\0';

  return strdup(p);
}

int yame_ui_confirm(const char *question, int default_yes) {
  for (;;) {
    fprintf(stderr, "%s%s%s %s [%s] ",
            yame_ui_bold(), question, yame_ui_reset(),
            yame_ui_dim(), default_yes ? "Y/n" : "y/N");
    fputs(yame_ui_reset(), stderr);
    fflush(stderr);

    char *ans = read_line();
    if (!ans) { fputc('\n', stderr); return 0; }   /* EOF declines */

    if (!*ans) { free(ans); return default_yes; }
    int y = (strcasecmp(ans, "y") == 0 || strcasecmp(ans, "yes") == 0);
    int n = (strcasecmp(ans, "n") == 0 || strcasecmp(ans, "no") == 0);
    free(ans);

    if (y) return 1;
    if (n) return 0;
    fprintf(stderr, "  %sPlease answer y or n.%s\n",
            yame_ui_yellow(), yame_ui_reset());
  }
}

char *yame_ui_ask(const char *question, const char *def) {
  fprintf(stderr, "%s%s%s", yame_ui_bold(), question, yame_ui_reset());
  if (def && *def)
    fprintf(stderr, "\n  %s[%s]%s ", yame_ui_dim(), def, yame_ui_reset());
  else
    fputs(" ", stderr);
  fflush(stderr);

  char *ans = read_line();
  if (!ans) return NULL;
  if (!*ans && def) { free(ans); return strdup(def); }
  return ans;
}

/** Render a numbered list, one item per line. */
static void print_items(const char **items, const char **notes, size_t n) {
  int width = 1;
  for (size_t t = n; t >= 10; t /= 10) ++width;

  for (size_t i = 0; i < n; ++i) {
    fprintf(stderr, "  %s%*zu%s  %s",
            yame_ui_cyan(), width, i + 1, yame_ui_reset(), items[i]);
    if (notes && notes[i] && *notes[i])
      fprintf(stderr, "  %s%s%s", yame_ui_dim(), notes[i], yame_ui_reset());
    fputc('\n', stderr);
  }
}

/* ------------------------------------------------- the in-place selector */

/**
 * Shared engine for single- and multi-select.
 *
 * Draws a scrolling viewport with a cursor, redrawn in place on every
 * keypress. `flags` carries the selection in and out for multi-select mode and
 * is ignored for single-select. Returns the cursor index on accept, -1 on
 * cancel, or -2 when the terminal cannot do this and the caller should fall
 * back to the numbered prompt.
 *
 * The filter (/) matters more than it looks: a populated store lists dozens of
 * sets, and scrolling to "TFBSrm" past thirty neighbours is worse than typing
 * three characters.
 */
static long select_widget(const char *title, const char **items,
                          const char **notes, size_t n, int *flags, int multi) {
  if (!n) return -1;
  if (!yame_ui_fancy() || raw_enter() != 0) return -2;

  size_t *view = malloc(n * sizeof(size_t));   /* indices passing the filter */
  if (!view) { raw_leave(); return -2; }

  char filter[128] = {0};
  int filtering = 0;
  size_t cur = 0, top = 0, nview = 0;
  frame_t f = {0};
  long result = -1;
  int done = 0;

  while (!done) {
    /* Rebuild the visible subset. */
    nview = 0;
    for (size_t i = 0; i < n; ++i)
      if (!filter[0] || strcasestr(items[i], filter)) view[nview++] = i;

    if (cur >= nview) cur = nview ? nview - 1 : 0;

    int rows = term_rows();
    int cols = term_cols();
    /* Fixed height: the viewport is the screen, not the content. A window
     * that resized itself as the list grew and shrank was the source of the
     * leftover-line bug, and it makes scrolling feel unanchored. */
    int avail = rows - 4;                       /* title, filter, footer */
    if (avail < 3) avail = 3;

    if (cur < top) top = cur;
    if (avail > 0 && cur >= top + (size_t)avail) top = cur - (size_t)avail + 1;
    if (top + (size_t)avail > nview) top = nview > (size_t)avail
                                          ? nview - (size_t)avail : 0;

    frame_begin(&f);

    size_t chosen = 0;
    if (multi) for (size_t i = 0; i < n; ++i) chosen += (size_t)(flags[i] != 0);

    frame_line(&f, "%s%s%s", yame_ui_bold(), title, yame_ui_reset());

    if (filtering || filter[0])
      frame_line(&f, "  %ssearch:%s %s%s%s%s", yame_ui_dim(), yame_ui_reset(),
                 yame_ui_cyan(), filter, yame_ui_reset(),
                 filtering ? "_" : "");
    else
      frame_line(&f, "");

    for (size_t k = 0; k < (size_t)avail; ++k) {
      size_t vi = top + k;
      if (vi >= nview) { frame_line(&f, ""); continue; }
      size_t i = view[vi];

      const char *mark = "  ";
      if (multi) mark = flags[i] ? "[x]" : "[ ]";

      char label[512];
      int budget = cols - 10 - (notes && notes[i] ? (int)strlen(notes[i]) + 2 : 0);
      if (budget < 12) budget = 12;
      fit(items[i], budget, label, sizeof(label));

      if (vi == cur)
        frame_line(&f, "%s%s %s %s%s%s%s",
                   yame_ui_cyan(), yame_ui_unicode() ? "❯" : ">",
                   mark, yame_ui_bold(), label, yame_ui_reset(),
                   notes && notes[i] ? "" : "");
      else
        frame_line(&f, "  %s %s%s%s", mark,
                   yame_ui_reset(), label, yame_ui_reset());
    }

    /* Footer: what the keys do, and how much is selected. */
    if (multi)
      frame_line(&f, "%s  %zu/%zu shown  %s  %zu selected  %s  "
                     "arrows move  space toggles  a/n all/none  / search  "
                     "enter accept  esc cancel%s",
                 yame_ui_dim(), nview, n, yame_ui_bullet(), chosen,
                 yame_ui_bullet(), yame_ui_reset());
    else
      frame_line(&f, "%s  %zu/%zu shown  %s  arrows move  / search  "
                     "enter select  esc cancel%s",
                 yame_ui_dim(), nview, n, yame_ui_bullet(), yame_ui_reset());

    frame_finish(&f);
    fflush(stderr);

    char ch = 0;
    keycode_t k = read_key(&ch);
    if (interrupted) { result = -1; break; }

    if (filtering) {
      /* While filtering, printable keys extend the pattern. */
      switch (k) {
      case K_CHAR: {
        size_t l = strlen(filter);
        if (l + 1 < sizeof(filter)) { filter[l] = ch; filter[l+1] = '\0'; }
        cur = top = 0;
        continue;
      }
      case K_BACKSPACE: {
        size_t l = strlen(filter);
        if (l) filter[l-1] = '\0';
        cur = top = 0;
        continue;
      }
      case K_ESC:   filter[0] = '\0'; filtering = 0; cur = top = 0; continue;
      case K_ENTER: filtering = 0; continue;
      default: break;      /* arrows still navigate while filtering */
      }
    }

    switch (k) {
    case K_UP:   if (cur) --cur; break;
    case K_DOWN: if (cur + 1 < nview) ++cur; break;
    case K_PGUP: cur = (cur > (size_t)avail) ? cur - (size_t)avail : 0; break;
    case K_PGDN: cur += (size_t)avail; if (cur >= nview) cur = nview ? nview-1 : 0; break;
    case K_HOME: cur = 0; break;
    case K_END:  cur = nview ? nview - 1 : 0; break;

    case K_SPACE:
      if (multi && nview) { size_t i = view[cur]; flags[i] = !flags[i]; }
      if (cur + 1 < nview) ++cur;
      break;

    case K_ENTER:
      if (!nview) break;
      result = (long)view[cur];
      done = 1;
      break;

    case K_ESC:
      result = -1;
      done = 1;
      break;

    case K_CHAR:
      if (ch == 'j') { if (cur + 1 < nview) ++cur; }
      else if (ch == 'k') { if (cur) --cur; }
      else if (ch == '/') { filtering = 1; }
      else if (ch == 'q') { result = -1; done = 1; }
      else if (multi && ch == 'a') { for (size_t i = 0; i < n; ++i) flags[i] = 1; }
      else if (multi && ch == 'n') { for (size_t i = 0; i < n; ++i) flags[i] = 0; }
      else if (multi && ch == ' ') { /* handled above */ }
      break;

    case K_NONE:
      result = -1;
      done = 1;
      break;

    default: break;
    }
  }

  /* raw_leave() drops the alternate screen, so the widget vanishes and the
   * caller's summary lands on the user's original prompt line. */
  free(view);
  raw_leave();
  return result;
}

long yame_ui_choose(const char *title, const char **items, const char **notes,
                    size_t n) {
  if (!n) return -1;

  long r = select_widget(title, items, notes, n, NULL, 0);
  if (r != -2) {
    if (r >= 0)
      fprintf(stderr, "%s%s%s %s%s%s\n", yame_ui_green(), yame_ui_check(),
              yame_ui_reset(), yame_ui_bold(), items[r], yame_ui_reset());
    return r;
  }

  /* Fallback: no usable terminal, so ask in one line. */
  fprintf(stderr, "\n%s%s%s\n", yame_ui_bold(), title, yame_ui_reset());
  print_items(items, notes, n);

  for (;;) {
    fprintf(stderr, "  %sselect 1-%zu%s ", yame_ui_dim(), n, yame_ui_reset());
    fflush(stderr);

    char *ans = read_line();
    if (!ans) { fputc('\n', stderr); return -1; }

    char *end = NULL;
    long v = strtol(ans, &end, 10);
    int clean = (end && end != ans && *end == '\0');
    free(ans);

    if (clean && v >= 1 && v <= (long)n) return v - 1;
    fprintf(stderr, "  %sEnter a number between 1 and %zu.%s\n",
            yame_ui_yellow(), n, yame_ui_reset());
  }
}

/**
 * Parse "all" / "none" / "1-5,8,12" into the flag array.
 * Returns 0 on success, -1 if the spec is malformed.
 */
static int parse_selection(const char *spec, int *flags, size_t n) {
  if (strcasecmp(spec, "all") == 0 || strcmp(spec, "*") == 0) {
    for (size_t i = 0; i < n; ++i) flags[i] = 1;
    return 0;
  }
  if (strcasecmp(spec, "none") == 0) {
    for (size_t i = 0; i < n; ++i) flags[i] = 0;
    return 0;
  }

  for (size_t i = 0; i < n; ++i) flags[i] = 0;

  const char *p = spec;
  while (*p) {
    while (*p == ',' || isspace((unsigned char)*p)) ++p;
    if (!*p) break;

    char *end = NULL;
    long lo = strtol(p, &end, 10);
    if (end == p) return -1;
    long hi = lo;
    p = end;

    if (*p == '-') {
      ++p;
      hi = strtol(p, &end, 10);
      if (end == p) return -1;
      p = end;
    }
    if (lo > hi) { long t = lo; lo = hi; hi = t; }
    if (lo < 1 || hi > (long)n) return -1;

    for (long i = lo; i <= hi; ++i) flags[i-1] = 1;

    while (isspace((unsigned char)*p)) ++p;
    if (*p && *p != ',') return -1;
  }
  return 0;
}

int *yame_ui_multiselect(const char *title, const char **items,
                         const char **notes, size_t n, int default_all) {
  if (!n) return NULL;

  int *flags = calloc(n, sizeof(int));
  if (!flags) return NULL;
  if (default_all) for (size_t i = 0; i < n; ++i) flags[i] = 1;

  long r = select_widget(title, items, notes, n, flags, 1);
  if (r != -2) {
    if (r < 0) { free(flags); return NULL; }
    size_t chosen = 0;
    for (size_t i = 0; i < n; ++i) chosen += (size_t)(flags[i] != 0);
    if (!chosen) { free(flags); return NULL; }
    fprintf(stderr, "%s%s%s %s%zu selected%s\n", yame_ui_green(),
            yame_ui_check(), yame_ui_reset(), yame_ui_bold(), chosen,
            yame_ui_reset());
    return flags;
  }

  /* Fallback: no usable terminal, so ask in one line. */
  for (size_t i = 0; i < n; ++i) flags[i] = 0;
  fprintf(stderr, "\n%s%s%s\n", yame_ui_bold(), title, yame_ui_reset());
  print_items(items, notes, n);

  for (;;) {
    fprintf(stderr,
            "  %sselect: 'all', 'none', or a list like 1-5,8  [%s]%s ",
            yame_ui_dim(), default_all ? "all" : "none", yame_ui_reset());
    fflush(stderr);

    char *ans = read_line();
    if (!ans) { fputc('\n', stderr); free(flags); return NULL; }

    const char *spec = *ans ? ans : (default_all ? "all" : "none");
    int rc = parse_selection(spec, flags, n);
    free(ans);

    if (rc == 0) {
      size_t chosen = 0;
      for (size_t i = 0; i < n; ++i) chosen += (size_t)(flags[i] != 0);
      if (chosen) return flags;
      fprintf(stderr, "  %sNothing selected.%s\n",
              yame_ui_yellow(), yame_ui_reset());
      continue;
    }
    fprintf(stderr, "  %sCould not read that. Use 'all', 'none', or "
                    "numbers between 1 and %zu.%s\n",
            yame_ui_yellow(), n, yame_ui_reset());
  }
}

/* Colour for a row's style. Applied around already-measured text, so it can
 * never disturb the column arithmetic. */
static const char *style_color(const unsigned char *styles, size_t i) {
  if (!styles) return "";
  switch (styles[i]) {
  case YAME_ROW_HAVE:     return yame_ui_green();
  case YAME_ROW_MISSING:  return yame_ui_dim();
  case YAME_ROW_REQUIRED: return yame_ui_red();
  case YAME_ROW_PARTIAL:  return yame_ui_yellow();
  default:               return "";
  }
}

/* ---------------------------------------------------------- the browser */

/**
 * Scrolling viewer for tab-separated rows, with the header pinned.
 *
 * Column widths are measured over the whole table once, so columns do not
 * jitter as the viewport scrolls -- the thing that makes a naive pager
 * unreadable for tabular data.
 */
int yame_ui_browse(const char *title, const char *header,
                   const char **rows, const unsigned char *styles, size_t n) {
  if (!n || !yame_ui_fancy()) return -1;
  if (raw_enter() != 0) return -1;

  /* Measure: at most 8 columns, width = widest cell, capped. */
  enum { MAXCOL = 8 };
  int w[MAXCOL] = {0};
  int ncol = 0;

  const char *all_first = header ? header : rows[0];
  for (const char *p = all_first; ; ) {
    const char *tab = strchr(p, '\t');
    if (ncol < MAXCOL) ++ncol;
    if (!tab) break;
    p = tab + 1;
  }

  for (size_t r = 0; r <= n; ++r) {
    const char *line = (r == 0) ? (header ? header : rows[0])
                                : rows[r - 1];
    if (r == 0 && !header) continue;
    int c = 0;
    for (const char *p = line; c < ncol; ++c) {
      const char *tab = strchr(p, '\t');
      int len = (int)(tab ? (size_t)(tab - p) : strlen(p));
      if (len > w[c]) w[c] = len;
      if (!tab) break;
      p = tab + 1;
    }
  }
  for (int c = 0; c < ncol; ++c) { if (w[c] > 34) w[c] = 34; w[c] += 2; }

  /*
   * The viewport tracks a cursor rather than a scroll offset.
   *
   * Scrolling alone is not enough: whenever the table fits on screen -- nine
   * targets in a `fetch` browser, or a 44-row listing in a tall terminal -- there is
   * nothing to scroll, so every arrow key becomes a no-op while the footer
   * still advertises them. A cursor always moves, and the view follows it only
   * when it has to.
   */
  size_t top = 0, cur = 0;
  frame_t f = {0};
  char filter[128] = {0};
  int filtering = 0;

  size_t *view = malloc(n * sizeof(size_t));
  if (!view) { raw_leave(); return -1; }

  for (;;) {
    size_t nview = 0;
    for (size_t i = 0; i < n; ++i)
      if (!filter[0] || strcasestr(rows[i], filter)) view[nview++] = i;

    if (cur >= nview) cur = nview ? nview - 1 : 0;

    int rowsz = term_rows();
    int cols = term_cols();
    int avail = rowsz - 4;
    if (avail < 3) avail = 3;

    /* Keep the cursor in view, scrolling only as far as needed. */
    if (avail > 0) {
      if (cur < top) top = cur;
      if (cur >= top + (size_t)avail) top = cur - (size_t)avail + 1;
      if (top + (size_t)avail > nview)
        top = nview > (size_t)avail ? nview - (size_t)avail : 0;
    }

    frame_begin(&f);
    frame_line(&f, "%s%s%s", yame_ui_bold(), title, yame_ui_reset());

    if (header) {
      char buf[1024]; size_t o = 0;
      int c = 0;
      for (const char *p = header; c < ncol && o + 40 < sizeof(buf); ++c) {
        const char *tab = strchr(p, '\t');
        int len = (int)(tab ? (size_t)(tab - p) : strlen(p));
        if (len > w[c] - 2) len = w[c] - 2;
        o += (size_t)snprintf(buf + o, sizeof(buf) - o, "%-*.*s", w[c], len, p);
        if (!tab) break;
        p = tab + 1;
      }
      char cut[1024];
      fit(buf, cols - 4, cut, sizeof(cut));
      /* Indented by the width of the cursor marker, so the header sits over
       * the columns it names rather than two cells to their left. */
      frame_line(&f, "  %s%s%s", yame_ui_dim(), cut, yame_ui_reset());
    } else {
      frame_line(&f, "");
    }

    for (int k = 0; k < avail; ++k) {
      size_t vi = top + (size_t)k;
      if (vi >= nview) { frame_line(&f, ""); continue; }

      char buf[1024]; size_t o = 0;
      int c = 0;
      for (const char *p = rows[view[vi]]; c < ncol && o + 40 < sizeof(buf); ++c) {
        const char *tab = strchr(p, '\t');
        int len = (int)(tab ? (size_t)(tab - p) : strlen(p));
        if (len > w[c] - 2) len = w[c] - 2;
        o += (size_t)snprintf(buf + o, sizeof(buf) - o, "%-*.*s", w[c], len, p);
        if (!tab) break;
        p = tab + 1;
      }
      char cut[1024];
      fit(buf, cols - 4, cut, sizeof(cut));

      const char *col = style_color(styles, view[vi]);
      if (vi == cur)
        frame_line(&f, "%s%s%s %s%s%s%s", yame_ui_cyan(),
                   yame_ui_unicode() ? "❯" : ">", yame_ui_reset(),
                   yame_ui_bold(), col, cut, yame_ui_reset());
      else
        frame_line(&f, "  %s%s%s", col, cut, yame_ui_reset());
    }

    /* Say "scroll" only when there is something off screen to scroll to. */
    const char *motion = ((size_t)avail < nview) ? "arrows scroll" : "arrows move";

    if (filtering || filter[0])
      frame_line(&f, "%s  search: %s%s%s%s   %zu/%zu   %s  "
                     "esc clear  q quit%s",
                 yame_ui_dim(), yame_ui_cyan(), filter,
                 filtering ? "_" : "", yame_ui_dim(), nview, n, motion,
                 yame_ui_reset());
    else
      frame_line(&f, "%s  row %zu of %zu   %s  / search  q quit%s",
                 yame_ui_dim(), nview ? cur + 1 : 0, nview, motion,
                 yame_ui_reset());

    frame_finish(&f);
    fflush(stderr);

    char ch = 0;
    keycode_t k = read_key(&ch);
    if (interrupted) break;

    if (filtering) {
      if (k == K_CHAR) {
        size_t l = strlen(filter);
        if (l + 1 < sizeof(filter)) { filter[l] = ch; filter[l+1] = '\0'; }
        cur = top = 0; continue;
      }
      if (k == K_BACKSPACE) {
        size_t l = strlen(filter);
        if (l) filter[l-1] = '\0';
        cur = top = 0; continue;
      }
      if (k == K_ENTER) { filtering = 0; continue; }
      if (k == K_ESC)   { filter[0] = '\0'; filtering = 0; cur = top = 0; continue; }
    }

    if (k == K_DOWN) { if (cur + 1 < nview) ++cur; }
    else if (k == K_UP) { if (cur) --cur; }
    else if (k == K_PGDN) { cur += (size_t)avail; if (cur >= nview) cur = nview ? nview - 1 : 0; }
    else if (k == K_PGUP) { cur = (cur > (size_t)avail) ? cur - (size_t)avail : 0; }
    else if (k == K_HOME) cur = 0;
    else if (k == K_END)  cur = nview ? nview - 1 : 0;
    else if (k == K_ESC) { if (filter[0]) { filter[0] = '\0'; cur = top = 0; } else break; }
    else if (k == K_NONE) break;
    else if (k == K_CHAR) {
      if (ch == 'q') break;
      else if (ch == 'j') { if (cur + 1 < nview) ++cur; }
      else if (ch == 'k') { if (cur) --cur; }
      else if (ch == '/') filtering = 1;
    }
  }

  free(view);
  raw_leave();
  return 0;
}

/* ------------------------------------------------------------- the tree */

/**
 * A node that owns nodes.
 *
 * This was two levels by construction -- an array of roots, each holding one
 * flat list of children -- which put a ceiling on what a catalogue could say.
 * A source, its platforms, and each platform's knowledgebase are three levels
 * of nesting, and they had to be flattened into one row per directory:
 * twenty strings like "InfiniumAnnotation/EPIC/KYCG" in a wall, with the
 * shared prefix repeated on every line and nothing to fold away. Nesting is
 * what makes a catalogue navigable rather than merely listed.
 *
 * Children arrive from expand() one level at a time and are kept, so opening
 * a row twice costs one call.
 */
typedef struct tnode_s {
  char *row;              /* display text, owned */
  char *key;              /* the caller's opaque key, owned; NULL at a root */
  char *path;             /* '/'-joined keys from the root down, owned */
  unsigned char style;
  int    checked;         /* only ever set on a leaf */
  int    expanded;
  int    loaded;          /* expand() has been asked already */
  int    branch;          /* holds children -- a claim until loaded, then fact */
  int    depth;           /* roots are 0; the unrendered forest head is -1 */
  struct tnode_s  *parent;
  struct tnode_s **kids;
  size_t n_kids;
} tnode_t;

/**
 * A root that cannot open is a heading: a label for the rows beneath it.
 *
 * It gets neither a fold marker nor a checkbox -- it is not a thing to open
 * or take -- so it is drawn flush left instead of indented past columns that
 * say nothing about it, and the cursor steps over it rather than landing on a
 * row where every key is a no-op.
 */
static int tn_is_heading(const tnode_t *n) {
  return n->depth == 0 && !n->branch;
}

/**
 * Move `cur` off a heading, preferring `dir` (+1 down, -1 up).
 *
 * A heading is a label: landing on it means every key does nothing, and
 * arrowing through the list would stop three times for no reason. Falls back
 * to the other direction at the ends of the list.
 */
static size_t cursor_skip(tnode_t **flat, size_t nflat, size_t cur, int dir) {
  if (!nflat) return 0;
  if (cur >= nflat) cur = nflat - 1;
  if (!tn_is_heading(flat[cur])) return cur;

  size_t c = cur;
  while (tn_is_heading(flat[c])) {
    if (dir >= 0) { if (c + 1 >= nflat) break; ++c; }
    else          { if (c == 0) break; --c; }
  }
  if (!tn_is_heading(flat[c])) return c;

  c = cur;                            /* ran out that way; try the other */
  while (tn_is_heading(flat[c])) {
    if (dir >= 0) { if (c == 0) break; --c; }
    else          { if (c + 1 >= nflat) break; ++c; }
  }
  return c;
}

/** Index of the action bound to `ch`, or -1. */
static int find_action(const yame_ui_tree_t *spec, char ch) {
  for (size_t i = 0; i < spec->n_actions; ++i)
    if (spec->actions[i].key == ch) return (int)i;
  return -1;
}

/* Same as style_color, for a node carrying its own style rather than sitting
 * in an array of them. */
static const char *style_color1(unsigned char s) { return style_color(&s, 0); }

static int tn_is_last(const tnode_t *n) {
  const tnode_t *p = n->parent;
  return p && p->n_kids && p->kids[p->n_kids - 1] == n;
}

/* Before loading, the caller's claim; after, what it turned out to be. */
static int tn_is_leaf(const tnode_t *n) {
  return n->loaded ? (n->n_kids == 0) : !n->branch;
}

static char *tn_join(const char *parent_path, const char *key) {
  if (!key) key = "";
  if (!parent_path || !*parent_path) return strdup(key);
  size_t n = strlen(parent_path) + strlen(key) + 2;
  char *p = malloc(n);
  if (p) snprintf(p, n, "%s/%s", parent_path, key);
  return p;
}

/* A root's path segment is its first tab-separated field: the rest of the row
 * is presentation, and a path has to be something the caller can parse back. */
static char *tn_root_path(const char *row) {
  const char *tab = strchr(row, '\t');
  size_t len = tab ? (size_t)(tab - row) : strlen(row);
  char *p = malloc(len + 1);
  if (p) { memcpy(p, row, len); p[len] = '\0'; }
  return p;
}

static void tn_free_kids(tnode_t *n) {
  for (size_t i = 0; i < n->n_kids; ++i) {
    tn_free_kids(n->kids[i]);
    free(n->kids[i]->row);
    free(n->kids[i]->key);
    free(n->kids[i]->path);
    free(n->kids[i]);
  }
  free(n->kids);
  n->kids = NULL;
  n->n_kids = 0;
  n->loaded = 0;
}

/**
 * Ask the caller for one node's children.
 *
 * The strings in yame_ui_kids_t are taken over rather than copied -- the
 * widget frees them, which is the contract the expand callback is already
 * written against.
 */
static void tn_load(tnode_t *n, const yame_ui_tree_t *spec) {
  if (n->loaded) return;
  n->loaded = 1;
  if (!spec->expand) { n->branch = 0; return; }

  yame_ui_kids_t k = {0};
  spec->expand(spec->ctx, n->path, &k);

  size_t taken = 0;
  if (k.n && (n->kids = calloc(k.n, sizeof(tnode_t *))) != NULL) {
    for (; taken < k.n; ++taken) {
      tnode_t *c = calloc(1, sizeof(tnode_t));
      if (!c) break;
      c->row    = k.rows ? k.rows[taken] : NULL;
      c->key    = k.keys ? k.keys[taken] : NULL;
      c->style  = k.styles ? k.styles[taken] : (unsigned char)YAME_ROW_PLAIN;
      c->branch = (k.branch && k.branch[taken]) ? 1 : 0;
      c->depth  = n->depth + 1;
      c->parent = n;
      c->path   = tn_join(n->path, c->key ? c->key : c->row);
      if (!c->row) c->row = strdup("");
      n->kids[n->n_kids++] = c;
    }
  }
  /* Whatever could not be taken over is still ours to release. */
  for (size_t j = taken; j < k.n; ++j) {
    if (k.rows) free(k.rows[j]);
    if (k.keys) free(k.keys[j]);
  }
  free(k.rows); free(k.keys); free(k.styles); free(k.branch);
  n->branch = (n->n_kids > 0);
}

static void tn_load_deep(tnode_t *n, const yame_ui_tree_t *spec) {
  tn_load(n, spec);
  for (size_t i = 0; i < n->n_kids; ++i) tn_load_deep(n->kids[i], spec);
}

/* Present or mandatory: either way there is nothing to ask for, so it is not
 * offered as a choice. */
static int tn_locked(const tnode_t *n, const yame_ui_tree_t *spec) {
  if (n->style == YAME_ROW_REQUIRED) return 1;
  if (!spec->have_selectable && n->style == YAME_ROW_HAVE) return 1;
  return 0;
}

/* Only leaves are checkable. Checking a branch means checking everything
 * under it, which keeps the accept walk a matter of collecting leaves and
 * spares every caller a rule about what a half-selected directory means. */
static int tn_selectable(const tnode_t *n, const yame_ui_tree_t *spec) {
  return tn_is_leaf(n) && n->key && !tn_locked(n, spec);
}

static size_t tn_count_checked(const tnode_t *n) {
  size_t c = n->checked ? 1 : 0;
  for (size_t i = 0; i < n->n_kids; ++i) c += tn_count_checked(n->kids[i]);
  return c;
}

static void tn_set_checked(tnode_t *n, const yame_ui_tree_t *spec, int val) {
  if (tn_selectable(n, spec)) n->checked = val;
  for (size_t i = 0; i < n->n_kids; ++i) tn_set_checked(n->kids[i], spec, val);
}

static void tn_apply_pred(tnode_t *n, const yame_ui_tree_t *spec,
                          yame_ui_preselect_fn pred) {
  /* Normalize to 0/1: checked is a flag, and a predicate returning, say, 256
   * to mean yes would truncate to 0. */
  if (tn_selectable(n, spec))
    n->checked = pred(spec->ctx, n->path, n->key) ? 1 : 0;
  for (size_t i = 0; i < n->n_kids; ++i) tn_apply_pred(n->kids[i], spec, pred);
}

static size_t tn_count_selectable(const tnode_t *n,
                                  const yame_ui_tree_t *spec) {
  size_t c = tn_selectable(n, spec) ? 1 : 0;
  for (size_t i = 0; i < n->n_kids; ++i)
    c += tn_count_selectable(n->kids[i], spec);
  return c;
}

/**
 * The four-cell marker at the head of a row.
 *
 * Four cells whatever it says, so the text after it lines up and the check
 * sits where the x would be in "[x]". A directory gets a box of its own
 * because space takes everything under it -- without one, that a folder is
 * itself selectable is something you have to be told rather than see. Partly
 * selected is "[*]": a folder showing "[x]" with two of its thirty sets
 * checked would be a lie, and "[ ]" would be a different one.
 */
static void node_box(const tnode_t *n, const yame_ui_tree_t *spec,
                     char *out, size_t cap) {
  int uni = yame_ui_unicode();

  if (!tn_is_leaf(n)) {
    /**
     * A folder's marker answers "how much of this do I have", then "how much
     * of the rest have I asked for" -- in that order, because the first is
     * true whether or not anyone has touched the row, and the colour beside
     * it says which of the two the [*] is about.
     *
     * Both readings come from the caller's counts, so they hold before the
     * row has ever been opened.
     */
    int all_here  = (n->style == YAME_ROW_HAVE && !spec->have_selectable);
    int some_here = (n->style == YAME_ROW_PARTIAL);

    if (all_here) { snprintf(out, cap, "%s", uni ? " ✓  " : " ok "); return; }
    if (!n->loaded) {
      snprintf(out, cap, "%s", some_here ? "[*] " : "[ ] ");
      return;
    }

    size_t sel = tn_count_selectable(n, spec);
    size_t chk = tn_count_checked(n);
    if (!sel)                  snprintf(out, cap, "%s",
                                        some_here ? "[*] " : "    ");
    else if (chk >= sel)       snprintf(out, cap, "[x] ");
    else if (chk || some_here) snprintf(out, cap, "[*] ");
    else                       snprintf(out, cap, "[ ] ");
    return;
  }

  if (n->style == YAME_ROW_REQUIRED)
    snprintf(out, cap, "%s", uni ? " ●  " : " ** ");
  else if (!spec->have_selectable && n->style == YAME_ROW_HAVE)
    snprintf(out, cap, "%s", uni ? " ✓  " : " ok ");
  else if (n->key)
    snprintf(out, cap, "%s", n->checked ? "[x] " : "[ ] ");
  else
    snprintf(out, cap, "    ");
}

static void tn_accept_walk(tnode_t *n, const yame_ui_action_t *act, void *ctx) {
  if (n->checked && n->key && act->accept) act->accept(ctx, n->path, n->key);
  for (size_t i = 0; i < n->n_kids; ++i) tn_accept_walk(n->kids[i], act, ctx);
}

/* Open every ancestor of a checked leaf, so a preselection is visible rather
 * than merely counted in the footer. Returns 1 when the subtree holds one. */
static int tn_open_to_checked(tnode_t *n) {
  int any = n->checked ? 1 : 0;
  for (size_t i = 0; i < n->n_kids; ++i)
    if (tn_open_to_checked(n->kids[i])) any = 1;
  if (any && n->n_kids && n->depth >= 0) n->expanded = 1;
  return any;
}

/**
 * The indent for one node: a stem for itself, and for each ancestor a rule
 * when that ancestor still has siblings to come.
 *
 * Without the rules a three-level tree reads as a stack of indented lines
 * with no way to tell what belongs to what -- which is most of the value of
 * nesting the catalogue in the first place. Every segment is four cells, so
 * the width a caller has left for text is 4 * depth.
 */
static void tn_prefix(const tnode_t *n, char *out, size_t cap) {
  const tnode_t *chain[24];
  size_t d = 0;
  for (const tnode_t *p = n; p && p->depth >= 1 && d < 24; p = p->parent)
    chain[d++] = p;

  int uni = yame_ui_unicode();
  size_t o = 0;
  out[0] = '\0';
  /* chain[d-1] is the shallowest ancestor; chain[0] is the node itself. */
  for (size_t i = d; i-- > 0; ) {
    const char *seg;
    if (i == 0)
      seg = tn_is_last(chain[i]) ? (uni ? "└── " : "`-- ")
                                 : (uni ? "├── " : "|-- ");
    else
      seg = tn_is_last(chain[i]) ? "    " : (uni ? "│   " : "|   ");
    size_t l = strlen(seg);
    if (o + l + 1 >= cap) break;
    memcpy(out + o, seg, l);
    o += l;
    out[o] = '\0';
  }
}

/* ---- search ----
 *
 * Reading order, ignoring folds: a match inside a closed node is still a
 * match, and landing on one opens whatever it takes to see it. A search that
 * found only what was already on screen would be a search for nothing.
 */

static int flat_push(tnode_t ***arr, size_t *n, size_t *cap, tnode_t *t);

static void tn_collect_matches(tnode_t *n, const char *q, tnode_t ***hits,
                               size_t *nh, size_t *cap) {
  if (n->depth >= 0 && n->row && strcasestr(n->row, q))
    flat_push(hits, nh, cap, n);
  for (size_t i = 0; i < n->n_kids; ++i)
    tn_collect_matches(n->kids[i], q, hits, nh, cap);
}

/* Open everything above a node, so the cursor can be put on it. */
static void tn_reveal(tnode_t *n) {
  for (tnode_t *p = n->parent; p && p->depth >= 0; p = p->parent)
    p->expanded = 1;
}

/* Recompute the hit list and aim the cursor at hit `*hit_i`. */
static void search_apply(tnode_t *forest, const char *q, tnode_t ***hits,
                         size_t *nh, size_t *cap, size_t *hit_i,
                         tnode_t **want_cur) {
  *nh = 0;
  if (!*q) return;
  tn_collect_matches(forest, q, hits, nh, cap);
  if (!*nh) return;
  if (*hit_i >= *nh) *hit_i = 0;
  tn_reveal((*hits)[*hit_i]);
  *want_cur = (*hits)[*hit_i];
}

/**
 * The key cheatsheet, as its own screen.
 *
 * The footer used to carry every binding, which made it the longest line on
 * screen and still incomplete -- it had no room to say what a marker or a
 * colour meant. One line that says where the rest is beats a wrapped line
 * that says most of it.
 */
/* One cheatsheet line: keys, then what they do, in a fixed column measured
 * in cells -- "↑ ↓   j k" is nine of them and thirteen bytes. */
static void help_row(frame_t *f, const char *keys, const char *what) {
  enum { KEYW = 20 };
  char pad[KEYW + 1];
  int n = KEYW - cells_of(keys);
  if (n < 1) n = 1;
  memset(pad, ' ', (size_t)n);
  pad[n] = '\0';
  frame_line(f, "    %s%s%s", keys, pad, what);
}

static void draw_help(const yame_ui_tree_t *spec, int picking) {
  frame_t f = {0};
  int uni = yame_ui_unicode();
  const char *dim = yame_ui_dim(), *rst = yame_ui_reset();
  const char *cyan = yame_ui_cyan(), *bold = yame_ui_bold();
  char buf[256];

  frame_begin(&f);
  frame_line(&f, "%s%s%s", bold, spec->title ? spec->title : "keys", rst);
  frame_line(&f, "");

  frame_line(&f, "  %sMOVING%s", cyan, rst);
  help_row(&f, uni ? "↑ ↓   j k" : "up down   j k", "move");
  help_row(&f, uni ? "→   enter   l" : "right   enter   l", "open");
  help_row(&f, uni ? "←" : "left",
           "close, or close the parent when already closed");
  help_row(&f, "PgUp PgDn", "page");
  help_row(&f, "Home End", "first, last");
  frame_line(&f, "");

  if (picking) {
    frame_line(&f, "  %sCHOOSING%s", cyan, rst);
    help_row(&f, "space   x", "select the row, or everything under a folder");
    help_row(&f, "a", "select everything, or clear the selection");
    if (spec->recommend)
      help_row(&f, "r", "the recommended selection for what the cursor is on");
    for (size_t a = 0; a < spec->n_actions; ++a) {
      char k[2] = { spec->actions[a].key, '\0' };
      snprintf(buf, sizeof(buf), "%s what is selected",
               spec->actions[a].verb ? spec->actions[a].verb : "use");
      help_row(&f, k, buf);
    }
    frame_line(&f, "");

    frame_line(&f, "  %sMARKERS%s", cyan, rst);
    help_row(&f, "[ ]", "none of it here or chosen");
    help_row(&f, "[*]", "some of it here, or some of it chosen");
    help_row(&f, "[x]", "everything still missing is chosen");
    help_row(&f, uni ? "✓" : "ok", "all of it already in the store");
    help_row(&f, uni ? "●" : "**", "comes with any fetch, not a choice");
    snprintf(buf, sizeof(buf), "%sall%s / %ssome%s / %snone%s of it in the store",
             yame_ui_green(), rst, yame_ui_yellow(), rst, dim, rst);
    help_row(&f, "colour", buf);
    frame_line(&f, "");
  }

  frame_line(&f, "  %sFINDING%s", cyan, rst);
  help_row(&f, "/", "search everything, opened or not");
  help_row(&f, "n   N", "next / previous match");
  frame_line(&f, "");

  frame_line(&f, "  %sOTHER%s", cyan, rst);
  if (spec->detail && spec->detail_key) {
    char k[2] = { spec->detail_key, '\0' };
    snprintf(buf, sizeof(buf), "%s pane",
             spec->detail_verb ? spec->detail_verb : "info");
    help_row(&f, k, buf);
  }
  help_row(&f, "h   ?", "this list");
  help_row(&f, "q", "quit");
  if (spec->hint) {
    frame_line(&f, "");
    frame_line(&f, "  %s%s%s", dim, spec->hint, rst);
  }
  frame_line(&f, "");
  frame_line(&f, "%s  any key to go back%s", dim, rst);
  frame_finish(&f);
  fflush(stderr);
}

/**
 * The width the row layout uses, which is not always the terminal's.
 *
 * A field right-aligned at the margin of a 200-column terminal sits half a
 * screen from the name it belongs to, and nothing pairs them by eye. Past
 * this the layout simply stops growing and the rest of the line is empty
 * space, which is what a wide terminal is for.
 */
enum { LAYOUT_MAX = 92 };

static int layout_cols(int cols) {
  return cols < LAYOUT_MAX ? cols : LAYOUT_MAX;
}

/**
 * Build one row's line.
 *
 * The LAST tab-separated field of a row is held back and right-aligned at the
 * margin. That one rule is what lets a folder's gauge, a file's size and a
 * root's column line up in a single column no matter how deep they sit --
 * padding each row to a fixed width could not, since every level of nesting
 * spends four more cells on the indent and a long filename spends the rest.
 *
 * `tail` overrides that field, which is how a transfer paints a progress bar
 * exactly where the size will be once it finishes.
 */
static void build_row(char *out, size_t cap, const tnode_t *nd, int is_cur,
                      int picking, const yame_ui_tree_t *spec, int cols,
                      const int *w, int ncol,
                      const char *tail, const char *tail_color) {
  int uni = yame_ui_unicode();
  const char *arrow = is_cur ? (uni ? "❯" : ">") : " ";
  const char *fold = tn_is_leaf(nd)
                       ? " "
                       : (nd->expanded ? (uni ? "▾" : "-") : (uni ? "▸" : "+"));

  char box[8] = "";
  if (picking) node_box(nd, spec, box, sizeof(box));

  char pre[256] = "";
  if (nd->depth > 0) tn_prefix(nd, pre, sizeof(pre));

  const char *row = nd->row ? nd->row : "";
  const char *lastt = strrchr(row, '\t');
  if (!tail) tail = lastt ? lastt + 1 : "";

  char content[1024];
  if (nd->depth == 0) {
    /* Roots are columnar; the held-back field is not one of the columns. */
    columns(row, w, ncol > 1 ? ncol - 1 : ncol, content, sizeof(content));
  } else {
    size_t n = lastt ? (size_t)(lastt - row) : strlen(row);
    if (n >= sizeof(content)) n = sizeof(content) - 1;
    memcpy(content, row, n);
    content[n] = '\0';
  }

  int heading = tn_is_heading(nd);
  int lead = heading ? 2 : 2 + 4 * nd->depth + 2 + (picking ? 4 : 0);
  int tailc = cells_of(tail);
  cols = layout_cols(cols);
  int room = cols - lead - tailc - 3;
  if (room < 8) room = 8;

  char cut[1024];
  fit(content, room, cut, sizeof(cut));
  int gap = cols - lead - cells_of(cut) - tailc - 2;
  if (gap < 1) gap = 1;

  if (heading)
    snprintf(out, cap, "%s%s%s %s%s%s%*s%s%s%s",
             yame_ui_cyan(), arrow, yame_ui_reset(),
             is_cur ? yame_ui_bold() : "", style_color1(nd->style), cut,
             gap, "", tail_color ? tail_color : yame_ui_dim(), tail,
             yame_ui_reset());
  else
    snprintf(out, cap, "%s%s%s %s%s%s%s%s%s %s%s%s%s%*s%s%s%s",
             yame_ui_cyan(), arrow, yame_ui_reset(),
             yame_ui_dim(), pre, yame_ui_reset(),
             yame_ui_cyan(), fold, yame_ui_reset(),
             is_cur ? yame_ui_bold() : "", style_color1(nd->style), box, cut,
             gap, "", tail_color ? tail_color : yame_ui_dim(), tail,
             yame_ui_reset());
}

/* ---- painting one row while an action runs ----
 *
 * The visible window, republished on every redraw so a commit callback can
 * name a row and have the widget find it. File-static like the panel state
 * above, and for the same reason: the callback runs inside the caller's code,
 * which has no handle on the frame.
 */
static tnode_t **vis_flat = NULL;
static size_t    vis_n = 0, vis_top = 0;
static int       vis_avail = 0;
static const int VIS_ROW0 = 3;   /* screen line of the first list row */
static const yame_ui_tree_t *vis_spec = NULL;
static int       vis_picking = 0;
static const int *vis_w = NULL;
static int       vis_ncol = 0;

int yame_ui_tree_progress(const char *key, uint64_t now, uint64_t total) {
  if (!vis_flat || !key || !raw_active) return 0;

  for (size_t i = vis_top; i < vis_n && i < vis_top + (size_t)vis_avail; ++i) {
    tnode_t *nd = vis_flat[i];
    if (!nd->key || strcmp(nd->key, key) != 0) continue;

    int row = VIS_ROW0 + (int)(i - vis_top);
    int uni = yame_ui_unicode();
    char line[1400];

    if (!total) {
      build_row(line, sizeof(line), nd, 0, vis_picking, vis_spec,
                term_cols(), vis_w, vis_ncol, NULL, NULL);
    } else {
      /* A bar exactly where the size sits, so a transfer reads as the row
       * filling in rather than as a second thing happening elsewhere. */
      enum { BAR = 10 };
      int pct = (int)((now * 100) / total);
      if (pct > 100) pct = 100;
      int on = (pct * BAR + 99) / 100;

      char tail[80];
      size_t o = 0;
      tail[0] = '\0';
      for (int c = 0; c < BAR; ++c) {
        const char *g = (c < on) ? (uni ? "▰" : "#") : (uni ? "▱" : ".");
        size_t l = strlen(g);
        if (o + l + 1 >= sizeof(tail)) break;
        memcpy(tail + o, g, l);
        o += l;
        tail[o] = '\0';
      }
      snprintf(tail + o, sizeof(tail) - o, " %3d%%", pct);

      build_row(line, sizeof(line), nd, 0, vis_picking, vis_spec,
                term_cols(), vis_w, vis_ncol, tail, yame_ui_cyan());
    }

    /* Save the cursor, paint the one line, put it back: the panel below is
     * drawing too, and it positions relative to where it left off. */
    fprintf(stderr, "\033[s\033[%d;1H\033[2K%s\033[u", row, line);
    fflush(stderr);
    return 1;
  }
  return 0;
}

int yame_ui_tree_settle(const char *key, int now_present) {
  if (!vis_flat || !key || !raw_active) return 0;

  for (size_t i = vis_top; i < vis_n && i < vis_top + (size_t)vis_avail; ++i) {
    tnode_t *nd = vis_flat[i];
    if (!nd->key || strcmp(nd->key, key) != 0) continue;

    if (now_present) {
      /* What the reload would have said, said now: it is here, so it is no
       * longer something to ask for. */
      nd->style = YAME_ROW_HAVE;
      nd->checked = 0;
    }

    char line[1400];
    build_row(line, sizeof(line), nd, 0, vis_picking, vis_spec, term_cols(),
              vis_w, vis_ncol, NULL, NULL);
    fprintf(stderr, "\033[s\033[%d;1H\033[2K%s\033[u",
            VIS_ROW0 + (int)(i - vis_top), line);
    fflush(stderr);
    return 1;
  }
  return 0;
}

/* ---- the flattened view ---- */

static int flat_push(tnode_t ***arr, size_t *n, size_t *cap, tnode_t *t) {
  if (*n == *cap) {
    size_t nc = *cap ? *cap * 2 : 64;
    tnode_t **na = realloc(*arr, nc * sizeof(tnode_t *));
    if (!na) return -1;
    *arr = na;
    *cap = nc;
  }
  (*arr)[(*n)++] = t;
  return 0;
}

/* Every node whose ancestors are all open, in reading order. The forest head
 * is structure, not a row, so it is walked through and never listed. */
static int tn_flatten(tnode_t *n, tnode_t ***arr, size_t *cnt, size_t *cap) {
  if (n->depth >= 0) {
    if (flat_push(arr, cnt, cap, n) != 0) return -1;
    if (!n->expanded) return 0;
  }
  for (size_t i = 0; i < n->n_kids; ++i)
    if (tn_flatten(n->kids[i], arr, cnt, cap) != 0) return -1;
  return 0;
}

/* ---- reloading ---- */

static void tn_collect_open(const tnode_t *n, char ***paths, size_t *np,
                            size_t *cap) {
  if (n->depth >= 0 && n->expanded && n->path) {
    if (*np == *cap) {
      size_t nc = *cap ? *cap * 2 : 16;
      char **na = realloc(*paths, nc * sizeof(char *));
      if (!na) return;
      *paths = na;
      *cap = nc;
    }
    char *dup = strdup(n->path);
    if (dup) (*paths)[(*np)++] = dup;
  }
  for (size_t i = 0; i < n->n_kids; ++i)
    tn_collect_open(n->kids[i], paths, np, cap);
}

/* Walk each saved path down from the forest, expanding as it goes. A node can
 * only have been open while its ancestors were, so re-opening along the path
 * restores exactly the set that was saved. */
static void tn_reopen(tnode_t *forest, const yame_ui_tree_t *spec,
                      char *const *paths, size_t np) {
  for (size_t i = 0; i < np; ++i) {
    tnode_t *n = forest;
    for (;;) {
      tn_load(n, spec);
      tnode_t *next = NULL;
      for (size_t j = 0; j < n->n_kids; ++j) {
        const char *cp = n->kids[j]->path;
        size_t l = cp ? strlen(cp) : 0;
        if (l && strncmp(paths[i], cp, l) == 0 &&
            (paths[i][l] == '\0' || paths[i][l] == '/')) {
          next = n->kids[j];
          break;
        }
      }
      if (!next) break;
      tn_load(next, spec);
      if (next->n_kids) next->expanded = 1;
      if (strcmp(next->path, paths[i]) == 0) break;
      n = next;
    }
  }
}

/**
 * Throw away everything the tree holds and ask for it again, keeping the
 * folds.
 *
 * After an action runs, every loaded row describes the store as it was
 * beforehand. Losing your place while re-reading it is the one thing a
 * stay-open commit was meant to avoid, so the fold state survives even though
 * the nodes do not.
 */
static void tn_refresh(tnode_t *forest, const yame_ui_tree_t *spec,
                       char **roots, const unsigned char *root_styles,
                       size_t n_roots, int keep_folds) {
  char **open = NULL;
  size_t nopen = 0, cap = 0;
  if (keep_folds) tn_collect_open(forest, &open, &nopen, &cap);

  for (size_t i = 0; i < forest->n_kids && i < n_roots; ++i) {
    tnode_t *r = forest->kids[i];
    tn_free_kids(r);
    r->expanded = 0;
    r->checked = 0;
    r->branch = !spec->expand ? 0
              : spec->root_branch ? (spec->root_branch[i] != 0) : 1;
    /* A commit may rewrite a root in place to show what it just did. */
    free(r->row);
    free(r->path);
    r->row = strdup(roots[i]);
    r->path = tn_root_path(roots[i]);
    r->style = root_styles ? root_styles[i] : (unsigned char)YAME_ROW_PLAIN;
  }

  tn_reopen(forest, spec, open, nopen);
  for (size_t i = 0; i < nopen; ++i) free(open[i]);
  free(open);
}

int yame_ui_tree(const yame_ui_tree_t *spec) {
  const char *title = spec->title;
  const char *header = spec->header;
  char **roots = spec->roots;
  const unsigned char *root_styles = spec->root_styles;
  size_t n_roots = spec->n_roots;
  yame_ui_key_fn on_key = spec->on_key;
  void *ctx = spec->ctx;

  if (!n_roots || !yame_ui_fancy()) return -1;
  if (raw_enter() != 0) return -1;

  const int picking = (spec->n_actions > 0);
  int accepted = 0;
  /* On by default: the descriptions are the reason the catalogue is browsable
   * rather than just listable, and a key you have to discover is a key most
   * people never press. `detail_key` turns it off. */
  int detail_open = (spec->detail != NULL);
  yame_ui_detail_t det = {0};

  /* Column widths over the root rows only; deeper rows arrive preformatted. */
  enum { MAXCOL = 8 };
  int w[MAXCOL] = {0};
  int ncol = 0;
  for (const char *p = header ? header : roots[0]; ; ) {
    const char *tab = strchr(p, '\t');
    if (ncol < MAXCOL) ++ncol;
    if (!tab) break;
    p = tab + 1;
  }
  for (size_t r = 0; r < n_roots + (header ? 1u : 0u); ++r) {
    const char *line = (header && r == 0) ? header
                                          : roots[header ? r - 1 : r];
    int c = 0;
    for (const char *p = line; c < ncol; ++c) {
      const char *tab = strchr(p, '\t');
      size_t blen = tab ? (size_t)(tab - p) : strlen(p);
      char field[512];
      size_t n = blen < sizeof(field) - 1 ? blen : sizeof(field) - 1;
      memcpy(field, p, n);
      field[n] = '\0';
      int len = cells_of(field);      /* cells: see columns() */
      if (len > w[c]) w[c] = len;
      if (!tab) break;
      p = tab + 1;
    }
  }
  for (int c = 0; c < ncol; ++c) { if (w[c] > 34) w[c] = 34; w[c] += 2; }

  /* An unrendered head above the roots, so every walk below is one recursion
   * rather than a loop over roots wrapping a recursion. */
  tnode_t forest;
  memset(&forest, 0, sizeof(forest));
  forest.depth = -1;
  forest.loaded = 1;
  forest.expanded = 1;
  forest.kids = calloc(n_roots, sizeof(tnode_t *));
  if (!forest.kids) { raw_leave(); return -1; }

  for (size_t i = 0; i < n_roots; ++i) {
    tnode_t *r = calloc(1, sizeof(tnode_t));
    if (!r) break;
    r->row    = strdup(roots[i]);
    r->path   = tn_root_path(roots[i]);
    r->style  = root_styles ? root_styles[i] : (unsigned char)YAME_ROW_PLAIN;
    r->branch = !spec->expand ? 0
              : spec->root_branch ? (spec->root_branch[i] != 0) : 1;
    r->parent = &forest;
    forest.kids[forest.n_kids++] = r;
  }

  size_t cur = 0, top = 0;
  frame_t f = {0};
  tnode_t **flat = NULL;
  size_t flat_cap = 0;
  tnode_t *want_cur = NULL;      /* land here once the view is flattened */

  /* Search state. The hit list holds node pointers, so anything that frees
   * nodes must drop it. */
  int      nav_dir = 1;          /* which way to step off a heading */
  char     query[128] = "";
  int      searching = 0;
  tnode_t *search_home = NULL;   /* where the cursor was when / was pressed */
  tnode_t **hits = NULL;
  size_t   n_hits = 0, hits_cap = 0, hit_i = 0;

  /* A named target arrives open and chosen: the user sees exactly what will
   * happen and presses f, or unpicks first. That is the same screen the
   * catalogue uses, rather than a second confirmation flow beside it. */
  if (spec->open_root) {
    for (size_t i = 0; i < forest.n_kids; ++i) {
      tnode_t *r = forest.kids[i];
      if (!r->path || strcmp(r->path, spec->open_root) != 0) continue;
      tn_load_deep(r, spec);
      if (r->n_kids) r->expanded = 1;
      if (picking && spec->preselect) {
        tn_apply_pred(r, spec, spec->preselect);
        tn_open_to_checked(r);
      }
      want_cur = r;
      break;
    }
  }

  for (;;) {
    size_t nflat = 0;
    if (tn_flatten(&forest, &flat, &nflat, &flat_cap) != 0) break;

    if (want_cur) {
      for (size_t i = 0; i < nflat; ++i)
        if (flat[i] == want_cur) { cur = i; break; }
      want_cur = NULL;
    }
    if (cur >= nflat) cur = nflat ? nflat - 1 : 0;
    cur = cursor_skip(flat, nflat, cur, nav_dir);

    /* Rebuild the detail pane before laying out, so its height is known and
     * the list shrinks to fit rather than the pane running off the screen. */
    for (size_t d = 0; d < det.n; ++d) free(det.rows[d]);
    free(det.rows);
    det.rows = NULL; det.n = 0;

    int rowsz = term_rows();
    int cols = term_cols();
    int avail = rowsz - 4;

    if (detail_open && spec->detail && nflat) {
      tnode_t *dn = flat[cur];
      spec->detail(ctx, dn->path, dn->key, cols - 2, &det);
      /* Never let the pane crowd the list below a few rows: seeing which row
       * the text describes is what makes it useful at all. */
      int cap = rowsz / 2;
      if (cap < 0) cap = 0;
      if ((int)det.n > cap && cap > 0) {
        /* Free what is being dropped first: shrinking det.n is also the loss
         * of the only pointers to those rows, and this runs on every redraw
         * -- so a pane taller than half the screen leaked a line per
         * keystroke. */
        for (size_t d = (size_t)cap; d < det.n; ++d) free(det.rows[d]);
        /* Say that it was cut. Silently stopping mid-sentence reads as the
         * whole of what is recorded, which for provenance is the wrong
         * impression to leave. */
        det.n = (size_t)cap;
        char *mark = malloc(64);
        if (mark) {
          snprintf(mark, 64, "  %s...%s", yame_ui_dim(), yame_ui_reset());
          free(det.rows[det.n - 1]);
          det.rows[det.n - 1] = mark;
        }
      }
      /* No lines, no space: a callback returns nothing for a row it has
       * nothing to say about, and reserving a separator for it would leave a
       * stray rule under the list. */
      if (det.n) avail -= (int)det.n + 1;
    }
    if (avail < 3) avail = 3;

    if (avail > 0) {
      if (cur < top) top = cur;
      if (cur >= top + (size_t)avail) top = cur - (size_t)avail + 1;
      if (top + (size_t)avail > nflat)
        top = nflat > (size_t)avail ? nflat - (size_t)avail : 0;
    }

    /* Republish the window a commit callback may want to paint into. */
    vis_flat = flat; vis_n = nflat; vis_top = top; vis_avail = avail;
    vis_spec = spec; vis_picking = picking; vis_w = w; vis_ncol = ncol;

    frame_begin(&f);
    frame_line(&f, "%s%s%s", yame_ui_bold(), title, yame_ui_reset());

    if (header) {
      /* Same right-aligned last field as the rows, so the header sits over
       * the column it names rather than beside it. */
      char buf[1024];
      columns(header, w, ncol > 1 ? ncol - 1 : ncol, buf, sizeof(buf));
      const char *ht = strrchr(header, '\t');
      ht = ht ? ht + 1 : "";

      int lead = 4 + (picking ? 4 : 0);
      int hcols = layout_cols(cols);
      char cut[1024];
      fit(buf, hcols - lead - cells_of(ht) - 3, cut, sizeof(cut));
      int gap = hcols - lead - cells_of(cut) - cells_of(ht) - 2;
      if (gap < 1) gap = 1;

      frame_line(&f, "    %s%s%s%*s%s%s", picking ? "    " : "",
                 yame_ui_dim(), cut, gap, "", ht, yame_ui_reset());
    } else {
      frame_line(&f, "");
    }

    for (int r = 0; r < avail; ++r) {
      size_t fi = top + (size_t)r;
      if (fi >= nflat) { frame_line(&f, ""); continue; }

      tnode_t *nd = flat[fi];
      char line[1400];
      build_row(line, sizeof(line), nd, fi == cur, picking, spec, cols,
                w, ncol, NULL, NULL);
      frame_line(&f, "%s", line);
    }

    const char *motion = ((size_t)avail < nflat) ? "arrows scroll"
                                                 : "arrows move";

    /* The pane sits between the list and the footer, inside the same frame,
     * so a cursor move redraws both together and never flickers. */
    if (detail_open && det.n) {
      frame_line(&f, "%s  %s%s", yame_ui_dim(),
                 "────────────────────────────────────────", yame_ui_reset());
      for (size_t d = 0; d < det.n; ++d) {
        char cutd[1024];
        fit(det.rows[d], cols - 2, cutd, sizeof(cutd));
        frame_line(&f, "%s", cutd);
      }
    }

    if (searching) {
      /* One line, like every other footer: a search box that added a row
       * would shift the whole list the moment you started typing. */
      char found[64];
      if (!query[0])   snprintf(found, sizeof(found), "type to search");
      else if (n_hits) snprintf(found, sizeof(found), "%zu of %zu",
                                hit_i + 1, n_hits);
      else             snprintf(found, sizeof(found), "no match");

      frame_line(&f, "%s  /%s%s%s%s  %s  %s  %s next  enter keep  esc cancel%s",
                 yame_ui_dim(), yame_ui_bold(), query, yame_ui_reset(),
                 yame_ui_dim(), yame_ui_bullet(), found,
                 yame_ui_unicode() ? "↑↓" : "up/down", yame_ui_reset());
    } else if (picking) {
      size_t nsel = tn_count_checked(&forest);
      char acts[160];
      size_t ao = 0;
      if (spec->detail && spec->detail_key)
        ao += (size_t)snprintf(acts + ao, sizeof(acts) - ao, "%c %s  ",
                               spec->detail_key,
                               detail_open ? "close"
                                           : (spec->detail_verb
                                              ? spec->detail_verb : "info"));
      for (size_t a = 0; a < spec->n_actions && ao + 24 < sizeof(acts); ++a)
        ao += (size_t)snprintf(acts + ao, sizeof(acts) - ao, "%c %s  ",
                               spec->actions[a].key, spec->actions[a].verb);
      /* `acts` already ends in two spaces, so nothing is added here. */
      frame_line(&f, "%s  row %zu of %zu  %s  %zu selected  %s  %sh keys  "
                     "q quit%s",
                 yame_ui_dim(), nflat ? cur + 1 : 0, nflat,
                 yame_ui_bullet(), nsel, yame_ui_bullet(),
                 acts, yame_ui_reset());
    } else {
      frame_line(&f, "%s  row %zu of %zu  %s  %s  %s  h keys  q quit%s",
                 yame_ui_dim(), nflat ? cur + 1 : 0, nflat, motion,
                 yame_ui_bullet(), yame_ui_unicode() ? "→ ← open close"
                                                     : "right left open close",
                 yame_ui_reset());
    }

    frame_finish(&f);
    fflush(stderr);

    char ch = 0;
    keycode_t key = read_key(&ch);
    if (interrupted) break;

    tnode_t *sel = nflat ? flat[cur] : NULL;

    /* While the search box is open it takes every printable key, so typing a
     * name cannot trip over j/k/q/f. */
    if (searching) {
      if (key == K_CHAR) {
        size_t l = strlen(query);
        if (l + 1 < sizeof(query)) { query[l] = ch; query[l + 1] = '\0'; }
        hit_i = 0;
        search_apply(&forest, query, &hits, &n_hits, &hits_cap, &hit_i,
                     &want_cur);
        continue;
      }
      if (key == K_BACKSPACE) {
        size_t l = strlen(query);
        if (l) query[l - 1] = '\0';
        hit_i = 0;
        search_apply(&forest, query, &hits, &n_hits, &hits_cap, &hit_i,
                     &want_cur);
        continue;
      }
      if (key == K_DOWN || key == K_UP) {
        if (n_hits) {
          hit_i = (key == K_DOWN) ? (hit_i + 1) % n_hits
                                  : (hit_i ? hit_i - 1 : n_hits - 1);
          tn_reveal(hits[hit_i]);
          want_cur = hits[hit_i];
        }
        continue;
      }
      if (key == K_ENTER) { searching = 0; continue; }
      if (key == K_ESC) {
        /* Put the cursor back where it was. What the search opened along the
         * way stays open -- collapsing it again would hide a row the user has
         * now seen, and folds are cheap to undo by hand. */
        searching = 0;
        query[0] = '\0';
        n_hits = 0;
        want_cur = search_home;
        continue;
      }
      if (key == K_NONE) break;
      searching = 0;              /* anything else leaves the box and acts */
    }

    int want_open = (key == K_RIGHT || key == K_ENTER ||
                     (key == K_CHAR && ch == 'l'));
    int want_close = (key == K_LEFT);

    /* Opening is lazy: children are asked for the first time a row unfolds. */
    if (want_open && sel) {
      tn_load(sel, spec);
      if (sel->n_kids) sel->expanded = 1;
    } else if (want_close && sel) {
      if (sel->expanded) {
        sel->expanded = 0;
      } else if (sel->parent && sel->parent->depth >= 0) {
        /* Closing from inside folds the parent and jumps to it, the way a
         * file tree behaves -- otherwise the cursor lands on whatever row
         * replaced it. */
        sel->parent->expanded = 0;
        want_cur = sel->parent;
      }
    }
    else if (picking && sel &&
             (key == K_SPACE || (key == K_CHAR && ch == 'x'))) {
      if (tn_is_leaf(sel)) {
        if (tn_selectable(sel, spec)) sel->checked = !sel->checked;
        if (cur + 1 < nflat) ++cur;
      } else {
        /* Space on a branch takes everything under it that is not already
         * present, at whatever depth, opening it so the effect is visible. */
        tn_load_deep(sel, spec);
        if (sel->n_kids) sel->expanded = 1;
        tn_set_checked(sel, spec, tn_count_checked(sel) ? 0 : 1);
      }
    }
    else if (picking && key == K_CHAR && find_action(spec, ch) >= 0) {
      int ai = find_action(spec, ch);
      const yame_ui_action_t *act = &spec->actions[ai];
      if (tn_count_checked(&forest)) {
        tn_accept_walk(&forest, act, ctx);

        if (!act->commit) {    /* caller wants the selection, not the work */
          accepted = ai + 1;
          break;
        }

        /* The panel simply paints over the bottom rows, pane included. The
         * catalogue stays visible above it, which is the part that matters
         * while a download runs; the description underneath is context you
         * have already read. */
        int changed = act->commit(ctx);
        yame_ui_panel_close();
        if (changed) {
          tn_refresh(&forest, spec, roots, root_styles, n_roots, 1);
          /* The hit list points at nodes that no longer exist. */
          n_hits = 0; query[0] = '\0'; search_home = NULL;
        }
        continue;
      }
    }
    else if (key == K_DOWN) { nav_dir = 1; if (cur + 1 < nflat) ++cur; }
    else if (key == K_UP)   { nav_dir = -1; if (cur) --cur; }
    else if (key == K_PGDN) { nav_dir = 1; cur += (size_t)avail; if (cur >= nflat) cur = nflat ? nflat-1 : 0; }
    else if (key == K_PGUP) { nav_dir = -1; cur = (cur > (size_t)avail) ? cur - (size_t)avail : 0; }
    else if (key == K_HOME) { nav_dir = 1; cur = 0; }
    else if (key == K_END)  { nav_dir = -1; cur = nflat ? nflat - 1 : 0; }
    /* ESC closes the pane first, then leaves the widget. */
    else if (key == K_ESC)  { if (detail_open) detail_open = 0; else break; }
    else if (key == K_NONE) break;
    else if (key == K_CHAR) {
      if (ch == 'q') break;
      /* The detail key toggles the pane and nothing else: the cursor keeps
       * moving underneath it, which is the point -- comparing two sets used
       * to mean closing the pane, moving, and reopening it. */
      else if (spec->detail && spec->detail_key && ch == spec->detail_key)
        detail_open = !detail_open;
      else if (picking && find_action(spec, ch) >= 0) { /* above */ }
      else if (picking && ch == 'x') { /* handled above, with space */ }
      else if (ch == 'j') { nav_dir = 1; if (cur + 1 < nflat) ++cur; }
      else if (ch == 'k') { nav_dir = -1; if (cur) --cur; }
      else if (picking && ch == 'a') {
        /* Everything, or nothing. Scoped selection is already what space on a
         * folder does, so the key worth having is the one that reaches rows
         * no cursor is near -- which means loading the tree, as search does.
         * The plan shown by `f` is what keeps that from being dangerous. */
        tn_load_deep(&forest, spec);
        tn_set_checked(&forest, spec, tn_count_checked(&forest) ? 0 : 1);
      }
      else if (ch == 'h' || ch == '?') {
        /* Takes 'h' from vim-style close; the left arrow still does that, and
         * a cheatsheet nobody can find is not one. */
        draw_help(spec, picking);
        char c2 = 0;
        read_key(&c2);
      }
      else if (ch == '/') {
        /* Load the whole tree first. Searching only what happens to be open
         * would miss the thing being looked for nearly every time -- the
         * reason to search a catalogue is that you do not know where in it to
         * look. Expansion is a table scan per node, so this is affordable at
         * catalogue scale; it would not be over a filesystem. */
        tn_load_deep(&forest, spec);
        searching = 1;
        query[0] = '\0';
        n_hits = 0;
        hit_i = 0;
        search_home = sel;
      }
      else if (query[0] && n_hits && (ch == 'n' || ch == 'N')) {
        hit_i = (ch == 'n') ? (hit_i + 1) % n_hits
                            : (hit_i ? hit_i - 1 : n_hits - 1);
        tn_reveal(hits[hit_i]);
        want_cur = hits[hit_i];
      }
      else if (picking && spec->recommend && ch == 'r' && sel) {
        /* Scoped to the collection under the cursor, which for a leaf means
         * the one holding it. Recommending across the whole catalogue would
         * check sets for platforms the user is not working on, and open every
         * collection to do it. */
        tnode_t *scope = (tn_is_leaf(sel) && sel->parent &&
                          sel->parent->depth >= 0) ? sel->parent : sel;
        tn_load_deep(scope, spec);
        if (scope->n_kids) scope->expanded = 1;
        tn_apply_pred(scope, spec, spec->recommend);
        tn_open_to_checked(scope);
      }
      else if (on_key && on_key(ctx, ch, sel ? sel->path : NULL,
                                sel ? sel->key : NULL)) {
        /* The caller changed what the roots describe -- a different store,
         * say -- so every expansion is about something else now, folds
         * included. */
        yame_ui_panel_close();
        tn_refresh(&forest, spec, roots, root_styles, n_roots, 0);
        n_hits = 0; query[0] = '\0'; search_home = NULL;
        cur = top = 0;
      }
    }
  }

  vis_flat = NULL; vis_n = vis_top = 0; vis_avail = 0;
  vis_spec = NULL; vis_w = NULL; vis_ncol = 0;

  tn_free_kids(&forest);
  free(flat);
  free(hits);
  raw_leave();
  for (size_t d = 0; d < det.n; ++d) free(det.rows[d]);
  free(det.rows);

  return accepted;
}
