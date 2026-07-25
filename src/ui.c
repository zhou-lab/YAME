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

typedef struct {
  yame_ui_kids_t kids;
  unsigned char *checked;   /* one per child, when the tree is a picker */
  int expanded;
  int loaded;
} treenode_t;

/* One line of the flattened view: a root, or a child of one. */
typedef struct {
  size_t root;
  long   child;      /* -1 for the root row itself */
} flatrow_t;

/** Index of the action bound to `ch`, or -1. */
static int find_action(const yame_ui_tree_t *spec, char ch) {
  for (size_t i = 0; i < spec->n_actions; ++i)
    if (spec->actions[i].key == ch) return (int)i;
  return -1;
}

int yame_ui_tree(const yame_ui_tree_t *spec) {
  const char *title = spec->title;
  const char *header = spec->header;
  char **roots = spec->roots;
  const unsigned char *root_styles = spec->root_styles;
  size_t n_roots = spec->n_roots;
  yame_ui_expand_fn expand = spec->expand;
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

  /* Column widths over the root rows only; children arrive preformatted. */
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
      int len = (int)(tab ? (size_t)(tab - p) : strlen(p));
      if (len > w[c]) w[c] = len;
      if (!tab) break;
      p = tab + 1;
    }
  }
  for (int c = 0; c < ncol; ++c) { if (w[c] > 34) w[c] = 34; w[c] += 2; }

  treenode_t *node = calloc(n_roots, sizeof(treenode_t));
  flatrow_t *flat = NULL;
  size_t flat_m = 0;
  if (!node) { raw_leave(); return -1; }

  size_t cur = 0, top = 0;
  frame_t f = {0};

  /* A named target arrives open and chosen: the user sees exactly what will
   * happen and presses f, or unpicks first. That is the same screen the
   * catalogue uses, rather than a second confirmation flow beside it. */
  if (spec->open_root) {
    for (size_t i = 0; i < n_roots; ++i) {
      const char *tab = strchr(roots[i], '\t');
      size_t len = tab ? (size_t)(tab - roots[i]) : strlen(roots[i]);
      if (strlen(spec->open_root) != len ||
          strncmp(roots[i], spec->open_root, len) != 0) continue;

      if (expand) expand(ctx, roots[i], &node[i].kids);
      if (picking && node[i].kids.n) node[i].checked = calloc(node[i].kids.n, 1);
      node[i].loaded = 1;
      node[i].expanded = node[i].kids.n ? 1 : 0;

      if (picking && spec->preselect && node[i].checked && node[i].kids.keys) {
        for (size_t j = 0; j < node[i].kids.n; ++j) {
          int req = node[i].kids.styles &&
                    (node[i].kids.styles[j] == YAME_ROW_HAVE ||
                     node[i].kids.styles[j] == YAME_ROW_REQUIRED);
          if (!req && node[i].kids.keys[j] &&
              spec->preselect(ctx, roots[i], node[i].kids.keys[j]))
            node[i].checked[j] = 1;
        }
      }
      /* Land on the row that was named, not on the top of the list. */
      for (size_t k = 0; k < i; ++k)
        cur += 1 + (node[k].expanded ? node[k].kids.n : 0);
      break;
    }
  }

  for (;;) {
    /* Flatten: every root, plus the children of the open ones. */
    size_t nflat = 0;
    for (size_t i = 0; i < n_roots; ++i)
      nflat += 1 + (node[i].expanded ? node[i].kids.n : 0);

    if (nflat > flat_m) {
      flatrow_t *nf = realloc(flat, nflat * sizeof(flatrow_t));
      if (!nf) break;
      flat = nf; flat_m = nflat;
    }

    size_t k = 0;
    for (size_t i = 0; i < n_roots; ++i) {
      flat[k].root = i; flat[k].child = -1; ++k;
      if (!node[i].expanded) continue;
      for (size_t j = 0; j < node[i].kids.n; ++j) {
        flat[k].root = i; flat[k].child = (long)j; ++k;
      }
    }

    if (cur >= nflat) cur = nflat ? nflat - 1 : 0;

    /* Rebuild the detail pane before laying out, so its height is known and
     * the list shrinks to fit rather than the pane running off the screen. */
    for (size_t d = 0; d < det.n; ++d) free(det.rows[d]);
    free(det.rows);
    det.rows = NULL; det.n = 0;

    int rowsz = term_rows();
    int cols = term_cols();
    int avail = rowsz - 4;

    if (detail_open && spec->detail && nflat) {
      size_t dri = flat[cur].root;
      long dci = flat[cur].child;
      spec->detail(ctx, roots[dri],
                   (dci >= 0 && node[dri].kids.keys)
                     ? node[dri].kids.keys[dci] : NULL,
                   cols - 2, &det);
      /* Never let the pane crowd the list below a few rows: seeing which row
       * the text describes is what makes it useful at all. */
      int cap = rowsz / 2;
      if (cap < 0) cap = 0;
      if ((int)det.n > cap && cap > 0) {
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
      fit(buf, cols - 6, cut, sizeof(cut));
      /* Indent matches the cursor marker plus the fold arrow. */
      frame_line(&f, "    %s%s%s", yame_ui_dim(), cut, yame_ui_reset());
    } else {
      frame_line(&f, "");
    }

    for (int r = 0; r < avail; ++r) {
      size_t fi = top + (size_t)r;
      if (fi >= nflat) { frame_line(&f, ""); continue; }

      const flatrow_t *fr = &flat[fi];
      int is_cur = (fi == cur);
      const char *arrow = is_cur ? (yame_ui_unicode() ? "❯" : ">") : " ";

      if (fr->child < 0) {
        /* A root: fold marker, then aligned columns. */
        treenode_t *nd = &node[fr->root];
        const char *fold = nd->expanded ? (yame_ui_unicode() ? "▾" : "-")
                                        : (yame_ui_unicode() ? "▸" : "+");
        char buf[1024]; size_t o = 0;
        int c = 0;
        for (const char *p = roots[fr->root]; c < ncol && o + 40 < sizeof(buf); ++c) {
          const char *tab = strchr(p, '\t');
          int len = (int)(tab ? (size_t)(tab - p) : strlen(p));
          if (len > w[c] - 2) len = w[c] - 2;
          o += (size_t)snprintf(buf + o, sizeof(buf) - o, "%-*.*s", w[c], len, p);
          if (!tab) break;
          p = tab + 1;
        }
        char cut[1024];
        fit(buf, cols - 6, cut, sizeof(cut));

        frame_line(&f, "%s%s%s %s%s%s %s%s%s%s",
                   yame_ui_cyan(), arrow, yame_ui_reset(),
                   yame_ui_cyan(), fold, yame_ui_reset(),
                   is_cur ? yame_ui_bold() : "",
                   style_color(root_styles, fr->root), cut, yame_ui_reset());
      } else {
        /* A child: preformatted, indented under its parent. */
        const treenode_t *nd = &node[fr->root];
        int last = ((size_t)fr->child + 1 == nd->kids.n);
        const char *stem = yame_ui_unicode() ? (last ? "└" : "├") : "|";
        char cut[1024];
        fit(nd->kids.rows[fr->child], cols - 14, cut, sizeof(cut));

        /* Neither a row already present nor a required one needs a checkbox:
         * in both cases there is nothing to ask for, so each shows what it is
         * instead of pretending to be a choice. */
        int have = nd->kids.styles && !spec->have_selectable &&
                   nd->kids.styles[fr->child] == YAME_ROW_HAVE;
        int required = nd->kids.styles &&
                       nd->kids.styles[fr->child] == YAME_ROW_REQUIRED;
        /* Every marker is four cells wide so the text after it lines up, and
         * the check sits where the x would be in "[x]". */
        const char *box = "";
        if (picking) {
          if (required)  box = yame_ui_unicode() ? " ●  " : " ** ";
          else if (have) box = yame_ui_unicode() ? " ✓  " : " ok ";
          else if (nd->kids.keys && nd->checked)
            box = nd->checked[fr->child] ? "[x] " : "[ ] ";
          else box = "    ";
        }

        frame_line(&f, "%s%s%s  %s%s%s %s%s%s%s%s",
                   yame_ui_cyan(), arrow, yame_ui_reset(),
                   yame_ui_dim(), stem, yame_ui_reset(),
                   is_cur ? yame_ui_bold() : "",
                   style_color(nd->kids.styles, (size_t)fr->child),
                   box, cut, yame_ui_reset());
      }
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

    if (picking) {
      size_t nsel = 0;
      for (size_t i = 0; i < n_roots; ++i)
        for (size_t j = 0; j < node[i].kids.n; ++j)
          if (node[i].checked && node[i].checked[j]) ++nsel;
      char acts[160];
      size_t ao = 0;
      if (spec->recommend)
        ao += (size_t)snprintf(acts + ao, sizeof(acts) - ao, "r recommended  ");
      if (spec->detail && spec->detail_key)
        ao += (size_t)snprintf(acts + ao, sizeof(acts) - ao, "%c %s  ",
                               spec->detail_key,
                               detail_open ? "close"
                                           : (spec->detail_verb
                                              ? spec->detail_verb : "info"));
      for (size_t a = 0; a < spec->n_actions && ao + 24 < sizeof(acts); ++a)
        ao += (size_t)snprintf(acts + ao, sizeof(acts) - ao, "%c %s  ",
                               spec->actions[a].key, spec->actions[a].verb);
      frame_line(&f, "%s  row %zu of %zu  %s  %zu selected  %s  %s open  "
                     "%s close  space select  %s%s%s q quit%s",
                 yame_ui_dim(), nflat ? cur + 1 : 0, nflat,
                 yame_ui_bullet(), nsel, yame_ui_bullet(),
                 yame_ui_unicode() ? "→" : "right",
                 yame_ui_unicode() ? "←" : "left",
                 acts, spec->hint ? spec->hint : "", spec->hint ? "  " : "",
                 yame_ui_reset());
    } else {
      frame_line(&f, "%s  row %zu of %zu   %s  %s open  %s close  q quit%s",
                 yame_ui_dim(), nflat ? cur + 1 : 0, nflat, motion,
                 yame_ui_unicode() ? "→" : "right",
                 yame_ui_unicode() ? "←" : "left",
                 yame_ui_reset());
    }

    frame_finish(&f);
    fflush(stderr);

    char ch = 0;
    keycode_t key = read_key(&ch);
    if (interrupted) break;

    /* Opening is lazy: children are asked for the first time a row unfolds. */
    size_t ri = nflat ? flat[cur].root : 0;
    long   ci = nflat ? flat[cur].child : -1;

    int want_open = (key == K_RIGHT || key == K_ENTER ||
                     (key == K_CHAR && ch == 'l'));
    int want_close = (key == K_LEFT || (key == K_CHAR && ch == 'h'));

    if (want_open && nflat && ci < 0) {
      treenode_t *nd = &node[ri];
      if (!nd->loaded) {
        if (expand) expand(ctx, roots[ri], &nd->kids);
        if (picking && nd->kids.n) nd->checked = calloc(nd->kids.n, 1);
        nd->loaded = 1;
      }
      nd->expanded = nd->kids.n ? 1 : 0;
    } else if (want_close && nflat) {
      if (ci >= 0) {
        /* Closing from inside jumps back to the parent, the way a file tree
         * behaves -- otherwise the cursor lands on whatever replaced it. */
        node[ri].expanded = 0;
        size_t back = 0;
        for (size_t i = 0; i < n_roots; ++i) {
          if (i == ri) break;
          back += 1 + (node[i].expanded ? node[i].kids.n : 0);
        }
        cur = back;
      } else {
        node[ri].expanded = 0;
      }
    }
    else if (picking && key == K_SPACE && nflat) {
      treenode_t *nd = &node[ri];
      if (ci >= 0) {
        int have = nd->kids.styles &&
                   ((!spec->have_selectable &&
                     nd->kids.styles[ci] == YAME_ROW_HAVE) ||
                    nd->kids.styles[ci] == YAME_ROW_REQUIRED);
        if (nd->checked && nd->kids.keys && !have)
          nd->checked[ci] = !nd->checked[ci];
        if (cur + 1 < nflat) ++cur;
      } else {
        /* Space on a parent takes everything under it that is not already
         * present, opening it first so the effect is visible. */
        if (!nd->loaded) {
          if (expand) expand(ctx, roots[ri], &nd->kids);
          if (nd->kids.n) nd->checked = calloc(nd->kids.n, 1);
          nd->loaded = 1;
        }
        if (nd->kids.n) nd->expanded = 1;
        int any = 0;
        for (size_t j = 0; j < nd->kids.n; ++j)
          if (nd->checked && nd->checked[j]) any = 1;
        for (size_t j = 0; j < nd->kids.n; ++j) {
          int have = nd->kids.styles &&
                     ((!spec->have_selectable &&
                       nd->kids.styles[j] == YAME_ROW_HAVE) ||
                      nd->kids.styles[j] == YAME_ROW_REQUIRED);
          if (nd->checked && nd->kids.keys && !have) nd->checked[j] = !any;
        }
      }
    }
    else if (picking && key == K_CHAR && find_action(spec, ch) >= 0) {
      const yame_ui_action_t *act = &spec->actions[find_action(spec, ch)];
      size_t nsel = 0;
      for (size_t i = 0; i < n_roots; ++i)
        for (size_t j = 0; j < node[i].kids.n; ++j)
          if (node[i].checked && node[i].checked[j]) ++nsel;
      if (nsel) {
        for (size_t i = 0; i < n_roots; ++i)
          for (size_t j = 0; j < node[i].kids.n; ++j)
            if (node[i].checked && node[i].checked[j] && node[i].kids.keys &&
                act->accept)
              act->accept(ctx, roots[i], node[i].kids.keys[j]);

        if (!act->commit) {    /* caller wants the selection, not the work */
          accepted = find_action(spec, ch) + 1;
          break;
        }

        /* The panel simply paints over the bottom rows, pane included. The
         * catalogue stays visible above it, which is the part that matters
         * while a download runs; the description underneath is context you
         * have already read. */
        act->commit(ctx);
        yame_ui_panel_close();

        /* What was fetched is now present, so every loaded expansion is
         * stale. Reload the open ones in place, keeping folds and cursor. */
        for (size_t i = 0; i < n_roots; ++i) {
          if (!node[i].loaded) continue;
          for (size_t j = 0; j < node[i].kids.n; ++j) {
            free(node[i].kids.rows[j]);
            if (node[i].kids.keys) free(node[i].kids.keys[j]);
          }
          free(node[i].kids.rows);
          free(node[i].kids.keys);
          free(node[i].kids.styles);
          free(node[i].checked);
          memset(&node[i].kids, 0, sizeof(node[i].kids));
          node[i].checked = NULL;
          node[i].loaded = 0;

          if (node[i].expanded && expand) {
            expand(ctx, roots[i], &node[i].kids);
            if (node[i].kids.n) node[i].checked = calloc(node[i].kids.n, 1);
            node[i].loaded = 1;
            if (!node[i].kids.n) node[i].expanded = 0;
          }
        }
        continue;
      }
    }
    else if (key == K_DOWN) { if (cur + 1 < nflat) ++cur; }
    else if (key == K_UP)   { if (cur) --cur; }
    else if (key == K_PGDN) { cur += (size_t)avail; if (cur >= nflat) cur = nflat ? nflat-1 : 0; }
    else if (key == K_PGUP) { cur = (cur > (size_t)avail) ? cur - (size_t)avail : 0; }
    else if (key == K_HOME) cur = 0;
    else if (key == K_END)  cur = nflat ? nflat - 1 : 0;
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
      else if (ch == 'j') { if (cur + 1 < nflat) ++cur; }
      else if (ch == 'k') { if (cur) --cur; }
      else if (picking && spec->recommend && ch == 'r' && nflat) {
        /* Scoped to the collection under the cursor. Recommending across
         * every collection at once would check sets for platforms the user is
         * not working on, and open all of them to do it. */
        size_t i = ri;
        if (!node[i].loaded) {
          if (expand) expand(ctx, roots[i], &node[i].kids);
          if (node[i].kids.n) node[i].checked = calloc(node[i].kids.n, 1);
          node[i].loaded = 1;
        }
        if (node[i].kids.n) node[i].expanded = 1;
        for (size_t j = 0; j < node[i].kids.n; ++j) {
          int req = node[i].kids.styles &&
                    ((!spec->have_selectable &&
                      node[i].kids.styles[j] == YAME_ROW_HAVE) ||
                     node[i].kids.styles[j] == YAME_ROW_REQUIRED);
          if (!req && node[i].checked && node[i].kids.keys &&
              node[i].kids.keys[j])
            /* Normalize to 0/1: checked[] is a flag, and a predicate that
             * returned, say, 256 to mean "yes" would truncate to 0 here. */
            node[i].checked[j] = spec->recommend(ctx, roots[i],
                                                 node[i].kids.keys[j]) ? 1 : 0;
        }
      }
      else if (on_key && on_key(ctx, ch, nflat ? roots[ri] : NULL,
                                (nflat && ci >= 0 && node[ri].kids.keys)
                                  ? node[ri].kids.keys[ci] : NULL)) {
        /* The caller changed what the roots describe -- a different store,
         * say -- so every expansion is about something else now. */
        yame_ui_panel_close();
        for (size_t i = 0; i < n_roots; ++i) {
          if (!node[i].loaded) continue;
          for (size_t j = 0; j < node[i].kids.n; ++j) {
            free(node[i].kids.rows[j]);
            if (node[i].kids.keys) free(node[i].kids.keys[j]);
          }
          free(node[i].kids.rows);
          free(node[i].kids.keys);
          free(node[i].kids.styles);
          free(node[i].checked);
          memset(&node[i].kids, 0, sizeof(node[i].kids));
          node[i].checked = NULL;
          node[i].loaded = 0;
          node[i].expanded = 0;
        }
        cur = top = 0;
      }
    }
  }

  for (size_t i = 0; i < n_roots; ++i) {
    for (size_t j = 0; j < node[i].kids.n; ++j) {
      free(node[i].kids.rows[j]);
      if (node[i].kids.keys) free(node[i].kids.keys[j]);
    }
    free(node[i].kids.rows);
    free(node[i].kids.keys);
    free(node[i].kids.styles);
    free(node[i].checked);
  }
  free(node);
  free(flat);
  raw_leave();
  for (size_t d = 0; d < det.n; ++d) free(det.rows[d]);
  free(det.rows);

  return accepted;
}
