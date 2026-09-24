#!/bin/bash
## The interactive catalogue browser (`yame fetch` with no arguments, and
## `summary -b`), driven through a pseudo-terminal. It reads single keys in
## raw mode and refuses to start unless stdin and stderr are terminals, so a
## pipe cannot reach it; a pty from Python's standard library can.
##
## Nothing is downloaded: browsing walks the compiled-in registry and the
## local store, and every run ends with q or Escape before a fetch is
## confirmed. What is asserted is that it starts, renders, takes every key it
## documents without crashing, and exits 0 -- the crash-in-front-of-the-user
## class, which no other test can reach.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
## The browser is driven through a pty, so it is the one test whose result
## depends on TIMING rather than on output. Run under load it reads keys late
## and a swallowed one changes what the next key means. The release tests run
## it ONCE, alone, after the build lanes finish, and set this in the lanes.
if [ -n "${YAME_SKIP_UI:-}" ]; then
  echo "skip: t_ui runs alone (YAME_SKIP_UI set)"; exit 0
fi
command -v python3 >/dev/null || { echo "skip: no python3 for the pty" >&2; exit 0; }

d=$(mktemp -d); srv=; trap '[ -n "$srv" ] && kill $srv 2>/dev/null; rm -rf "$d"' EXIT
export YAME_DATA_HOME="$d/store"; mkdir -p "$YAME_DATA_HOME"

## A local mirror serving the one real set the repo carries, so a picker run
## can choose it and have it fetched, as t_fetch.sh does, with no network.
here=$(cd "$(dirname "$0")" && pwd)
tree="$d/mirror/zhou-lab/InfiniumAnnotation/v8.2/EPIC/KYCG"
mkdir -p "$tree" && cp "$here/fixtures/Blacklist.20220304.cm" "$tree/"
port=$(python3 -c 'import socket; s=socket.socket(); s.bind(("127.0.0.1",0)); print(s.getsockname()[1])')
( cd "$d/mirror" && exec python3 -m http.server --bind 127.0.0.1 "$port" ) >"$d/server.log" 2>&1 &
srv=$!
## Up to 30 s: under the release lanes a loaded box has taken over 5 s to
## start python, and a test that went on anyway failed its first fetch with
## "couldn't connect", which read as a fetch bug. Never up is its own error.
up=0
for i in $(seq 1 300); do
  python3 -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:$port/', timeout=1)" 2>/dev/null && { up=1; break; }
  sleep 0.1
done
[ "$up" = 1 ] || { echo "the local mirror on port $port never started"; cat "$d/server.log" 2>/dev/null; exit 1; }
export YAME_ASSETS_MIRROR="http://127.0.0.1:$port"
export FIXTURES="$here/fixtures"

python3 - "$YAME" <<'PY'
import os, pty, sys, time, select, signal

YAME = sys.argv[1]
fails = 0

## Quiet-detection thresholds. A finished frame goes quiet in a few ms on an
## idle box, but the release tests run six lanes on four cores, and there a
## frame can stall mid-draw for longer than the old 25 ms. The pump then
## returned early, the next key went out before the browser had finished
## reading the last one, and a swallowed `/` turned the following Enter into
## "open the info pane". Scaled by YAME_UI_SLOW for a loaded machine.
_SLOW = float(os.environ.get("YAME_UI_SLOW", "1"))
QUIET = 0.08 * _SLOW
SETTLE = 1.0 * _SLOW

def drive(args, keys, label, want_exit=0, settle=SETTLE, timeout=15.0):
    """Run yame under a pty, feed keys with a pause between them, collect
    everything it draws, and return (exit status, output)."""
    global fails
    pid, fd = pty.fork()
    if pid == 0:
        os.environ["TERM"] = "xterm"
        os.environ["COLUMNS"] = "100"; os.environ["LINES"] = "30"
        os.environ.pop("NO_COLOR", None)  # NO_COLOR or TERM=dumb means "plain
        os.execv(YAME, ["yame"] + args)    # terminal": a table, not the browser
    out = b""
    def pump(t, quiet=QUIET):
        """Read until the output has been quiet for `quiet` seconds, or `t`
        has elapsed -- whichever is first. A finished frame goes quiet in a
        few milliseconds, so this is the wait a keypress actually needs; a
        fixed sleep of the largest plausible time made this file take 24 s."""
        nonlocal out
        end = time.time() + t; last = time.time()
        while time.time() < end:
            r, _, _ = select.select([fd], [], [], 0.01)
            if r:
                try: out += os.read(fd, 65536); last = time.time()
                except OSError: return False
            elif time.time() - last >= quiet and out: return True
        return True
    pump(2.0 * _SLOW)                      # first frame, however long it takes
    alive = True
    for k in keys:
        if not alive: break
        try: os.write(fd, k)
        except OSError: alive = False; break
        alive = pump(settle, quiet=QUIET)
        ## A bare ESC is told apart from the start of an escape sequence by a
        ## 40 ms poll for a following byte (read_key in src/ui.c). A key that
        ## lands inside that window is taken for the sequence's tail and
        ## DISCARDED, so a `q` after ESC never arrives and the session runs to
        ## its timeout. 200 ms rather than 60: on a loaded CI runner the pump
        ## above can return while the frame is still coming, which shortens the
        ## real gap, and this failed exactly once that way on macOS.
        if k == b"\x1b": time.sleep(0.20)
    # wait for exit, but never forever
    t0 = time.time(); status = None
    while time.time() - t0 < timeout:
        p, st = os.waitpid(pid, os.WNOHANG)
        if p: status = st; break
        pump(0.05, quiet=0.01)
    if status is None:
        os.kill(pid, signal.SIGKILL); os.waitpid(pid, 0)
        print(f"  FAIL {label}: still running after the keys; killed")
        ## What it was sitting on. Without this a hang here is only reproducible
        ## by pushing again and reading the next run's log.
        tail = out[-400:].decode("utf8", "replace").replace("\x1b", "<ESC>")
        print(f"   keys sent: {b' '.join(keys)!r}")
        print(f"   last output: {tail}")
        fails += 1
        return None, out
    pump(0.1)
    if os.WIFSIGNALED(status):
        print(f"  FAIL {label}: died on signal {os.WTERMSIG(status)}"); fails += 1
        return None, out
    code = os.WEXITSTATUS(status)
    if code != want_exit:
        print(f"  FAIL {label}: exit {code}, want {want_exit}")
        print("   last output:", out[-300:].decode("utf8", "replace").replace("\x1b", "<ESC>"))
        fails += 1
    return code, out

def settle(fd, t=2.0 * _SLOW, quiet=QUIET):
    """Read from fd until it has been quiet for `quiet` s or `t` s pass; return
    the bytes. The wait every interactive check needs: it ends when the frame
    is finished, not on a timer."""
    out = b""; end = time.time() + t; last = time.time()
    while time.time() < end:
        r, _, _ = select.select([fd], [], [], 0.01)
        if r:
            try: out += os.read(fd, 65536); last = time.time()
            except OSError: break
        elif out and time.time() - last >= quiet: break
    return out

def frame_has(out, text, label):
    global fails
    if text.encode() not in out:
        print(f"  FAIL {label}: expected {text!r} somewhere in what was drawn"); fails += 1

UP, DOWN, PGDN, PGUP, HOME, END = b"\x1b[A", b"\x1b[B", b"\x1b[6~", b"\x1b[5~", b"\x1b[H", b"\x1b[F"
ESC, ENTER, BS = b"\x1b", b"\r", b"\x7f"

# 1. open, render the header, quit with q
code, out = drive(["fetch"], [b"q"], "open and quit")
frame_has(out, "YAME_DATA_HOME: ", "header")
frame_has(out, "\x1b[?1049h", "alternate screen entered")
frame_has(out, "\x1b[?1049l", "alternate screen left on exit")

# 2. every navigation key it documents, then q
keys = [b"j", b"j", b"k", DOWN, DOWN, UP, PGDN, PGUP, END, HOME, b" ", b"a", b"n", b"q"]
drive(["fetch"], keys, "navigation keys")

# 3. Escape means BACK or CLEAR, never quit: at the top level it redraws and
#    waits, and q is the way out. Encoded here so a change in that convention
#    is a visible test change rather than a surprise.
code, out = drive(["fetch"], [ESC, b"q"], "escape at top, then q")
frame_has(out, "YAME_DATA_HOME: ", "redrawn after escape")

# 4. filtering: type, backspace, apply with Enter, clear with Escape, quit
drive(["fetch"], [b"/", b"E", b"P", b"I", b"C", BS, ENTER, b"/", b"z", b"z", ESC, b"q"], "filter")

# 5. Enter on an entry descends (or selects); Escape backs out; q leaves.
#    Whatever the entry is, the browser must not crash on the transition.
drive(["fetch"], [ENTER, ESC, b"j", ENTER, ESC, b"q"], "enter and back")

# 6. a very narrow terminal must not crash the renderer
pid, fd = pty.fork()
if pid == 0:
    os.environ["TERM"] = "xterm"; os.environ["COLUMNS"] = "20"; os.environ["LINES"] = "6"
    os.execv(YAME, ["yame", "fetch"])
settle(fd)
os.write(fd, b"jjq"); settle(fd, t=1.0)
_, st = os.waitpid(pid, 0)
if os.WIFSIGNALED(st):
    print(f"  FAIL narrow terminal: died on signal {os.WTERMSIG(st)}"); fails += 1

# 7. a plain terminal gets a table instead of the browser, with no keys needed
pid, fd = pty.fork()
if pid == 0:
    os.environ["TERM"] = "dumb"; os.environ["COLUMNS"] = "100"; os.environ["LINES"] = "30"
    os.execv(YAME, ["yame", "fetch"])
out = b""; t0 = time.time()
while time.time() - t0 < 5:
    r, _, _ = select.select([fd], [], [], 0.1)
    if r:
        try: out += os.read(fd, 65536)
        except OSError: break
    p, st = os.waitpid(pid, os.WNOHANG)
    if p: break
else:
    os.kill(pid, signal.SIGKILL); os.waitpid(pid, 0)
    print("  FAIL TERM=dumb: did not exit on its own (it opened the browser?)"); fails += 1
    st = None
if st is not None and os.WIFEXITED(st) and os.WEXITSTATUS(st) != 0:
    print(f"  FAIL TERM=dumb: exit {os.WEXITSTATUS(st)}"); fails += 1
if b"\x1b[?1049h" in out:
    print("  FAIL TERM=dumb: opened the alternate screen anyway"); fails += 1
frame_has(out, "target\tsource", "TERM=dumb table header")

# 8. summary -b uses the same browser to choose masks
code, out = drive(["summary", "-b", "/dev/null"], [b"q"], "summary -b opens", want_exit=1)

# 9. deeper walks: open a subtree, mark several entries, search and clear,
#    page through, then leave. The tree is built lazily, so descending is what
#    reaches the flatten/reopen/refresh paths that a single frame never does.
deep = [ENTER, DOWN, ENTER, DOWN, b" ", DOWN, b" ", PGDN, PGUP, ESC, ESC, b"q"]
drive(["fetch"], deep, "descend, mark, back out")

# 10. select-all and select-none inside a subtree
drive(["fetch"], [ENTER, b"a", b"n", b"a", ESC, b"q"], "select all and none")

# 11. a search that matches, then one that matches nothing, then clear
code, out = drive(["fetch"], [b"/", b"h", b"g", b"3", b"8", ENTER, b"j", ENTER, ESC,
                              b"/", b"q", b"q", b"z", b"z", b"z", ESC, b"q"], "search hits and misses")

# 12. resize mid-session: the renderer reads COLUMNS/LINES, and SIGWINCH
#     arrives while it is waiting for a key
pid, fd = pty.fork()
if pid == 0:
    os.environ["TERM"] = "xterm"; os.environ["COLUMNS"] = "120"; os.environ["LINES"] = "40"
    os.execv(YAME, ["yame", "fetch"])
settle(fd)
import fcntl, termios, struct
try:
    fcntl.ioctl(fd, termios.TIOCSWINSZ, struct.pack("HHHH", 12, 40, 0, 0))
    os.kill(pid, signal.SIGWINCH)
except Exception:
    pass
settle(fd, t=1.0)                      # the redraw the signal triggers, if any
os.write(fd, b"j"); settle(fd, t=1.0)
os.write(fd, b"q"); settle(fd, t=1.0)
t0 = time.time(); st = None
while time.time() - t0 < 5:
    p, s_ = os.waitpid(pid, os.WNOHANG)
    if p: st = s_; break
    settle(fd, t=0.05, quiet=0.01)
if st is None:
    os.kill(pid, signal.SIGKILL); os.waitpid(pid, 0)
    print("  FAIL resize: did not exit after q"); fails += 1
elif os.WIFSIGNALED(st):
    print(f"  FAIL resize: died on signal {os.WTERMSIG(st)}"); fails += 1

# 13. the browser also drives mask selection for summary, over a real file
import subprocess, tempfile
tmpd = os.environ["YAME_DATA_HOME"]
cg = os.path.join(tmpd, "q.cg")
subprocess.run(f"awk 'BEGIN{{for(i=0;i<8;i++) print (i%2)}}' | {YAME} pack -f b - > {cg}",
               shell=True, check=True)
drive(["summary", "-b", cg], [b"j", b"q"], "summary -b over a real query", want_exit=1)
# 13b. -m beside -b names what arrives checked. A query in a row space the
#      build knows (HM27's 27,722 rows) opens the picker on that unit, and
#      CGI -- a set name, any case -- is ticked inside its knowledgebase.
cg27 = os.path.join(tmpd, "q27.cg")
subprocess.run(f"awk 'BEGIN{{for(i=0;i<27722;i++) print (i%2)}}' | {YAME} pack -f b - > {cg27}",
               shell=True, check=True)
code, out = drive(["summary", "-b", "-m", "cgi", cg27], [b"q"], "summary -b -m preselects", want_exit=1)
frame_has(out, "[x] CGI.", "the -m set arrives checked")
## (the platform mask is offered too, below the 15 KYCG sets -- past the
## bottom of this 30-line window, so not asserted as drawn)
if b"[x] HM27.hg38.mask" in out:
    print("  FAIL summary -b checked a mask -m did not name"); fails += 1
if b"HM27.ordering" in out:
    print("  FAIL summary -b offered the ordering, which is not a mask"); fails += 1
# 13c. choose with u: the picker fetches what is missing (from the mirror),
#      then summary runs once per chosen mask under a single header. An
#      EPIC-sized query opens EPIC; -m blacklist arrives checked.
cgE = os.path.join(tmpd, "qE.cg")
subprocess.run(f"awk 'BEGIN{{for(i=0;i<866553;i++) print (i%3?1:0)}}' | {YAME} pack -f b - > {cgE}",
               shell=True, check=True)
code, out = drive(["summary", "-b", "-m", "blacklist", cgE], [b"u"], "summary -b: u fetches and summarizes",
                  settle=SETTLE * 3)
frame_has(out, "QFile", "the summary table was printed")
if out.count(b"QFile") != 1:
    print(f"  FAIL summary -b printed the header {out.count(b'QFile')} times, want once"); fails += 1
frame_has(out, "Blacklist.20220304.cm", "the chosen mask was summarized")
if not os.path.exists(os.path.join(tmpd, "EPIC", "KYCG", "Blacklist.20220304.cm")):
    print("  FAIL summary -b: the chosen set was not fetched into the store"); fails += 1

# 14. the help screen, reached with h and listing the key groups it documents
code, out = drive(["fetch"], [b"h", ESC, b"q"], "help screen")
for section in ("MOVING", "CHOOSING"):
    frame_has(out, section, f"help screen lists {section}")
drive(["fetch"], [b"h", b"h", ESC, b"q"], "help twice")

# 15. the info pane: enter opens it, i closes it (the status bar says "i close")
drive(["fetch"], [ENTER, b"i", ENTER, b"i", b"q"], "info pane open and close")

# 16. the remaining keys the help screen documents: l/left to open and close,
#     x to select
RIGHT, LEFT = b"\x1b[C", b"\x1b[D"
drive(["fetch"], [b"l", DOWN, RIGHT, LEFT, LEFT, b"x", b"a", b"a", b"q"],
      "open, close, select")

# 17. f with nothing selected must not start a fetch
code, out = drive(["fetch"], [b"f", b"q"], "fetch with an empty selection")

# 18. off a terminal the browser must refuse cleanly rather than hang or crash
import subprocess
p = subprocess.run([YAME, "fetch"], stdin=subprocess.DEVNULL, capture_output=True, timeout=10)
if p.returncode < 0:
    print(f"  FAIL no-tty: died on signal {-p.returncode}"); fails += 1
if p.returncode == 0 and not p.stderr and not p.stdout:
    print("  FAIL no-tty: exited 0 with nothing said"); fails += 1

# 19. stale files: the browser ASKS before it opens -- which files, under
#     their directory, and whether to replace them now -- then, declined,
#     counts them on the unit and knowledgebase rows, marks the file row
#     and explains in the info pane. Staged as present files whose manifest
#     line records another digest. (`y` would fetch; not driven here.)
for sub, name in (("hg38", "cpg_nocontig.cr"), ("EPIC/KYCG", "Blacklist.20220304.cm")):
    sd = os.path.join(tmpd, sub); os.makedirs(sd, exist_ok=True)
    open(os.path.join(sd, name), "wb").close()
    open(os.path.join(sd, "SHA256SUMS"), "w").write("f" * 64 + "  " + name + "\n")
code, out = drive(["fetch"], [b"n\r", b"/", b"cpg_nocontig", ENTER, ENTER, b"q"], "stale files: dialog, then browser")
frame_has(out, "2 directories hold files from an earlier release", "dialog names the count")
frame_has(out, "      cpg_nocontig.cr", "dialog lists the file under its directory")
frame_has(out, "Replace them now?", "dialog asks")
frame_has(out, "1 stale: -f", "unit row counts its stale file")
frame_has(out, "stale ", "stale file row is marked beside its size")
frame_has(out, "earlier release", "info pane explains the stale file")   # wrapped, so two
frame_has(out, "zhou-lab/genomes", "info pane names the upstream")           # pieces, not one
if out.find(b"Replace them now?") > out.find(b"\x1b[2J"):
    print("  FAIL dialog: asked after the browser's first frame, not before"); fails += 1
code, out = drive(["fetch"], [b"n\r", b"/", b"EPIC", ENTER, b"l", b"q"], "stale count on a knowledgebase row")
frame_has(out, "KYCG        sets", "the EPIC unit opened to its knowledgebase row")
frame_has(out, "1 stale: -f", "knowledgebase row counts its stale file")
# -q: no dialog, straight to the browser
code, out = drive(["fetch", "-q"], [b"q"], "-q skips the dialog")
if b"Replace them now?" in out:
    print("  FAIL -q: the dialog was still asked"); fails += 1

# 21. naming a directory on a terminal opens the browser on it, with exactly
#     what the plan would move already checked: the browser is the
#     confirmation. q leaves without fetching.
code, out = drive(["fetch", "-g", "CGI", "EPIC/KYCG"], [b"q"], "a named directory opens checked")
frame_has(out, "[x] CGI.", "the planned file arrives checked")
frame_has(out, "[ ] ChromHMM.", "a file the -g filter left out stays unchecked")
if b"Proceed?" in out:
    print("  FAIL a named directory still asked Proceed? on a terminal"); fails += 1
# names in two units cannot open as one tree: the prompt stays, n declines
code, out = drive(["fetch", "EPIC", "MSA"], [b"n\r"], "two units keep the prompt", want_exit=1)
frame_has(out, "Proceed?", "two units are confirmed with the prompt")

# 26. each platform and genome row shows its row-space size, so a knowledge-
#     base can be matched to a query before anything is fetched
code, out = drive(["fetch", "-q"], [b"q"], "ROWS column")
frame_has(out, "ROWS", "the header names the ROWS column")
frame_has(out, "866,553", "EPIC's row shows its 866,553 probes")
frame_has(out, "29,401,795", "hg38's row shows its CpG count")

# 22. y to the stale dialog replaces those files -- here from the local
#     mirror -- and then the browser opens over a current store
import shutil
shutil.rmtree(os.path.join(tmpd, "hg38"), ignore_errors=True)   # 19's stale file: not mirrored
sd = os.path.join(tmpd, "EPIC", "KYCG"); os.makedirs(sd, exist_ok=True)
open(os.path.join(sd, "Blacklist.20220304.cm"), "wb").write(b"stale bytes")
open(os.path.join(sd, "SHA256SUMS"), "w").write("0" * 64 + "  Blacklist.20220304.cm\n")
code, out = drive(["fetch"], [b"y\r", b"q"], "stale dialog: y replaces", settle=SETTLE * 3)
frame_has(out, "Replace them now?", "the dialog asked")
fixture = open(os.path.join(os.environ["FIXTURES"], "Blacklist.20220304.cm"), "rb").read()
if open(os.path.join(sd, "Blacklist.20220304.cm"), "rb").read() != fixture:
    print("  FAIL y to the stale dialog did not replace the file"); fails += 1

# 23. a one-file fetch on a terminal draws a progress line
os.remove(os.path.join(sd, "Blacklist.20220304.cm"))
code, out = drive(["fetch", "EPIC/KYCG/Blacklist.20220304.cm"], [], "CLI fetch progress", settle=SETTLE * 3)
frame_has(out, "[1/1]", "the progress line counts the file")
frame_has(out, "100%", "the progress line reaches 100%")

# 25. a store set only by a sibling tool's variable is named as inherited
env_backup = dict(os.environ)
os.environ["METHSCOPE_DATA_HOME"] = os.environ.pop("YAME_DATA_HOME")
try:
    code, out = drive(["fetch", "-q"], [b"q"], "an inherited store")
    frame_has(out, "METHSCOPE_DATA_HOME (inherited)", "the title says whose variable set the store")
finally:
    os.environ.clear(); os.environ.update(env_backup)

# 20. d: point the browser at another store without leaving it. The second
#     store holds one stale file, so the switch shows in the rows as well as
#     in the title. Backspaces clear the current path from the prompt.
other = os.path.join(os.path.dirname(tmpd), "other")
od = os.path.join(other, "hg38"); os.makedirs(od, exist_ok=True)
open(os.path.join(od, "cpg_nocontig.cr"), "wb").close()
open(os.path.join(od, "SHA256SUMS"), "w").write("f" * 64 + "  cpg_nocontig.cr\n")
clear = b"\x7f" * (len(tmpd) + 8)
code, out = drive(["fetch", "-q"], [b"d", clear, other.encode(), b"\r", b"?", b" ", b"q"], "d switches the store")
frame_has(out, "store (d): " + other, "title names the new store")
frame_has(out, "1 stale: -f", "rows describe the new store")
frame_has(out, "change the store", "the help screen lists d")
# a path that is a file is refused, and the store stays as it was
notdir = os.path.join(os.path.dirname(tmpd), "afile"); open(notdir, "w").close()
code, out = drive(["fetch", "-q"], [b"d", clear, notdir.encode(), b"\r", ENTER, b"q"], "d refuses a file")
frame_has(out, "is not a directory", "d says why it refused")
if ("store (d): " + notdir).encode() in out:
    print("  FAIL d: switched to a path that is a file"); fails += 1

if fails:
    print(f"{fails} browser assertion(s) failed"); sys.exit(1)
PY

## ---- the parts of yame_ui.h only a downstream caller reaches -------------
## yame never asks for text inside a widget and never takes a signal in one,
## so yame_ui_panel_ask and the raw-mode signal handler had no coverage --
## and kycg asks for its store with the first and relies on the second to
## hand a Ctrl-C'd terminal back. test/probe_ui.c is the smallest caller of
## both, built against libyame.a as t_probe.sh builds its probe.
root=$(cd "$(dirname "$0")/.." && pwd)
if [ ! -f "$root/libyame.a" ] || [ ! -x "$root/yame-config" ]; then
  echo "skip: probe_ui needs libyame.a (run 'make lib')" >&2; exit 0
fi
## ${CC:-cc}: the coverage run instruments libyame.a; see t_probe.sh
${CC:-cc} -O1 -g -std=gnu99 $("$root/yame-config" --cflags) -o "$d/probe_ui" \
  "$root/test/probe_ui.c" $("$root/yame-config" --libs) 2>"$d/cc.err" ||
  { echo "probe_ui did not build"; cat "$d/cc.err"; exit 1; }

python3 - "$d/probe_ui" "$d" <<'PY'
import os, pty, sys, time, select, signal, termios
PROBE, D = sys.argv[1], sys.argv[2]
_SLOW = float(os.environ.get("YAME_UI_SLOW", "1"))
fails = 0

def run(keys, arg=None, cols="100", interrupt=False):
    """Run the probe on a pty, feed keys, and return (status, output, log,
    the pty's local modes after the probe is gone)."""
    log = os.path.join(D, "ask.log")
    if os.path.exists(log): os.remove(log)
    pid, fd = pty.fork()
    if pid == 0:
        os.environ.update(TERM="xterm", COLUMNS=cols, LINES="20", PROBE_UI_LOG=log)
        os.environ.pop("NO_COLOR", None)
        os.execv(PROBE, [PROBE] + ([arg] if arg else []))
    out = b""
    def pump(t):
        nonlocal out
        end = time.time() + t * _SLOW
        while time.time() < end:
            r, _, _ = select.select([fd], [], [], 0.05)
            if r:
                try: out += os.read(fd, 65536)
                except OSError: return
    pump(1.0)                                   # the tree's first frame
    for k in keys:
        os.write(fd, k); pump(0.3)
    if interrupt:
        os.kill(pid, signal.SIGINT); pump(0.5)
    lflag = termios.tcgetattr(fd)[3]
    _, st = os.waitpid(pid, 0)
    return st, out, (open(log).read() if os.path.exists(log) else ""), lflag

def check(cond, what):
    global fails
    if not cond: print("  FAIL " + what); fails += 1

## edit: three backspaces take /old to /, then type new; Enter accepts.
## Then a second ask, typed into and cancelled with Escape: rc 0.
st, out, log, _ = run([b"s", b"\x7f\x7f\x7f", b"new", b"\r", b"s", b"xx", b"\x1b", b"q"])
check(os.waitstatus_to_exitcode(st) == 0, "probe_ui did not exit 0")
check("RESULT rc=1 buf=[/new]" in log, "panel_ask: backspace and typing did not give /new: " + log)
check("RESULT rc=0 buf=" in log, "panel_ask: Escape did not return 0: " + log)
check(b"store:" in out, "panel_ask: the prompt was never drawn")

## a value longer than the line shows its tail, where the cursor is
long = "a" * 60 + "TAILEND"
st, out, log, _ = run([b"s", b"\r", b"q"], arg=long, cols="40")
check("buf=[" + long + "]" in log, "panel_ask: a long value did not come back whole")
check(b"TAILEND" in out and (b"a" * 60) not in out, "panel_ask: a long value was not shown by its tail")

## Ctrl-C while the widget holds the terminal: the process still dies of the
## signal, but only after leaving the alternate screen and restoring cooked
## mode -- or the user's shell is left raw and blank
st, out, log, lflag = run([], interrupt=True)
check(os.WIFSIGNALED(st) and os.WTERMSIG(st) == signal.SIGINT, "SIGINT did not end the probe by SIGINT")
check(b"\x1b[?1049l" in out, "SIGINT: the alternate screen was not left")
check(bool(lflag & termios.ICANON) and bool(lflag & termios.ECHO), "SIGINT: the terminal was left in raw mode")

## the picker as a downstream tool opens it (yame_browse_pick_opt): only the
## units it offers, its own title and verb, a unit open with a set checked
pid, fd = pty.fork()
if pid == 0:
    os.environ.update(TERM="xterm", COLUMNS="110", LINES="40")
    os.environ.pop("NO_COLOR", None)
    os.execv(PROBE, [PROBE, "pick"])
out = b""
end = time.time() + 1.5 * _SLOW
while time.time() < end:
    r, _, _ = select.select([fd], [], [], 0.05)
    if r:
        try: out += os.read(fd, 65536)
        except OSError: break
os.write(fd, b"q")
end = time.time() + 1.0 * _SLOW
while time.time() < end:
    r, _, _ = select.select([fd], [], [], 0.05)
    if r:
        try: out += os.read(fd, 65536)
        except OSError: break
_, st = os.waitpid(pid, 0)
check(os.waitstatus_to_exitcode(st) == 0, "pick: probe_ui did not exit 0")
check(b"probe pick" in out, "pick: the caller's title was not shown")
check(b"t test" in out, "pick: the caller's verb was not offered")
check(b"[x] CGI." in out, "pick: the preselected set did not arrive checked")
check(b"MSA" in out, "pick: an offered unit is missing")
## offer = */KYCG/*.cm: the row list and the platform mask are not choices
check(b"HM27.ordering" not in out, "pick: the ordering is offered, though offer allows only KYCG sets")
check(b"HM27.hg38.mask" not in out, "pick: the platform mask is offered outside KYCG")
## a genome unit's row reads "genome"; HM27's own files carry hg38 in their
## names, so the unit name itself is no test
check(b"EPICv2" not in out and b"genome" not in out, "pick: a unit that was not offered is listed")

if fails:
    print(f"{fails} probe_ui assertion(s) failed"); sys.exit(1)
PY
