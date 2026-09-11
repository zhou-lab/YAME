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
command -v python3 >/dev/null || { echo "skip: no python3 for the pty" >&2; exit 0; }

d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
export YAME_DATA_HOME="$d/store"; mkdir -p "$YAME_DATA_HOME"

python3 - "$YAME" <<'PY'
import os, pty, sys, time, select, signal

YAME = sys.argv[1]
fails = 0

def drive(args, keys, label, want_exit=0, settle=0.15, timeout=8.0):
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
    def pump(t):
        nonlocal out
        end = time.time() + t
        while time.time() < end:
            r, _, _ = select.select([fd], [], [], 0.05)
            if not r: continue
            try: out += os.read(fd, 65536)
            except OSError: return False
        return True
    pump(0.6)                              # first frame
    alive = True
    for k in keys:
        if not alive: break
        try: os.write(fd, k)
        except OSError: alive = False; break
        alive = pump(settle)
    # wait for exit, but never forever
    t0 = time.time(); status = None
    while time.time() - t0 < timeout:
        p, st = os.waitpid(pid, os.WNOHANG)
        if p: status = st; break
        pump(0.05)
    if status is None:
        os.kill(pid, signal.SIGKILL); os.waitpid(pid, 0)
        print(f"  FAIL {label}: still running after the keys; killed"); fails += 1
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

def frame_has(out, text, label):
    global fails
    if text.encode() not in out:
        print(f"  FAIL {label}: expected {text!r} somewhere in what was drawn"); fails += 1

UP, DOWN, PGDN, PGUP, HOME, END = b"\x1b[A", b"\x1b[B", b"\x1b[6~", b"\x1b[5~", b"\x1b[H", b"\x1b[F"
ESC, ENTER, BS = b"\x1b", b"\r", b"\x7f"

# 1. open, render the header, quit with q
code, out = drive(["fetch"], [b"q"], "open and quit")
frame_has(out, "yame fetch", "header")
frame_has(out, "\x1b[?1049h", "alternate screen entered")
frame_has(out, "\x1b[?1049l", "alternate screen left on exit")

# 2. every navigation key it documents, then q
keys = [b"j", b"j", b"k", DOWN, DOWN, UP, PGDN, PGUP, END, HOME, b" ", b"a", b"n", b"q"]
drive(["fetch"], keys, "navigation keys")

# 3. Escape means BACK or CLEAR, never quit: at the top level it redraws and
#    waits, and q is the way out. Encoded here so a change in that convention
#    is a visible test change rather than a surprise.
code, out = drive(["fetch"], [ESC, b"q"], "escape at top, then q")
frame_has(out, "yame fetch", "redrawn after escape")

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
time.sleep(0.6)
try: os.read(fd, 65536)
except OSError: pass
os.write(fd, b"jjq"); time.sleep(0.3)
try: os.read(fd, 65536)
except OSError: pass
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

# 9. off a terminal the browser must refuse cleanly rather than hang or crash
import subprocess
p = subprocess.run([YAME, "fetch"], stdin=subprocess.DEVNULL, capture_output=True, timeout=10)
if p.returncode < 0:
    print(f"  FAIL no-tty: died on signal {-p.returncode}"); fails += 1
if p.returncode == 0 and not p.stderr and not p.stdout:
    print("  FAIL no-tty: exited 0 with nothing said"); fails += 1

if fails:
    print(f"{fails} browser assertion(s) failed"); sys.exit(1)
PY
