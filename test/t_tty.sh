#!/bin/bash
## What the tool does when its output IS a terminal, and when it is not.
##
## Both of these were reported by a reader on 2026-09-17, and neither could
## have been caught by a test that only ever redirects: every test here writes
## to a file or a pipe, which is exactly the case that behaved correctly.
## `script` gives us a pty, so the terminal case is testable.
##
##   1. a CX stream must never land on a terminal. Compressed binary garbles
##      the display. The page prints `yame pairwise -H 1 -c 5 -d 0.2 a.cg b.cg`
##      with no -o, so a reader copying it got exactly that.
##   2. hprint must not colour a pipe or a file, and must honour NO_COLOR.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}

command -v script >/dev/null || { echo "skip: no script(1) for a pty" >&2; exit 0; }
## script(1) exists in the conda test image but cannot get a pty there, so the
## terminal cases have nothing to test. Prove a pty works before relying on it.
script -qec true /dev/null >/dev/null 2>&1 ||
  { echo "skip: script(1) cannot allocate a pty here" >&2; exit 0; }
d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
cd "$d"

printf '5\t0\n0\t5\n3\t3\n4\t1\n' > mu.txt
"$YAME" pack -f 3 mu.txt > mu.cg
printf '0\n1\n0\n1\n' > b.txt
"$YAME" pack -f b b.txt > b.cg

## `script` runs a command with a pty on stdout. -q quiet, -e take its status.
onpty() { script -qec "$*" /dev/null 2>/dev/null; }

## ---- 1. binary to a terminal is refused, with a way out -------------------
for cmd in "binarize mu.cg" "mask mu.cg b.cg" "perturb b.cg" "dsample -N 2 mu.cg"; do
  out=$(onpty "$YAME $cmd" || true)
  printf '%s\n' "$out" | grep -i 'refusing to write compressed binary' >/dev/null ||
    { echo "yame $cmd wrote to a terminal without refusing"
      printf '%s\n' "$out" | head -3 | cat -v; exit 1; }
  printf '%s\n' "$out" | grep -- '-o' >/dev/null ||
    { echo "yame $cmd refused without naming -o"; exit 1; }
done

## ---- and the same commands still work when redirected --------------------
"$YAME" binarize mu.cg > piped.cg 2>/dev/null
[ -s piped.cg ] || { echo "binarize wrote nothing when redirected"; exit 1; }
"$YAME" info piped.cg >/dev/null 2>&1 || { echo "the redirected output is not readable"; exit 1; }
## and through a pipe, which is neither a file nor a terminal
n=$("$YAME" binarize mu.cg 2>/dev/null | "$YAME" info - 2>/dev/null | tail -n +2 | grep -c .)
[ "$n" -eq 1 ] || { echo "binarize through a pipe gave $n records, want 1"; exit 1; }

## ---- 2. hprint colours a terminal, and nothing else -----------------------
esc=$(printf '\033')
plain=$("$YAME" hprint b.cg 2>/dev/null | cat -v)
case $plain in
  *"^["*) echo "hprint coloured a PIPE: $plain"; exit 1 ;;
esac
"$YAME" hprint b.cg > file.txt 2>/dev/null
grep -q "$esc" file.txt && { echo "hprint coloured a REDIRECT"; exit 1; }

tty_out=$(onpty "$YAME hprint b.cg" | cat -v || true)
case $tty_out in
  *"^["*) ;;
  *) echo "hprint did not colour a terminal: $tty_out"; exit 1 ;;
esac

nc=$(NO_COLOR=1 script -qec "$YAME hprint b.cg" /dev/null 2>/dev/null | cat -v || true)
case $nc in
  *"^["*) echo "hprint ignored NO_COLOR on a terminal: $nc"; exit 1 ;;
esac

## -c still forces plain, terminal or not
fc=$(onpty "$YAME hprint -c b.cg" | cat -v || true)
case $fc in
  *"^["*) echo "hprint -c still coloured"; exit 1 ;;
esac

echo "ok: t_tty"
