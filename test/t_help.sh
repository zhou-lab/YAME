#!/bin/bash
## Every subcommand's -h runs, says something, and does not crash. Cheap, and
## it is the only thing that exercises the usage text at all -- a usage
## function that segfaults or reads past its option table fails here rather
## than in front of a user asking for help.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}

## The banner lists the subcommands; take them from it rather than a hand list,
## so a new subcommand is covered the day it is added.
## The bare banner IS a usage, so it exits non-zero by design. Capture it once
## with `|| true` -- piping it directly would fail the whole pipeline under
## `set -o pipefail` no matter what the downstream filter found.
banner=$("$YAME" </dev/null 2>&1) || true
cmds=$(printf '%s\n' "$banner" | awk '/^  [a-z][a-z0-9]+ +[A-Z]/ { print $1 }' | sort -u)
[ -n "$cmds" ] || { echo "could not read the subcommand list from the banner"; exit 1; }
n=0
for c in $cmds; do
  out=$("$YAME" "$c" -h </dev/null 2>&1) || true     # -h conventionally exits 1
  case "$out" in
    *Usage*) ;;
    *) echo "$c -h did not print a usage block"; printf '%s\n' "$out" | head -3; exit 1 ;;
  esac
  n=$((n + 1))
done
[ "$n" -ge 15 ] || { echo "only $n subcommands found; the banner parse is wrong"; exit 1; }

## The bare banner carries the compiled-in version, and an unknown subcommand
## is an error rather than a silent no-op.
printf '%s\n' "$banner" | grep 'v[0-9]\+\.[0-9]\+' >/dev/null ||
  { echo "banner has no version"; exit 1; }
if "$YAME" definitely-not-a-subcommand </dev/null >/dev/null 2>&1; then
  echo "an unknown subcommand exited 0"; exit 1
fi

## ---- the version, and a terminal-safe default -----------------------------
## `yame --version` was "unrecognized command": the version only appeared in
## the bare banner, which goes to stderr and exits 1, so no packaging check or
## script could read it. Reported by a reader, 2026-09-17.
for f in --version -V -v version; do
  out=$("$YAME" $f </dev/null 2>/dev/null) || { echo "yame $f exited nonzero"; exit 1; }
  printf '%s\n' "$out" | grep -E '^yame v[0-9]+\.[0-9]+' >/dev/null ||
    { echo "yame $f printed [$out], want a version on stdout"; exit 1; }
done
## it must be on STDOUT, not stderr, or `$(yame --version)` is empty
[ -n "$("$YAME" --version </dev/null 2>/dev/null)" ] ||
  { echo "yame --version wrote nothing to stdout"; exit 1; }

## --help reaches the same banner as a bare yame
h=$("$YAME" --help </dev/null 2>&1) || true
printf '%s\n' "$h" | grep 'Yet Another Methylation Encoder' >/dev/null ||
  { echo "yame --help is not the banner"; exit 1; }

## an unknown command still fails, and says where to look
if "$YAME" nosuchcmd </dev/null >/dev/null 2>&1; then
  echo "an unknown command exited 0"; exit 1
fi
