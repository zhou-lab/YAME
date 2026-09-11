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
printf '%s\n' "$banner" | grep -q 'v[0-9]\+\.[0-9]\+' ||
  { echo "banner has no version"; exit 1; }
if "$YAME" definitely-not-a-subcommand </dev/null >/dev/null 2>&1; then
  echo "an unknown subcommand exited 0"; exit 1
fi
