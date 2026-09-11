#!/bin/bash
## Run every test/t_*.sh against the built binary and tally the results.
## Each test is self-contained: it makes its own fixtures, exits non-zero on
## failure, and prints what went wrong. Nothing here needs the network or the
## data store, so this runs inside the conda build on every platform.
##
##   make test               # or: YAME=/path/to/yame bash test/run.sh
set -uo pipefail

here=$(cd "$(dirname "$0")" && pwd)
export YAME=${YAME:-$here/../yame}
if [ ! -x "$YAME" ]; then
  echo "no binary at $YAME (run make first)" >&2
  exit 2
fi

pass=0; fail=0
for t in "$here"/t_*.sh; do
  name=$(basename "$t" .sh)
  if out=$(bash "$t" 2>&1); then
    pass=$((pass + 1)); echo "ok    $name"
  else
    fail=$((fail + 1)); echo "FAIL  $name"
    printf '%s\n' "$out" | sed 's/^/      /'
  fi
done
echo "$pass passed, $fail failed"
[ "$fail" -eq 0 ]
