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

## Every test makes its own temp directory and touches nothing shared, so
## they run concurrently -- JOBS at a time, default 4 -- and the wall time is
## the slowest test rather than the sum. Output is captured per test and
## printed in name order once all are done, so the report reads the same as
## a serial run. JOBS=1 gives the serial run back.
JOBS=${JOBS:-4}
logs=$(mktemp -d); trap 'rm -rf "$logs"' EXIT
running=0
for t in "$here"/t_*.sh; do
  name=$(basename "$t" .sh)
  ( if bash "$t" > "$logs/$name.out" 2>&1; then : > "$logs/$name.ok"; fi ) &
  running=$((running + 1))
  if [ "$running" -ge "$JOBS" ]; then wait -n 2>/dev/null || wait; running=$((running - 1)); fi
done
wait

pass=0; fail=0
for t in "$here"/t_*.sh; do
  name=$(basename "$t" .sh)
  if [ -e "$logs/$name.ok" ]; then
    pass=$((pass + 1)); echo "ok    $name"
  else
    fail=$((fail + 1)); echo "FAIL  $name"
    sed 's/^/      /' "$logs/$name.out"
  fi
done
echo "$pass passed, $fail failed"
[ "$fail" -eq 0 ]
