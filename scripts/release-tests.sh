#!/bin/bash
## The release tests: the suite run every way a release needs, unattended and
## all at once. Step 3 of the release SOP in one command.
##
## Four of the five lanes run the SAME suite (test/run.sh). They differ only in
## how the binary was built, or which shell runs them. The first lane adds the
## one check that is not a test: the coverage number in the badge against the
## number the suite actually measures.
##
## Three build lanes, plus two that reuse the in-tree binary:
##
##   tree     the in-tree build: suite, then the coverage badge check
##   ndebug   -O3 -DNDEBUG, which is what conda-forge compiles and what
##            deletes every assert(); three seeks once vanished that way
##   ubtrap   -fsanitize=undefined -fsanitize-undefined-trap-on-error, the
##            trap form, because the lab boxes have no libubsan and no sudo
##   bash32   the same suite under bash 3.2 and with it first in PATH, which
##            is what macOS ships and the one dialect CI cannot show us early
##   layer5   t_store_info.sh against the real shared store (needs network or
##            a populated YAME_DATA_HOME); skips cleanly without one
##
## On a 4-core box the lanes mostly trade CPU rather than add throughput -- the
## suite already runs its tests JOBS-wide -- so the win here is that all five
## run from one command instead of five, and that none is forgotten. JOBS is
## split across the lanes so they do not oversubscribe.
##
## Usage: bash scripts/release-tests.sh [-q]   (-q: only the summary)
set -uo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root"

quiet=0
[ "${1:-}" = "-q" ] && quiet=1

ncpu=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
## Three build lanes share the cores; at least 1 apiece.
lane_jobs=$(( ncpu / 3 )); [ "$lane_jobs" -lt 1 ] && lane_jobs=1
B32=${B32:-$HOME/tmp/yame/bash32/bash-3.2/bash}

t_start=$(date +%s)
work=$(mktemp -d "${TMPDIR:-/tmp}/release-tests.XXXXXX")
logs="$work/logs"; mkdir -p "$logs"
trap 'rm -rf "$work"' EXIT

say() { [ "$quiet" = 1 ] || printf '%s\n' "$*"; }

## Each lane records its own seconds, so one run reports both the wall clock
## and what the same lanes would have cost one after another.
stamp() { date +%s > "$logs/$1.t0"; }
elapsed() { echo $(( $(date +%s) - $(cat "$logs/$1.t0") )) > "$logs/$1.s"; }

## A lane: copy the tree, build it with the given CC, run the suite.
build_lane() {
  ## Separate statements: every word of a `local` is expanded before the
  ## builtin runs, so `d="$work/$name"` on the same line reads an unset name
  ## and trips `set -u`.
  local name=$1
  local cc=$2
  local d="$work/$name"
  mkdir -p "$d"
  rsync -a --exclude '*.o' --exclude yame --exclude libyame.a \
        --exclude 'htslib/libhts.a' --exclude '.git' "$root/" "$d/" || return 1
  ( cd "$d" \
    && make -B CC="$cc" -j"$lane_jobs" \
    && make lib CC="$cc" \
    && JOBS="$lane_jobs" YAME="$d/yame" bash test/run.sh )
}

## ---- start the lanes --------------------------------------------------------
stamp ndebug
( build_lane ndebug "cc -O3 -DNDEBUG"; rc=$?; elapsed ndebug; exit $rc ) > "$logs/ndebug" 2>&1 &
p_ndebug=$!
stamp ubtrap
( build_lane ubtrap \
    "cc -O1 -g -fsanitize=undefined -fsanitize-undefined-trap-on-error"; \
  rc=$?; elapsed ubtrap; exit $rc ) > "$logs/ubtrap" 2>&1 &
p_ubtrap=$!

## The in-tree lane goes here in the foreground: the two lanes after it reuse
## the binary it produces, and rebuilding it in a copy would only duplicate
## work the developer has already done.
stamp tree
(
  make -j"$lane_jobs" && make lib &&
  JOBS="$lane_jobs" YAME="$root/yame" bash test/run.sh &&
  ./scripts/coverage.sh --check
) > "$logs/tree" 2>&1
rc_tree=$?; elapsed tree

## bash 3.2 and layer 5 both want the in-tree binary, so they start now and
## run against the lanes still finishing above.
stamp bash32; stamp layer5
(
  if [ ! -x "$B32" ]; then
    echo "skip: no bash 3.2 at $B32 (see the release SOP for how to build one)"
    elapsed bash32; exit 0
  fi
  bad=0
  for f in test/*.sh scripts/*.sh; do
    "$B32" -n "$f" || { echo "PARSE FAIL $f"; bad=1; }
  done
  [ "$bad" = 0 ] || { elapsed bash32; exit 1; }
  PATH="$(dirname "$B32")_bin:$PATH" JOBS="$lane_jobs" \
    YAME="$root/yame" "$B32" test/run.sh
  rc=$?; elapsed bash32; exit $rc
) > "$logs/bash32" 2>&1 &
p_bash32=$!

( YAME_TEST_LAYER5=1 YAME="$root/yame" bash test/t_store_info.sh
  rc=$?; elapsed layer5; exit $rc ) > "$logs/layer5" 2>&1 &
p_layer5=$!

wait $p_ndebug; rc_ndebug=$?
wait $p_ubtrap; rc_ubtrap=$?
wait $p_bash32; rc_bash32=$?
wait $p_layer5; rc_layer5=$?

## ---- report -----------------------------------------------------------------
fails=0; serial=0
for g in tree ndebug ubtrap bash32 layer5; do
  eval "rc=\$rc_$g"
  tail=$(grep -E '^[0-9]+ passed|^skip:|^ok:' "$logs/$g" | tail -1)
  secs=$(cat "$logs/$g.s" 2>/dev/null || echo 0)
  serial=$((serial + secs))
  if [ "$rc" = 0 ]; then
    printf '  ok    %-8s %4s s  %s\n' "$g" "$secs" "$tail"
  else
    fails=$((fails + 1))
    printf '  FAIL  %-8s %4s s  %s\n' "$g" "$secs" "$tail"
    ## The suite prints "FAIL <name>" followed by that test's own output, and
    ## it can be anywhere in the log -- a plain tail showed the end of a green
    ## run and named nothing. Show those lines first, then the tail.
    if grep -q '^FAIL ' "$logs/$g"; then
      grep -A12 '^FAIL ' "$logs/$g" | sed 's/^/          /'
    else
      sed 's/^/          /' "$logs/$g" | tail -25
    fi
  fi
done
printf '%d of 5 lanes passed in %d s; one after another they would be %d s\n' \
  $((5 - fails)) $(( $(date +%s) - t_start )) "$serial"
[ "$fails" = 0 ]
