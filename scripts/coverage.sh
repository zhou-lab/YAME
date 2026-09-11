#!/bin/sh
# scripts/coverage.sh — line coverage of test/run.sh over YAME's own sources.
#
#   scripts/coverage.sh            measure, print the table, rewrite docs/coverage.json
#   scripts/coverage.sh --print    measure and print only (leave the badge alone)
#   scripts/coverage.sh --check    measure and fail if docs/coverage.json is stale
#                                  by more than $TOLERANCE points (default 2.0)
#
# The number is gcov's own summary, counting EXECUTABLE lines, and it is LINE
# coverage rather than branch coverage (branch coverage is lower and is the
# more telling number for the format codecs -- measure it before quoting it).
#
# bgzf.c and kstring.c are excluded: they are vendored (htslib, klib), we do
# not maintain them, and their ~2000 lines would dominate a number meant to
# say how well OUR code is tested.
#
# Everything happens in a scratch copy. An instrumented binary must never
# become the one in the repo: it is built -O0 and writes .gcda files beside
# itself on every run.
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$here"
TOLERANCE=${TOLERANCE:-2.0}
mode=${1:---write}

command -v gcov >/dev/null || { echo "coverage: gcov not found (it ships with gcc)" >&2; exit 1; }

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT HUP INT TERM
cp -r src htslib test Makefile "$work/"
[ -d tools ] && cp -r tools "$work/"

# -O0 so line numbers map 1:1; optimisation merges and elides lines and the
# annotation stops meaning anything. The flags ride on CC because the Makefile
# ASSIGNS CFLAGS rather than appending to it -- a command-line CFLAGS would
# drop the -I flags the build needs. conda-recipe/build.sh does the same.
( cd "$work" && make -B CC="${CC:-cc} -O0 -g --coverage" LDFLAGS="--coverage" ) \
    >"$work/build.log" 2>&1
[ -x "$work/yame" ] || { tail -20 "$work/build.log" >&2; echo "coverage: build failed" >&2; exit 1; }

# Counters accumulate across processes, so clear them and let ONLY the suite
# run: a stray `yame -h` adds main.c's usage path and moves the total.
rm -f "$work"/src/*.gcda
YAME="$work/yame" sh "$work/test/run.sh" >"$work/suite.log" 2>&1 || {
    tail -20 "$work/suite.log" >&2
    echo "coverage: the suite failed; coverage of a red suite is meaningless" >&2
    exit 1
}

( cd "$work" && gcov -n -o src src/*.c ) >"$work/gcov.txt" 2>/dev/null || true

pct=$(python3 - "$work/gcov.txt" <<'PY'
import re, sys
SKIP = {"bgzf.c", "kstring.c"}          # vendored: htslib, klib
t = open(sys.argv[1]).read()
rows = re.findall(r"File '([^']+)'\nLines executed:([0-9.]+)% of (\d+)", t)
rows = [(f, float(p), int(n)) for f, p, n in rows
        if f.startswith("src/") and f.rsplit("/", 1)[-1] not in SKIP]
if not rows:
    sys.exit("coverage: gcov produced no rows for src/")
tot = sum(n for _, _, n in rows)
cov = sum(p / 100 * n for _, p, n in rows)
for f, p, n in sorted(rows, key=lambda r: r[1]):
    print(f"  {f:<24}{n:>7}{p:>8.1f}%", file=sys.stderr)
print(f"  {'TOTAL':<24}{tot:>7}{100*cov/tot:>8.1f}%", file=sys.stderr)
print(f"{100*cov/tot:.1f}")
PY
)

echo "coverage: ${pct}% of executable lines (gcov, line coverage: test/run.sh)"

badge="$here/docs/coverage.json"
case "$mode" in
  --print) exit 0 ;;
  --check)
      [ -f "$badge" ] || { echo "coverage: $badge is missing; run scripts/coverage.sh" >&2; exit 1; }
      old=$(sed -n 's/.*"message" *: *"\([0-9.]*\)%".*/\1/p' "$badge")
      awk -v a="$old" -v b="$pct" -v t="$TOLERANCE" 'BEGIN{d=a-b; if(d<0)d=-d; exit !(d>t)}' \
          && { echo "coverage: badge says ${old}% but the suite measures ${pct}% (> ${TOLERANCE} points); run scripts/coverage.sh and commit docs/coverage.json" >&2; exit 1; }
      echo "coverage: badge (${old}%) is within ${TOLERANCE} points"
      exit 0 ;;
esac

colour=$(awk -v p="$pct" 'BEGIN{
  print (p>=90)?"brightgreen":(p>=80)?"green":(p>=70)?"yellowgreen":(p>=60)?"yellow":(p>=50)?"orange":"red"}')
cat > "$badge" <<JSON
{
  "schemaVersion": 1,
  "label": "coverage",
  "message": "${pct}%",
  "color": "${colour}"
}
JSON
echo "coverage: wrote docs/coverage.json (${pct}%, ${colour})"
