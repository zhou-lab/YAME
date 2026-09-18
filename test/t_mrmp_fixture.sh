#!/bin/bash
## The .cm-exported-from-.mrmp fixture: one membership, two views, one answer.
##
## methscope feeds YAME's kernel membership RUNS walked out of a .mrmp; yame
## reads the exported .cm as RECORDS. Both come from one build in methscope,
## and nothing else enforces that the two views agree. test/fixtures/mrmp_export
## holds that build (6 cells x 400 CpGs, 7 KB): a flat set of 2 patterns, and a
## bank of 18 overlapping sets, one .cm record each, plus the per-pattern means
## methscope's own mrmp-summary printed for them. Dropped by methscope on
## 2026-09-18; build.sh to regenerate it lives with methscope.
##
## Two checks. probe_fixture.c decodes the .cm into runs and drives the
## enumerator path and the inverted index against the record path, bit for
## bit. Then `yame summary` over the same files is joined with methscope's
## table on (sample, set, pattern) and every beta must agree at yame's three
## printed decimals -- the agreement methscope's suite checks from its side.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/.." && pwd)
cfg="$root/yame-config"
fx="$here/fixtures/mrmp_export"

[ -f "$fx/bank.cm" ] || { echo "skip: no fixtures/mrmp_export in this tree" >&2; exit 0; }

d=$(mktemp -d); trap 'rm -rf "$d"' EXIT

## ---- runs against records, through the library ----------------------------
if [ -f "$root/libyame.a" ] && [ -x "$cfg" ]; then
  ${CC:-cc} -O1 -g -std=gnu99 $("$cfg" --cflags) -o "$d/probe_fixture" \
    "$here/probe_fixture.c" $("$cfg" --libs) 2>"$d/cc.err" ||
    { echo "probe_fixture did not build"; cat "$d/cc.err"; exit 1; }
  "$d/probe_fixture" "$fx/flat.cm" "$fx/ref.cg" ||
    { echo "flat: runs and records disagree"; exit 1; }
  "$d/probe_fixture" "$fx/bank.cm" "$fx/bankref.cg" ||
    { echo "bank: runs and records disagree"; exit 1; }
else
  echo "no libyame.a: skipping the library half" >&2
fi

## ---- yame summary against methscope's mrmp-summary --------------------------
## methscope prints (sample, set, pattern, beta) with beta at four decimals and
## no Pna row; yame prints one row per state with -T naming it "set-state".
## Joined on sample+set+pattern, every beta must agree within 5e-4, which is
## yame's third decimal.
## flat.cm is one record with no .idx, so yame names its set "1" where
## methscope's table says "root"; the 4th argument supplies that name.
check() {                       # <mask.cm> <query.cg> <expected.tsv> [setname]
  "$YAME" summary -T -m "$1" "$2" 2>/dev/null |
    /usr/bin/awk -F'\t' -v setname="${4:-}" 'NR > 1 && $4 !~ /Pna$/ {
      n = split($4, a, "-"); set = a[1]; for (i = 2; i < n; i++) set = set "-" a[i]
      if (setname != "") set = setname
      print $2 "\t" set "\t" a[n] "\t" $10 }' | sort > "$d/got.tsv"
  tail -n +2 "$3" | sort > "$d/want.tsv"
  [ "$(grep -c . "$d/got.tsv")" -eq "$(grep -c . "$d/want.tsv")" ] ||
    { echo "$1: yame printed $(grep -c . "$d/got.tsv") rows, methscope $(grep -c . "$d/want.tsv")"
      diff "$d/got.tsv" "$d/want.tsv" | head; exit 1; }
  paste "$d/got.tsv" "$d/want.tsv" |
    /usr/bin/awk -F'\t' '{
      if ($1 != $5 || $2 != $6 || $3 != $7) { print "key mismatch: " $0; bad = 1 }
      dv = $4 - $8; if (dv < 0) dv = -dv
      if (dv > 5e-4) { print "beta differs: " $0; bad = 1 } }
      END { exit bad }' ||
    { echo "$1: yame and mrmp-summary disagree"; exit 1; }
}
check "$fx/flat.cm" "$fx/ref.cg" "$fx/flat.expected.tsv" root
check "$fx/bank.cm" "$fx/bankref.cg" "$fx/bank.expected.tsv"

echo "ok: t_mrmp_fixture"
