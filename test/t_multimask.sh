#!/bin/bash
## The multi-mask kernel: identical to summarize1(), and faster.
##
## summarize1() takes ONE mask per call, so N masks walk the query N times, and
## about half the per-mask cost was that walk. The kernel reads the query once
## into a bitmap of covered rows, then measures each mask against it 64 rows at
## a time. Against the full 1359-set TFBS knowledgebase on a 29.4M-row
## methylome: 198 s -> 7.1 s for a fmt3 query, 111 s -> 5.8 s for a fmt6 one,
## byte-identical either way, peak memory unchanged at 45 MB. The fmt6 path is
## the faster of the two because it never reads a value: every count is a
## popcount of two or three bitmaps.
##
## Word-wise is half the point. A first version walked runs row by row and was
## only 2x, and at 128 overlapping sets it was SLOWER than the loop it replaced.
##
## Requested by methscope, which ships its own copy of the arithmetic and asked
## for the mask side to be an ENUMERATOR rather than records, because it walks
## membership runs inside a model file and materialising a .cm cost it 0.86 s
## per set. probe_multi.c drives the kernel directly and compares it against
## summarize1() mask by mask -- bit for bit, since within a mask both add the
## same doubles in the same row order.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/.." && pwd)
cfg="$root/yame-config"

[ -f "$root/libyame.a" ] && [ -x "$cfg" ] ||
  { echo "skip: no libyame.a for the kernel probe" >&2; exit 0; }

d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
cd "$d"

## ${CC:-cc}: the coverage run instruments libyame.a, and a probe built without
## the same flags cannot link against it.
${CC:-cc} -O1 -g -std=gnu99 $("$cfg" --cflags) -o probe_multi "$here/probe_multi.c" \
  $("$cfg" --libs) 2>cc.err || { echo "probe_multi did not build"; cat cc.err; exit 1; }
./probe_multi || { echo "the kernel disagrees with summarize1"; exit 1; }

## The ENUMERATOR path separately. It shares the arithmetic with the record path
## but not the route in, and it is the one that must survive a bank's slot count
## without allocating per slot: 40 masks x 8 states is 320 accumulators here.
${CC:-cc} -O1 -g -std=gnu99 $("$cfg" --cflags) -o probe_runs "$here/probe_runs.c" \
  $("$cfg" --libs) 2>cc2.err || { echo "probe_runs did not build"; cat cc2.err; exit 1; }
./probe_runs || { echo "the enumerator path disagrees with the record path"; exit 1; }

## ---- and through the CLI, where the numbers must not move -----------------
## A mask of several records against one query: the printed table has to be the
## same whether the kernel ran or the per-mask path did.
/usr/bin/awk 'BEGIN { for (i = 0; i < 4000; i++) {
    if (i % 5 == 0) print 0 "\t" 0
    else print (i % 9) "\t" (9 - (i % 9)) } }' > mu.txt
"$YAME" pack -f m mu.txt > q.cg

## three masks with different densities, so one is dense and one is sparse
for spec in "3 1" "7 3" "101 40"; do
  set -- $spec
  /usr/bin/awk -v e="$1" -v r="$2" 'BEGIN { for (i = 0; i < 4000; i++)
      print (i % e < r) ? 1 : 0 }' > m$1.txt
  "$YAME" pack -f b m$1.txt > m$1.cg
done
cat m3.cg m7.cg m101.cg > multi.cm
printf 'dense\nmid\nsparse\n' > mnames.txt
"$YAME" index -s mnames.txt multi.cm >/dev/null 2>&1

## Both paths: -M loads the masks, the default streams them, and the kernel is
## used either way. The two tables must be identical, and identical to what the
## per-mask path produced before the kernel existed.
"$YAME" summary -M -m multi.cm q.cg > got.txt 2>/dev/null
"$YAME" summary -m multi.cm q.cg > got_stream.txt 2>/dev/null
diff got.txt got_stream.txt ||
  { echo "the one-pass kernel (-M) and the per-mask path disagree"; exit 1; }
[ "$(tail -n +2 got.txt | grep -c .)" -eq 3 ] ||
  { echo "summary over a 3-record mask gave $(tail -n +2 got.txt | grep -c .) rows"; cat got.txt; exit 1; }

## each mask separately must give the same numbers as the combined run
for i in 3 7 101; do
  "$YAME" summary -m m$i.cg q.cg 2>/dev/null | tail -n +2 |
    cut -f5,6,7,8,10 >> sep.txt
done
tail -n +2 got.txt | cut -f5,6,7,8,10 > comb.txt
diff sep.txt comb.txt ||
  { echo "a combined mask run differs from the masks run one at a time"; exit 1; }

## ---- a mask the kernel DECLINES, among masks it takes ---------------------
## The fallback must be per mask. It used to be all-or-nothing: a flag was
## cleared and the whole list re-ran, so every mask already measured printed a
## SECOND time. `-M` over a binary mask followed by a state mask gave 5 rows
## where 4 were right, and only with -M, because the streamed branch always
## fell back one mask at a time. Both branches are checked, both orders.
/usr/bin/awk 'BEGIN { for (i = 0; i < 4000; i++)
    print (i % 3 == 0) ? "A" : ((i % 3 == 1) ? "B" : "C") }' > states.txt
"$YAME" pack -f 2 states.txt > state.cg 2>/dev/null
[ -s state.cg ] || { echo "fixture: could not pack a state mask"; exit 1; }

## 3 binary masks + 1 state mask of 3 terms = 3 + 3 = 6 rows, in any order
for order in "m3.cg m7.cg m101.cg state.cg" "state.cg m3.cg m7.cg m101.cg" \
             "m3.cg state.cg m7.cg m101.cg"; do
  cat $order > mixed.cm
  for opt in "-M" ""; do
    got=$("$YAME" summary $opt -m mixed.cm q.cg 2>/dev/null | tail -n +2 | grep -c .)
    [ "$got" -eq 6 ] ||
      { echo "mixed mask ($order, ${opt:-streamed}) gave $got rows, want 6"
        "$YAME" summary $opt -m mixed.cm q.cg 2>/dev/null | cut -f4; exit 1; }
  done
  ## and the two branches must agree with each other
  diff <("$YAME" summary -M -m mixed.cm q.cg 2>/dev/null) \
       <("$YAME" summary -m mixed.cm q.cg 2>/dev/null) ||
    { echo "mixed mask ($order): -M and streamed disagree"; exit 1; }
done

## ---- several query FILES with -M ------------------------------------------
## The masks are loaded once, before the query-file loop, but were freed INSIDE
## it, leaving the array dangling while its count stayed set. So
## `summary -M -m masks.cm a.cg b.cg` read freed records and printed nothing at
## all for b.cg -- a silent wrong answer, present since -M existed. Without -M
## it was always right, which is why nobody saw it.
cp q.cg q2.cg
for opt in "-M" ""; do
  rows=$("$YAME" summary $opt -m multi.cm q.cg q2.cg 2>/dev/null | tail -n +2 | grep -c .)
  [ "$rows" -eq 6 ] ||
    { echo "two query files (${opt:-streamed}) gave $rows rows, want 6"
      "$YAME" summary $opt -m multi.cm q.cg q2.cg 2>/dev/null | cut -f1,4; exit 1; }
  files=$("$YAME" summary $opt -m multi.cm q.cg q2.cg 2>/dev/null |
          tail -n +2 | cut -f1 | sort -u | grep -c .)
  [ "$files" -eq 2 ] ||
    { echo "two query files (${opt:-streamed}): only $files named in the output"; exit 1; }
done

## ---- the fast path must actually BE taken -----------------------------------
## Both paths give identical numbers, which is the point, so nothing else here
## can tell whether the kernel ran. A guard added during review rejected almost
## every mask and sent the whole knowledgebase down the fallback: the output
## stayed correct, all 29 tests passed, and only a stopwatch noticed -- 6 s
## became 107 s. YAME_SUMMARY_PATH reports the split so this can be asserted.
for opt in "-M" ""; do
  p=$(YAME_SUMMARY_PATH=1 "$YAME" summary $opt -m multi.cm q.cg 2>&1 >/dev/null |
      grep 'one-pass' | tail -1)
  [ -n "$p" ] || { echo "YAME_SUMMARY_PATH printed nothing (${opt:-streamed})"; exit 1; }
  k=$(printf '%s\n' "$p" | sed 's/.*one-pass \([0-9]*\).*/\1/')
  b=$(printf '%s\n' "$p" | sed 's/.*per-mask \([0-9]*\).*/\1/')
  [ "$k" -eq 3 ] && [ "$b" -eq 0 ] ||
    { echo "3 binary masks (${opt:-streamed}) took one-pass=$k per-mask=$b, want 3 and 0"
      exit 1; }
done
## a state mask must fall back, and say so
p=$(YAME_SUMMARY_PATH=1 "$YAME" summary -m state.cg q.cg 2>&1 >/dev/null |
    grep 'one-pass' | tail -1)
printf '%s\n' "$p" | grep 'one-pass 0, per-mask 1' >/dev/null ||
  { echo "a state mask did not fall back: $p"; exit 1; }

echo "ok: t_multimask"
