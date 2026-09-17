#!/bin/bash
## The multi-mask kernel: identical to summarize1(), and faster.
##
## summarize1() takes ONE mask per call, so N masks walk the query N times.
## Measured before the kernel, on a 29.4M-row methylome of 9 samples: 1 mask
## 1.6 s, 4 masks 4.8 s, 16 masks 16.8 s. The walk decompresses; the
## accumulation is a few adds. A 178-set bank is about half a minute per cell.
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

## -M is what switches the kernel on: it needs every mask in hand at once.
## Without it masks are streamed and the per-mask path runs, so BOTH are
## checked here and both must give the same table.
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

echo "ok: t_multimask"
