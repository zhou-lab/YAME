#!/bin/bash
## rowop beyond the fmt3 sums t_rowop_summary.sh checks: binasum over binary
## (fmt0) and byte (fmt1) inputs, threaded and not; the refusals for mixed
## formats and lengths on both paths; -d, -t with an op that cannot thread,
## the positional output and -v.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
cd "$d"

## ---- 1. binasum over binary and byte records -----------------------------
## three fmt0 samples: 101010, 010101, 101010 -> per row M = ones, U = zeros
for s in 1 2 3; do
  awk -v s=$s 'BEGIN { for (i = 0; i < 6; i++) print ((i + s) % 2) }' | "$YAME" pack -f b - b$s.cg
done
cat b1.cg b2.cg b3.cg > b.cg
want=$(printf '2\t1 1\t2 2\t1 1\t2 2\t1 1\t2')
for t in 1 2; do
  [ "$("$YAME" rowop -t $t -o binasum b.cg 2>/dev/null | "$YAME" unpack -f -1 - 2>/dev/null | paste -sd' ' -)" = "$want" ] ||
    { echo "binasum over fmt0 (-t $t)"; exit 1; }
done
## fmt1: any value but 0 counts as methylated. 120120 and 201201
for s in 1 2; do
  awk -v s=$s 'BEGIN { for (i = 0; i < 6; i++) print (i + s) % 3 }' | "$YAME" pack -f c - c$s.cg
done
cat c1.cg c2.cg > c.cg
want=$(printf '2\t0 1\t1 1\t1 2\t0 1\t1 1\t1')
for t in 1 2; do
  [ "$("$YAME" rowop -t $t -o binasum c.cg 2>/dev/null | "$YAME" unpack -f -1 - 2>/dev/null | paste -sd' ' -)" = "$want" ] ||
    { echo "binasum over fmt1 (-t $t)"; exit 1; }
done

## ---- 2. records that cannot be summed together -----------------------------
awk 'BEGIN { for (i = 0; i < 6; i++) print i "\t1" }' | "$YAME" pack -f m - m6.cg
awk 'BEGIN { for (i = 0; i < 5; i++) print i "\t1" }' | "$YAME" pack -f m - m5.cg
cat m6.cg b1.cg > mix.cg                  # fmt3 then fmt0
cat m6.cg m5.cg > dim.cg                  # 6 rows then 5
for t in 1 2; do
  for pair in "mix:formats are inconsistent" "dim:imensions are inconsistent"; do
    f=${pair%%:*}; msg=${pair#*:}
    if "$YAME" rowop -t $t -o musum $f.cg >/dev/null 2>err.txt; then echo "musum -t $t took $f.cg"; exit 1; fi
    grep -q "$msg" err.txt || { echo "musum -t $t on $f.cg failed without saying '$msg'"; exit 1; }
  done
done
for op in musum stat; do                  # fmt3 only
  if "$YAME" rowop -t 2 -o $op b.cg >/dev/null 2>err.txt; then echo "$op took fmt0"; exit 1; fi
  grep -q 'format 0 unsupported' err.txt || { echo "$op on fmt0 failed without saying why"; exit 1; }
done
if "$YAME" rowop -o nope m6.cg >/dev/null 2>err.txt; then echo "an unknown op was accepted"; exit 1; fi
grep -q 'Unsupported operation: nope' err.txt || { echo "an unknown op failed without naming it"; exit 1; }

## ---- 3. options ---------------------------------------------------------------
cat m6.cg m6.cg > mm.cg
## -d sets stat's decimals; outside 0-9 is refused
[ "$("$YAME" rowop -o stat -d 2 mm.cg 2>/dev/null | sed -n 2p | cut -f2)" = "0.00" ] ||
  { echo "-d 2 did not print two decimals"; exit 1; }
if "$YAME" rowop -o stat -d 12 mm.cg >/dev/null 2>err.txt; then echo "-d 12 was accepted"; exit 1; fi
grep -q 'takes 0-9 decimals' err.txt || { echo "-d 12 failed without saying why"; exit 1; }
## -t on an op that has no threaded form says so and still runs
"$YAME" rowop -t 2 -o binstring mm.cg > bs.txt 2>err.txt
grep -q 'runs single-threaded' err.txt || { echo "-t on binstring did not say it runs single-threaded"; exit 1; }
[ -s bs.txt ] || { echo "binstring under -t printed nothing"; exit 1; }
## the output may be named last
"$YAME" rowop -o musum mm.cg out.cg 2>/dev/null
[ "$("$YAME" unpack -f -1 out.cg 2>/dev/null | head -1)" = "$(printf '0\t2')" ] ||
  { echo "rowop in out did not write the sum to out"; exit 1; }
## -v reports the threaded plan
"$YAME" rowop -v -t 2 -o musum mm.cg >/dev/null 2>v.txt
grep -q 'reduce .* (2 threads)' v.txt || { echo "-v -t 2 did not report the reduce"; cat v.txt; exit 1; }
echo "ok: t_rowop"
