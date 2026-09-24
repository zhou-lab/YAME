#!/bin/bash
## dsample: keep N covered sites per sample at random, mask the rest.
##
## t_ops.sh counted records and checked the format; nothing checked WHICH
## sites survive, nothing ran -o (so the output index and its replicate names
## had no test at all), and three bugs lived in the gaps: `dsample in out`
## ignored `out` and wrote stdout, -b did nothing once N reached the
## coverage, and a store whose records differ in length overflowed the
## sampling buffers. Each is pinned below.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
cd "$d"

## two samples, 50 rows, every third row uncovered: 33 covered sites each,
## and M encodes the row so a moved value cannot pass for a kept one
for s in 1 2; do
  awk -v s=$s 'BEGIN { for (i = 0; i < 50; i++) print (i % 3 ? i + s : 0) "\t" (i % 3 ? 50 - i : 0) }' > s$s.txt
  "$YAME" pack -f m s$s.txt s$s.cg
done
cat s1.cg s2.cg > two.cg
printf 'alpha\nbeta\n' > nm.txt
"$YAME" index -s nm.txt two.cg

## ---- 1. fmt3: exactly N sites survive, each unchanged ----------------------
"$YAME" dsample -s 7 -N 10 s1.cg 2>/dev/null | "$YAME" unpack -f -1 - 2>/dev/null > ds.txt
[ "$(wc -l < ds.txt)" -eq 50 ] || { echo "dsample changed the row count"; exit 1; }
[ "$(awk '$1 + $2 > 0' ds.txt | wc -l)" -eq 10 ] ||
  { echo "-N 10 kept $(awk '$1 + $2 > 0' ds.txt | wc -l) covered sites"; exit 1; }
## a kept row carries the input's M/U; a dropped or uncovered row is 0 0
paste s1.txt ds.txt | awk -F'\t' '$3 + $4 > 0 && ($1 != $3 || $2 != $4) { bad = 1 }
  $1 + $2 == 0 && $3 + $4 > 0 { bad = 1 } END { exit bad }' ||
  { echo "a kept site's M/U changed, or an uncovered site gained coverage"; exit 1; }

## the seed decides the sites: same seed, same file; another seed, others
"$YAME" dsample -s 7 -N 10 s1.cg 2>/dev/null | "$YAME" unpack -f -1 - 2>/dev/null > ds2.txt
diff -q ds.txt ds2.txt >/dev/null || { echo "-s 7 twice gave two different samples"; exit 1; }
"$YAME" dsample -s 8 -N 10 s1.cg 2>/dev/null | "$YAME" unpack -f -1 - 2>/dev/null > ds3.txt
! diff -q ds.txt ds3.txt >/dev/null || { echo "-s 7 and -s 8 kept the same sites"; exit 1; }

## N at or past the coverage keeps everything as it was
"$YAME" dsample -s 7 -N 33 s1.cg 2>/dev/null | "$YAME" unpack -f -1 - 2>/dev/null > all.txt
diff s1.txt all.txt || { echo "-N 33 on 33 covered sites changed the data"; exit 1; }

## ---- 2. -b: the kept sites become one read, methylated or not --------------
for n in 10 1000; do                    # below and past the coverage
  "$YAME" dsample -s 7 -N $n -b s1.cg 2>/dev/null | "$YAME" unpack -f -1 - 2>/dev/null > b.txt
  want=$([ $n -lt 33 ] && echo $n || echo 33)
  [ "$(awk '$1 + $2 > 0' b.txt | wc -l)" -eq $want ] ||
    { echo "-b -N $n kept the wrong number of sites"; exit 1; }
  [ "$(awk '$1 + $2 > 0 && !(($1 == 1 && $2 == 0) || ($1 == 0 && $2 == 1))' b.txt | wc -l)" -eq 0 ] ||
    { echo "-b -N $n left a site that is not a single read"; awk '$1 + $2 > 1' b.txt | head -3; exit 1; }
done

## ---- 3. fmt6: N universe rows survive, their set bits intact ---------------
"$YAME" binarize s1.cg > b6.cg 2>/dev/null
"$YAME" unpack -f -1 b6.cg 2>/dev/null > b6.txt           # set <tab> universe
"$YAME" dsample -s 7 -N 10 b6.cg 2>/dev/null | "$YAME" unpack -f -1 - 2>/dev/null > d6.txt
[ "$(awk '$2 == 1' d6.txt | wc -l)" -eq 10 ] || { echo "fmt6 -N 10 kept the wrong universe"; exit 1; }
paste b6.txt d6.txt | awk -F'\t' '$4 == 1 && ($2 != 1 || $1 != $3) { bad = 1 } END { exit bad }' ||
  { echo "fmt6: a kept row is outside the input universe or changed its set bit"; exit 1; }

## ---- 4. -o writes an index; replicates are named from the input's ---------
"$YAME" dsample -s 7 -N 5 -r 2 -p rep -o o.cg two.cg 2>/dev/null
[ "$(cut -f1 o.cg.idx | paste -sd' ' -)" = "alpha-rep-0 alpha-rep-1 beta-rep-0 beta-rep-1" ] ||
  { echo "-r 2 -p rep named $(cut -f1 o.cg.idx | paste -sd' ' -)"; exit 1; }
## the index points at the right records: each replicate is its own draw
"$YAME" subset o.cg alpha-rep-1 2>/dev/null | "$YAME" unpack -f -1 - 2>/dev/null > r1.txt
paste s1.txt r1.txt | awk -F'\t' '$3 + $4 > 0 && ($1 != $3 || $2 != $4) { bad = 1 } END { exit bad }' ||
  { echo "the index entry alpha-rep-1 does not point at a sample of alpha"; exit 1; }
"$YAME" subset o.cg alpha-rep-0 2>/dev/null | "$YAME" unpack -f -1 - 2>/dev/null > r0.txt
! diff -q r0.txt r1.txt >/dev/null || { echo "two replicates drew the same sites"; exit 1; }
"$YAME" dsample -s 7 -N 5 -r 2 -o o2.cg two.cg 2>/dev/null
[ "$(cut -f1 o2.cg.idx | paste -sd' ' -)" = "alpha-0 alpha-1 beta-0 beta-1" ] ||
  { echo "-r 2 without -p named $(cut -f1 o2.cg.idx | paste -sd' ' -)"; exit 1; }
"$YAME" dsample -s 7 -N 5 -o o1.cg two.cg 2>/dev/null
[ "$(cut -f1 o1.cg.idx | paste -sd' ' -)" = "alpha beta" ] ||
  { echo "-r 1 renamed the samples: $(cut -f1 o1.cg.idx | paste -sd' ' -)"; exit 1; }
## no input index: 0-based numbers stand in for the names
cat s1.cg s2.cg > noidx.cg
"$YAME" dsample -s 7 -N 5 -r 2 -o o3.cg noidx.cg 2>/dev/null
[ "$(cut -f1 o3.cg.idx | paste -sd' ' -)" = "0-0 0-1 1-0 1-1" ] ||
  { echo "unindexed input named $(cut -f1 o3.cg.idx | paste -sd' ' -)"; exit 1; }

## ---- 5. the output may be named last; two outputs, or three files, are not -
"$YAME" dsample -s 7 -N 5 two.cg pos.cg > stdout.cg 2>/dev/null
[ -s pos.cg ] && [ -s pos.cg.idx ] || { echo "dsample in.cg out.cg did not write out.cg"; exit 1; }
cmp -s pos.cg o1.cg || { echo "dsample in out and -o out disagree"; exit 1; }
if "$YAME" dsample -o a.cg two.cg b.cg >/dev/null 2>err.txt; then
  echo "-o a.cg with a second output b.cg was accepted"; exit 1
fi
grep -q 'Two outputs' err.txt || { echo "-o plus a positional output failed without saying why"; exit 1; }
if "$YAME" dsample two.cg a.cg b.cg >/dev/null 2>err.txt; then
  echo "dsample accepted three files"; exit 1
fi

## ---- 6. records of different lengths in one store --------------------------
awk 'BEGIN { for (i = 0; i < 500; i++) print i + 1 "\t1" }' | "$YAME" pack -f m - long.cg
cat s1.cg long.cg > mixed.cg
"$YAME" dsample -s 7 -N 5 -o dm.cg mixed.cg 2>/dev/null   # unindexed: records 0 and 1
[ "$("$YAME" info dm.cg 2>/dev/null | tail -n +2 | cut -f4 | paste -sd' ' -)" = "50 500" ] ||
  { echo "a mixed-length store did not keep its record lengths"; "$YAME" info dm.cg; exit 1; }
[ "$("$YAME" subset dm.cg 1 2>/dev/null | "$YAME" unpack -f -1 - 2>/dev/null | awk '$1 + $2 > 0' | wc -l)" -eq 5 ] ||
  { echo "the longer record did not keep -N 5 sites"; exit 1; }

## ---- 7. only formats 3 and 6 can be downsampled ----------------------------
printf '1\n0\n1\n' | "$YAME" pack -f b - f0.cg
if "$YAME" dsample -N 1 f0.cg >/dev/null 2>err.txt; then
  echo "dsample accepted a format 0 file"; exit 1
fi
grep -q 'only 3 and 6' err.txt || { echo "a format 0 input failed without saying why"; exit 1; }
echo "ok: t_dsample"
