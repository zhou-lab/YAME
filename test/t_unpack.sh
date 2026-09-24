#!/bin/bash
## unpack's selection and coordinate options: -H, -T, -C, -R with each -r mode,
## and the errors for an unindexed name, a bad -u and no input.
##
## Values encode their sample (M = 10*s + row), so a column from the wrong
## record cannot pass for the right one.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
cd "$d"

for s in 1 2 3; do
  awk -v s=$s 'BEGIN { for (i = 0; i < 4; i++) print s * 10 + i "\t" i }' > s$s.txt
  "$YAME" pack -f m s$s.txt s$s.cg
done
cat s1.cg s2.cg s3.cg > t.cg
printf 'a\nb\nc\n' > nm.txt
"$YAME" index -s nm.txt t.cg
printf 'chr1\t10\t12\nchr1\t20\t22\nchr2\t5\t7\nchr2\t9\t11\n' > r.bed
"$YAME" pack -f r r.bed r.cr

## ---- 1. which records ------------------------------------------------------
diff <(paste s1.txt s2.txt) <("$YAME" unpack -H 2 -f -1 t.cg 2>/dev/null) ||
  { echo "-H 2 is not the first two records"; exit 1; }
diff s3.txt <("$YAME" unpack -T 1 -f -1 t.cg 2>/dev/null) || { echo "-T 1 is not the last record"; exit 1; }
diff s1.txt <("$YAME" unpack -f -1 t.cg 2>/dev/null | tail -n +1) || { echo "no name is not the first record"; exit 1; }
[ "$("$YAME" unpack -C -H 2 -f 1 t.cg 2>/dev/null | head -1)" = "$(printf 'a\tb')" ] ||
  { echo "-C did not name the columns a b"; exit 1; }

## a record printing two columns gets two names, so the header lines up
[ "$("$YAME" unpack -C -H 2 -f -1 t.cg 2>/dev/null | head -1)" = "$(printf 'a_M\ta_U\tb_M\tb_U')" ] ||
  { echo "-C -f -1 did not name both M and U columns"; exit 1; }
"$YAME" binarize s1.cg > b6.cg 2>/dev/null
printf 'x\n' > b6.txt; "$YAME" index -s b6.txt b6.cg
[ "$("$YAME" unpack -C -f -1 b6.cg 2>/dev/null | head -1)" = "$(printf 'x_set\tx_uni')" ] ||
  { echo "-C -f -1 on format 6 did not name set and universe"; exit 1; }

## ---- 2. coordinates, in each -r mode ---------------------------------------
## 0: chrm beg0 end1 (as packed); 1: chrm beg0 end0; else: chrm_beg1
diff <(paste r.bed s2.txt) <("$YAME" unpack -R r.cr -f -1 t.cg b 2>/dev/null) ||
  { echo "-r 0 is not chrm beg0 end1"; exit 1; }
diff <(awk -F'\t' -v OFS='\t' '{ print $1, $2, $3 - 1 }' r.bed | paste - s3.txt) \
     <("$YAME" unpack -R r.cr -r 1 -f -1 t.cg c 2>/dev/null) || { echo "-r 1 is not chrm beg0 end0"; exit 1; }
[ "$("$YAME" unpack -C -R r.cr -r 2 -f 1 t.cg a 2>/dev/null | head -2 | paste -sd' ' -)" = \
  "$(printf 'chrm_beg1\ta chr1_11\t1.000')" ] || { echo "-r 2 is not chrm_beg1, or its header is wrong"; exit 1; }
[ "$("$YAME" unpack -C -R r.cr -r 1 -f 1 t.cg a 2>/dev/null | head -1)" = "$(printf 'chrm\tbeg0\tend0\ta')" ] ||
  { echo "-C with -r 1 did not name the coordinate columns"; exit 1; }

## ---- 3. what unpack refuses --------------------------------------------------
cat s1.cg s2.cg > noidx.cg
for args in "noidx.cg b" "-T 1 noidx.cg"; do
  if "$YAME" unpack -f -1 $args >/dev/null 2>err.txt; then echo "unpack $args worked without an index"; exit 1; fi
  grep -q 'needs indexing' err.txt || { echo "unpack $args failed without saying why"; exit 1; }
done
if "$YAME" unpack -u 3 -f -1 t.cg >/dev/null 2>err.txt; then echo "-u 3 was accepted"; exit 1; fi
grep -q 'can only be 1,2,4,6,8' err.txt || { echo "-u 3 failed without saying why"; exit 1; }
if "$YAME" unpack >/dev/null 2>err.txt; then echo "unpack with no input exited 0"; exit 1; fi
grep -q 'supply input file' err.txt || { echo "unpack with no input said nothing useful"; exit 1; }
echo "ok: t_unpack"
