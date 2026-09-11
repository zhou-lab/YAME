#!/bin/bash
## rowsub: row selection by block, by index and by mask. Its row indices are
## uint64_t everywhere and were parsed with atoi, so 2^32+10 wrapped to 10 and
## the command answered a question nobody asked, rc 0 and no message.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
cd "$d"

## 50 rows whose beta encodes the row number, so a wrong slice is obvious
awk 'BEGIN { for (i = 0; i < 50; i++) print i "\t" (50 - i) }' > all.txt
"$YAME" pack -f m all.txt > all.cg

## ---- 1. -B: a contiguous block by absolute row range ----------------------
"$YAME" rowsub -B 10_20 all.cg > blk.cg 2>/dev/null
"$YAME" unpack -f -1 blk.cg 2>/dev/null > blk.txt
sed -n '11,20p' all.txt > blk.want          # 0-based [10,20)
diff blk.want blk.txt || { echo "-B 10_20 is not rows [10,20)"; exit 1; }

## the whole file, and a single row
"$YAME" rowsub -B 0_50 all.cg 2>/dev/null | "$YAME" unpack -f -1 - 2>/dev/null > whole.txt
diff all.txt whole.txt || { echo "-B 0_50 is not the whole file"; exit 1; }
"$YAME" rowsub -B 7_8 all.cg 2>/dev/null | "$YAME" unpack -f -1 - 2>/dev/null > one.txt
[ "$(wc -l < one.txt)" -eq 1 ] || { echo "-B 7_8 did not return one row"; exit 1; }
diff <(sed -n 8p all.txt) one.txt || { echo "-B 7_8 returned the wrong row"; exit 1; }

## ---- 2. row indices are 64-bit -------------------------------------------
## 2^32 + 10 must not wrap to 10. Either it is refused, or it is treated as
## out of range -- what it may NOT do is quietly return rows 10 to 20.
if "$YAME" rowsub -B 4294967306_4294967316 all.cg > wrap.cg 2>err.txt; then
  "$YAME" unpack -f -1 wrap.cg 2>/dev/null > wrap.txt || true
  if diff -q <(sed -n '11,20p' all.txt) wrap.txt >/dev/null 2>&1; then
    echo "-B 2^32+10_2^32+20 wrapped to rows 10-20 and exited 0"; exit 1
  fi
else
  [ -s err.txt ] || { echo "an out-of-range -B failed silently"; exit 1; }
fi

## non-numeric must not mean row 0
if "$YAME" rowsub -B abc_def all.cg >/dev/null 2>err.txt; then
  echo "-B abc_def was accepted"; exit 1
fi
[ -s err.txt ] || { echo "-B abc_def failed silently"; exit 1; }

## ---- 3. -I: block index x block size --------------------------------------
## The block index is 0-BASED, so -I 1_10 is rows [10,20) -- the same span
## -B 10_20 gives, which is the check worth making.
"$YAME" rowsub -I 1_10 all.cg 2>/dev/null | "$YAME" unpack -f -1 - 2>/dev/null > idx.txt
diff blk.txt idx.txt || { echo "-I 1_10 and -B 10_20 disagree"; exit 1; }
"$YAME" rowsub -I 0_10 all.cg 2>/dev/null | "$YAME" unpack -f -1 - 2>/dev/null > idx0.txt
diff <(sed -n '1,10p' all.txt) idx0.txt || { echo "-I 0_10 is not the first block"; exit 1; }
## a zero block size has no meaning
if "$YAME" rowsub -I 1_0 all.cg >/dev/null 2>err.txt; then
  echo "-I with block size 0 was accepted"; exit 1
fi

## ---- 4. -m: a binary mask keeps the rows whose bit is set -----------------
awk 'BEGIN { for (i = 0; i < 50; i++) print (i % 5 == 0) ? 1 : 0 }' > m.txt
"$YAME" pack -f b m.txt > m.cg
"$YAME" rowsub -m m.cg all.cg 2>/dev/null | "$YAME" unpack -f -1 - 2>/dev/null > msk.txt
paste m.txt all.txt | awk -F'\t' '$1 == 1 { print $2 "\t" $3 }' > msk.want
diff msk.want msk.txt || { echo "-m did not keep exactly the set rows"; exit 1; }
[ "$(wc -l < msk.txt)" -eq 10 ] || { echo "-m kept $(wc -l < msk.txt) rows, want 10"; exit 1; }

## ---- 5. a mask of the wrong length is refused -----------------------------
awk 'BEGIN { for (i = 0; i < 20; i++) print 1 }' | "$YAME" pack -f b - > short.cg
if "$YAME" rowsub -m short.cg all.cg >/dev/null 2>err.txt; then
  echo "-m accepted a mask shorter than the data"; exit 1
fi
[ -s err.txt ] || { echo "a short mask failed silently"; exit 1; }
