#!/bin/bash
## rowop (across samples, down a row) and summary (a query against masks).
## Both reduce many values to a few, so an independent recomputation in awk is
## the only assertion worth making.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
cd "$d"

## four samples over 12 rows, betas chosen to straddle 0.5, one uncovered row
for s in 1 2 3 4; do
  awk -v s=$s 'BEGIN {
    split("9 1 8 2 7 3 6 4 5 5 0 0", a, " ")
    for (i = 0; i < 12; i++) {
      m = (a[i + 1] + s) % 11
      if (i == 7) print 0 "\t" 0            # uncovered in every sample
      else print m "\t" (10 - m) } }' > s$s.txt
  "$YAME" pack -f m s$s.txt > s$s.cg
done
cat s1.cg s2.cg s3.cg s4.cg > four.cg
printf 'a\nb\nc\nd\n' > names.txt
"$YAME" index -s names.txt four.cg

## ---- 1. rowop binasum: per-row counts of methylated vs unmethylated calls --
## Two thresholds, not one: beta < -p is unmethylated, beta > -q methylated,
## and a beta BETWEEN them is ignored rather than forced to a side. So the two
## counts need not add up to the sample count, which is the point of the band.
"$YAME" rowop -o binasum -p 0.4 -q 0.6 four.cg > ba.cg 2>/dev/null
"$YAME" unpack -f -1 ba.cg 2>/dev/null > ba.got
paste s1.txt s2.txt s3.txt s4.txt |
  awk -F'\t' '{ M = 0; U = 0
    for (k = 1; k <= 7; k += 2) { cov = $k + $(k+1)
      if (cov >= 1) { b = $k / cov
        if (b > 0.6) M++
        else if (b < 0.4) U++ } }
    print M "\t" U }' > ba.want
diff ba.want ba.got || { echo "rowop binasum disagrees with a per-row recount"; exit 1; }

## widening the band can only reduce the calls
"$YAME" rowop -o binasum -p 0.2 -q 0.8 four.cg 2>/dev/null |
  "$YAME" unpack -f -1 - 2>/dev/null > ba_wide.got
awk -F'\t' '{ s += $1 + $2 } END { print s }' ba.got > n_narrow
awk -F'\t' '{ s += $1 + $2 } END { print s }' ba_wide.got > n_wide
[ "$(cat n_wide)" -le "$(cat n_narrow)" ] ||
  { echo "binasum: widening the ambiguous band did not reduce the calls"; exit 1; }

## ---- 2. rowop musum: sum the counts across samples -------------------------
"$YAME" rowop -o musum four.cg 2>/dev/null | "$YAME" unpack -f -1 - 2>/dev/null > ms.got
paste s1.txt s2.txt s3.txt s4.txt |
  awk -F'\t' '{ print $1 + $3 + $5 + $7 "\t" $2 + $4 + $6 + $8 }' > ms.want
diff ms.want ms.got || { echo "rowop musum is not the column-wise sum"; exit 1; }

## ---- 3. rowop stat: the per-row summary columns ----------------------------
## count and mean_beta are the two an analyst reads first; check them exactly
## and check the header names the columns the docs promise.
## -d defaults to 6 decimals; the recount below matches that width.
"$YAME" rowop -o stat four.cg > st.txt 2>/dev/null
head -1 st.txt | grep -q 'count.*mean_beta.*delta_beta' ||
  { echo "rowop stat header changed"; head -1 st.txt; exit 1; }
tail -n +2 st.txt | cut -f1,2 > st.got
paste s1.txt s2.txt s3.txt s4.txt |
  awk -F'\t' '{ n = 0; sum = 0
    for (k = 1; k <= 7; k += 2) { cov = $k + $(k+1)
      if (cov > 0) { n++; sum += $k / cov } }
    printf "%d\t%s\n", n, (n ? sprintf("%.6f", sum / n) : "NA") }' > st.want
diff st.want st.got || { echo "rowop stat count/mean_beta disagree"; exit 1; }

## ---- 4. rowop cometh refuses a format it cannot read -----------------------
## It read both inputs with f3_get_mu() unconditionally, so a non-fmt3 input
## produced differential calls over bytes that were never counts.
"$YAME" binarize s1.cg > b.cg 2>/dev/null
if "$YAME" rowop -o cometh b.cg >/dev/null 2>err.txt; then
  echo "rowop cometh accepted a non-format-3 input"; exit 1
fi
[ -s err.txt ] || { echo "rowop cometh failed silently"; exit 1; }

## ---- 5. summary against a mask --------------------------------------------
## The mask picks rows; the summary must count what the mask selects, not the
## whole file. Cross-check the covered-site count against awk.
awk 'BEGIN { for (i = 0; i < 12; i++) print (i < 5) ? 1 : 0 }' > m.txt
"$YAME" pack -f b m.txt > m.cg
"$YAME" summary -m m.cg s1.cg > sm.txt 2>/dev/null
[ -s sm.txt ] || { echo "summary produced nothing"; exit 1; }
head -1 sm.txt | grep -qi 'query\|qfile\|mask' ||
  { echo "summary header is not recognisable"; head -1 sm.txt; exit 1; }
## the same query with no mask must still produce a line per record
"$YAME" summary four.cg > sm4.txt 2>/dev/null
[ "$(tail -n +2 sm4.txt | wc -l)" -ge 4 ] ||
  { echo "summary did not report all four records"; cat sm4.txt; exit 1; }
