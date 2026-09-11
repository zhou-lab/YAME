#!/bin/bash
## Everything that needs COORDINATES: a format-7 .cr stream, the row finder
## that maps chrom+position to a row, rowsub's -l/-L/-R selection, and
## hprint's region and whole-genome views. These share one machinery and are
## the least exercised part of the tree.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
cd "$d"

## ---- a small reference: 3 chroms, 30 CpGs, positions 100 apart -----------
awk 'BEGIN {
  n = 0
  split("chr1 chr2 chr3", c, " ")
  for (k = 1; k <= 3; k++)
    for (i = 0; i < 10; i++) {
      beg = 100 + i * 100
      printf "%s\t%d\t%d\n", c[k], beg, beg + 2
      n++
    } }' > ref.bed
"$YAME" pack -f r ref.bed > ref.cr
[ "$("$YAME" info ref.cr 2>/dev/null | tail -1 | cut -f4)" -eq 30 ] ||
  { echo "reference is not 30 rows"; exit 1; }

## two samples over the same 30 rows
for s in 1 2; do
  awk -v s=$s 'BEGIN { for (i = 0; i < 30; i++) print ((i * (s + 1)) % 11) "\t" (10 - ((i * (s + 1)) % 11)) }' > s$s.txt
  "$YAME" pack -f m s$s.txt > s$s.cg
done
cat s1.cg s2.cg > two.cg
printf 'alpha\nbeta\n' > names.txt
"$YAME" index -s names.txt two.cg

## ---- 1. unpack -R prepends coordinates, in each -r mode -------------------
"$YAME" unpack -R ref.cr -f 1 s1.cg 2>/dev/null > withco.txt
[ "$(wc -l < withco.txt)" -eq 30 ] || { echo "-R changed the row count"; exit 1; }
[ "$(head -1 withco.txt | cut -f1)" = "chr1" ] || { echo "-R did not prepend the chrom"; head -1 withco.txt; exit 1; }
[ "$(head -1 withco.txt | cut -f2)" = "100" ] || { echo "-R begin is wrong"; head -1 withco.txt; exit 1; }
## the value column still matches the plain unpack
paste <("$YAME" unpack -f 1 s1.cg 2>/dev/null) <(cut -f4 withco.txt) |
  awk -F'\t' '$1 != $2 { bad++ } END { exit bad > 0 }' ||
  { echo "-R changed the values"; exit 1; }
for r in 0 1 2; do
  [ "$("$YAME" unpack -R ref.cr -r $r -f 1 s1.cg 2>/dev/null | wc -l)" -eq 30 ] ||
    { echo "unpack -R -r $r lost rows"; exit 1; }
done
## -C names the columns
"$YAME" unpack -C -R ref.cr -f 1 -a two.cg 2>/dev/null | head -1 > hdr.txt
grep -q 'alpha' hdr.txt || { echo "-C header does not name the samples"; cat hdr.txt; exit 1; }

## ---- 2. rowsub -L: pick rows by coordinate --------------------------------
## One chrm_beg1 per line, 1-based begin. The .cr begins at 100 (0-based), so
## the 1-based begin of the first CpG of chr2 is 101.
## chr3's CpGs are at 0-based 100..1000, so 1-based 101..1001: chr3_1001 is
## the LAST row of the file, which is the one worth naming here.
printf 'chr2_101\nchr2_301\nchr3_1001\n' > want.txt
"$YAME" rowsub -R ref.cr -L want.txt s1.cg 2>/dev/null > byco.cg
"$YAME" unpack -f -1 byco.cg 2>/dev/null > byco.txt
[ "$(wc -l < byco.txt)" -eq 3 ] || { echo "-L returned $(wc -l < byco.txt) rows, want 3"; exit 1; }
## rows 10, 12 and 29 of the original (0-based), in the order asked for
{ sed -n '11p' s1.txt; sed -n '13p' s1.txt; sed -n '30p' s1.txt; } > byco.want
diff byco.want byco.txt || { echo "-L picked the wrong rows"; exit 1; }

## a coordinate the reference does not carry
printf 'chr2_999999\n' > bad.txt
if "$YAME" rowsub -R ref.cr -L bad.txt s1.cg >/dev/null 2>err.txt; then
  [ -s err.txt ] || { echo "-L accepted an unknown coordinate silently"; exit 1; }
fi

## ---- 3. rowsub -l: pick rows by 1-based index, order preserved ------------
printf '5\n1\n30\n' > idx.txt
"$YAME" rowsub -l idx.txt s1.cg 2>/dev/null | "$YAME" unpack -f -1 - 2>/dev/null > byidx.txt
{ sed -n '5p' s1.txt; sed -n '1p' s1.txt; sed -n '30p' s1.txt; } > byidx.want
diff byidx.want byidx.txt || { echo "-l did not preserve the order asked for"; exit 1; }
## out of range
printf '31\n' > oor.txt
if "$YAME" rowsub -l oor.txt s1.cg >/dev/null 2>err2.txt; then
  [ -s err2.txt ] || { echo "-l accepted an out-of-range index silently"; exit 1; }
fi

## ---- 4. hprint region view -----------------------------------------------
## rows = samples, columns = CpG sites in the region.
"$YAME" hprint -c -R ref.cr -r chr2 two.cg 2>/dev/null > reg.txt
[ "$(wc -l < reg.txt)" -ge 2 ] || { echo "region view printed fewer than 2 sample rows"; cat reg.txt; exit 1; }
"$YAME" hprint -c -g -R ref.cr -r chr2:1-1000 two.cg 2>/dev/null > reg2.txt
[ -s reg2.txt ] || { echo "a chr:beg-end region printed nothing"; exit 1; }
## A region the reference has no CpGs in is an explicit error, not an empty
## view: it says so and exits non-zero, with nothing on stdout. Encoded here
## because silence would be the worse answer and a change either way should
## be deliberate.
if "$YAME" hprint -c -R ref.cr -r chr1:50000-60000 two.cg > empty.out 2> empty.err; then
  echo "an empty region exited 0"; exit 1
fi
[ ! -s empty.out ] || { echo "an empty region still printed a view"; exit 1; }
grep -q 'chr1:50000-60000' empty.err ||
  { echo "the empty-region message does not name the region"; cat empty.err; exit 1; }
## a chromosome the reference does not have at all
if "$YAME" hprint -c -R ref.cr -r chrZ two.cg >/dev/null 2>nochr.err; then
  echo "an unknown chromosome exited 0"; exit 1
fi
[ -s nochr.err ] || { echo "an unknown chromosome failed silently"; exit 1; }
## label column and ruler
"$YAME" hprint -c -R ref.cr -r chr2 -l 12 -t 5 two.cg 2>/dev/null > lab.txt
grep -q 'alpha' lab.txt || { echo "-l did not print sample labels"; head -3 lab.txt; exit 1; }
## Width: a region wider than -w is window-averaged rather than clipped, and
## the title says so. Compare a DATA row, not line 1 -- the narrow view's
## title is the LONGER of the two because it appends the window size.
wide=$("$YAME" hprint -c -g -R ref.cr -r chr2 -w 100 two.cg 2>/dev/null | sed -n 2p | wc -c)
narrow=$("$YAME" hprint -c -g -R ref.cr -r chr2 -w 4 two.cg 2>/dev/null | sed -n 2p | wc -c)
[ "$narrow" -lt "$wide" ] ||
  { echo "-w 4 did not narrow the data row ($narrow vs $wide)"; exit 1; }
"$YAME" hprint -c -g -R ref.cr -r chr2 -w 4 two.cg 2>/dev/null | head -1 | grep -q 'win=' ||
  { echo "a window-averaged view does not say so in its title"; exit 1; }

## A format-3 region is streamed compressed; every other format is inflated
## and drawn by a different printer, and only that one does the windowing. So
## the windowed region view needs a non-fmt3 record to be reached at all.
"$YAME" unpack -f 1 s1.cg 2>/dev/null > b1.txt
"$YAME" unpack -f 1 s2.cg 2>/dev/null > b2.txt
"$YAME" pack -f n b1.txt > b1.cg
"$YAME" pack -f n b2.txt > b2.cg
cat b1.cg b2.cg > beta.cg
"$YAME" index -s names.txt beta.cg
for extra in "" "-g"; do
  "$YAME" hprint -c $extra -R ref.cr -r chr2 -w 3 beta.cg 2>/dev/null > bw.txt
  [ -s bw.txt ] || { echo "windowed fmt4 region view ($extra) printed nothing"; exit 1; }
  head -1 bw.txt | grep -q 'win=' ||
    { echo "windowed fmt4 region ($extra) did not say win= in its title"; head -1 bw.txt; exit 1; }
done
## and a format-6 record through the same printer
"$YAME" binarize s1.cg > q6.cg 2>/dev/null
"$YAME" hprint -c -R ref.cr -r chr2 -w 3 q6.cg >/dev/null 2>&1 ||
  { echo "windowed fmt6 region view exited non-zero"; exit 1; }

## ---- 5. hprint whole-genome view ------------------------------------------
## One column per CpG window across all chroms; -R without -r.
"$YAME" hprint -c -R ref.cr two.cg 2>/dev/null > wg.txt
[ -s wg.txt ] || { echo "whole-genome view printed nothing"; exit 1; }
"$YAME" hprint -c -g -R ref.cr -w 10 two.cg 2>/dev/null > wg2.txt
[ -s wg2.txt ] || { echo "windowed whole-genome view printed nothing"; exit 1; }
## colour on (the default) must also run -- it is a different code path
"$YAME" hprint -R ref.cr -r chr2 two.cg >/dev/null 2>&1 ||
  { echo "the coloured region view exited non-zero"; exit 1; }

## ---- 6. a reference whose row count does not match the data ---------------
awk 'BEGIN { for (i = 0; i < 5; i++) printf "chr1\t%d\t%d\n", 10 + i, 12 + i }' |
  "$YAME" pack -f r - > short.cr
if "$YAME" unpack -R short.cr s1.cg >/dev/null 2>err3.txt; then
  echo "a reference of the wrong length was accepted"; exit 1
fi
[ -s err3.txt ] || { echo "mismatched reference failed silently"; exit 1; }
