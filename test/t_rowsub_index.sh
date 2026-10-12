#!/bin/bash
## rowsub through the row index: the .idx third column `yame index` writes
## for format-3 records. With it, rowsub inflates only the BGZF blocks that
## hold the rows asked for; without it, or when it no longer matches the
## file, every record is read in full. Either way the output must be the
## same bytes -- that is the whole contract, so most checks here are cmp.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
cd "$d"

## ---- fixtures ---------------------------------------------------------------
## Three format-3 records of 200000 rows, each several BGZF blocks long, with
## every token kind: zero runs, M/U <= 6 (1 byte), <= 126 (2 bytes) and
## large counts (8 bytes). "small" holds only 1-byte values, so a full decode
## infers a 1-byte unit where the row index decodes at 8: the bytes out must
## still agree. A format-0 record rides along to show a non-format-3 record
## is read in full in the same file.
N=200000
mk() {
  awk -v seed="$1" -v big="$2" -v N=$N 'BEGIN { srand(seed)
    for (i = 0; i < N; i++) {
      r = rand()
      if (r < 0.3)       print "0\t0"
      else if (big && r < 0.32) print int(rand() * 100000) + 200 "\t" int(rand() * 50000)
      else if (big && r < 0.6)  print int(rand() * 120) "\t" int(rand() * 120)
      else               print int(rand() * 6) "\t" int(rand() * 6)
    } }'
}
mk 1 1 | "$YAME" pack -f 3 - > a.cg
mk 2 1 | "$YAME" pack -f 3 - > b.cg
mk 3 0 | "$YAME" pack -f 3 - > small.cg
awk -v N=$N 'BEGIN { for (i = 0; i < N; i++) print (i % 3 == 0) }' |
  "$YAME" pack -f b - > bits.cg
cat a.cg bits.cg b.cg small.cg > s.cg
printf 'a\nbits\nb\nsmall\n' > names.txt
"$YAME" index -s names.txt s.cg

## every line carries the third column; the format-0 record holds "."
[ "$(awk -F'\t' 'NF != 3' s.cg.idx | wc -l)" -eq 0 ] ||
  { echo "index: a line without the third column"; cat s.cg.idx | cut -c1-80; exit 1; }
[ "$(awk -F'\t' '$1 == "bits" {print $3}' s.cg.idx)" = "." ] ||
  { echo "index: the format-0 record should carry '.'"; exit 1; }
for n in a b small; do
  awk -F'\t' -v n=$n '$1 == n && $3 !~ /^r1\./ {bad = 1} END {exit bad}' s.cg.idx ||
    { echo "index: $n has no row index"; exit 1; }
done
## only digits, lower-case letters, '.' and '_': nothing a separator-guessing
## reader (data.table::fread tries , ; : | and space) can take for a column
[ "$(cut -f3 s.cg.idx | tr -d '0-9a-z._\n' | wc -c)" -eq 0 ] ||
  { echo "index: the third column holds a character outside [0-9a-z._]"; exit 1; }
## and the first two columns are what they always were
cut -f1,2 s.cg.idx > two.idx
cp s.cg.idx three.idx

## the full-decode answer: the same file with a two-column index
full() { cp two.idx s.cg.idx; "$YAME" rowsub "$@" s.cg; cp three.idx s.cg.idx; }

## ---- 1. -l: rows in any order, repeated, first and last ---------------------
awk -v N=$N 'BEGIN { srand(7); for (i = 0; i < 300; i++) print int(rand() * N) + 1
  print 1; print N; print 5; print 5; print 65536; print 65537 }' > rows.txt
full -l rows.txt > want.cg
"$YAME" rowsub -l rows.txt s.cg > got.cg 2> err.txt
cmp want.cg got.cg || { echo "-l: row index and full decode differ"; exit 1; }
[ ! -s err.txt ] || { echo "-l: unexpected stderr"; cat err.txt; exit 1; }
## and the values are the ones packed: row 5 of a.cg, via unpack
"$YAME" unpack -f -1 got.cg 2>/dev/null | sed -n '304p' > r5.txt
mk 1 1 | sed -n '5p' > r5.want
diff r5.want r5.txt || { echo "-l: row 5 of a is not the packed value"; exit 1; }

## ---- 2. -B and -I: blocks that span BGZF blocks -----------------------------
for sel in "-B 0_1" "-B 1000_150000" "-B 199990_200000" "-I 3_40000" "-B 199999"; do
  full $sel > want.cg
  "$YAME" rowsub $sel s.cg > got.cg
  cmp want.cg got.cg || { echo "$sel: row index and full decode differ"; exit 1; }
done

## ---- 3. -S: chosen records, in the order asked -------------------------------
printf 'small\na\nbits\n' > pick.txt
full -S pick.txt -l rows.txt > want.cg
"$YAME" rowsub -S pick.txt -l rows.txt s.cg > got.cg
cmp want.cg got.cg || { echo "-S: row index and full decode differ"; exit 1; }
## the same as subset | rowsub, the pipe form, which cannot seek
"$YAME" subset -l pick.txt s.cg | "$YAME" rowsub -l rows.txt - > pipe.cg
cmp want.cg pipe.cg || { echo "-S and subset | rowsub differ"; exit 1; }
## -S needs a file and names the index knows
if "$YAME" rowsub -S pick.txt -l rows.txt - < s.cg > /dev/null 2> err.txt; then
  echo "-S on a stream was accepted"; exit 1
fi
printf 'nope\n' > bad.txt
if "$YAME" rowsub -S bad.txt -l rows.txt s.cg > /dev/null 2> err.txt; then
  echo "-S with an unknown name was accepted"; exit 1
fi
grep -q nope err.txt || { echo "-S unknown name: error does not name it"; exit 1; }

## ---- 4. the index is actually used ------------------------------------------
## Damage a BGZF block in the middle of record a. A full decode must fail on
## it; the row index, asked only for rows in a's first block, never reads it.
cp s.cg dmg.cg; cp three.idx dmg.cg.idx
python3 - dmg.cg <<'PY'
import struct, sys
f = sys.argv[1]; b = bytearray(open(f, "rb").read())
offs, o = [], 0
while o < len(b):
    offs.append(o); o += struct.unpack("<H", b[o+16:o+18])[0] + 1
mid = offs[3]                      # a's fourth block: past every row asked for
for i in range(mid + 40, mid + 80): b[i] ^= 0xFF
open(f, "wb").write(b)
PY
printf '1\n2\n10\n' > early.txt
printf 'a\n' > a.txt
"$YAME" rowsub -S a.txt -l early.txt dmg.cg > got.cg ||
  { echo "row index read a block it did not need"; exit 1; }
full -S a.txt -l early.txt > want.cg
cmp want.cg got.cg || { echo "damaged file: rows from the row index differ"; exit 1; }
cut -f1,2 three.idx > dmg.cg.idx
if "$YAME" rowsub -S a.txt -l early.txt dmg.cg > /dev/null 2>&1; then
  echo "fixture: the damaged block did not break a full decode"; exit 1
fi

## ---- 5. no index, a stale one, a garbled one: read in full ------------------
## two columns, the form every old writer produces: same bytes, no message
cp two.idx s.cg.idx
"$YAME" rowsub -l rows.txt s.cg > got.cg 2> err.txt
full -l rows.txt > want.cg
cmp want.cg got.cg && [ ! -s err.txt ] || { echo "two-column index: wrong"; exit 1; }
## a stale one: the file rewritten behind an index that described it
cp three.idx stale.cg.idx
cat b.cg bits.cg a.cg small.cg > stale.cg       # same names, other bytes
cp stale.cg.idx stale3.idx
"$YAME" rowsub -l rows.txt stale.cg > got.cg 2> err.txt
cut -f1,2 stale3.idx > stale.cg.idx
"$YAME" rowsub -l rows.txt stale.cg > want.cg
cmp want.cg got.cg || { echo "stale index: output differs from a full decode"; exit 1; }
cp stale3.idx stale.cg.idx
"$YAME" rowsub -l rows.txt stale.cg 2>&1 > /dev/null | grep -q 'no longer match' ||
  { echo "stale index: no message"; exit 1; }
## garbage in the column is ignored, not believed
awk -F'\t' 'BEGIN {OFS = "\t"} {print $1, $2, "r1.zz.q_"}' three.idx > s.cg.idx
"$YAME" rowsub -l rows.txt s.cg > got.cg
full -l rows.txt > want.cg
cmp want.cg got.cg || { echo "garbled index: output differs"; exit 1; }
cp three.idx s.cg.idx

## ---- 6. writers carry it ----------------------------------------------------
## subset -o copies records byte for byte, so it moves their row index too
"$YAME" subset -o sub.cg s.cg small a > /dev/null
[ "$(awk -F'\t' '$3 ~ /^r1\./' sub.cg.idx | wc -l)" -eq 2 ] ||
  { echo "subset -o dropped the row index"; cut -c1-60 sub.cg.idx; exit 1; }
"$YAME" rowsub -l rows.txt sub.cg > got.cg
cut -f1,2 sub.cg.idx > sub2.idx; cp sub.cg.idx sub3.idx; cp sub2.idx sub.cg.idx
"$YAME" rowsub -l rows.txt sub.cg > want.cg
cp sub3.idx sub.cg.idx
cmp want.cg got.cg || { echo "subset -o: moved row index gives other rows"; exit 1; }
## index -1 keeps the records already there and adds the new one
cp a.cg ap.cg; "$YAME" index -1 one ap.cg
cat small.cg >> ap.cg; "$YAME" index -1 two ap.cg
[ "$(awk -F'\t' '$3 ~ /^r1\./' ap.cg.idx | wc -l)" -eq 2 ] ||
  { echo "index -1 lost a row index"; cut -c1-60 ap.cg.idx; exit 1; }
"$YAME" rowsub -l rows.txt ap.cg > got.cg
cut -f1,2 ap.cg.idx > ap2.idx; cp ap2.idx ap.cg.idx
"$YAME" rowsub -l rows.txt ap.cg > want.cg
cmp want.cg got.cg || { echo "index -1: row index gives other rows"; exit 1; }

## ---- 7. readers of the two-column form see no change -----------------------
cp three.idx s.cg.idx
"$YAME" info -1 s.cg > i3.txt
cp two.idx s.cg.idx
"$YAME" info -1 s.cg > i2.txt
diff i2.txt i3.txt || { echo "info -1 reads the third column"; exit 1; }
cp three.idx s.cg.idx
"$YAME" subset s.cg b small > m3.cg
cp two.idx s.cg.idx
"$YAME" subset s.cg b small > m2.cg
cmp m2.cg m3.cg || { echo "subset by name reads the third column"; exit 1; }

echo "ok: rowsub through the row index (byte-identical to a full decode)"
