#!/bin/bash
## Unit widths.
##
## A format-3 record does not store M and U at a fixed size: the packer picks
## a width from the largest count it sees, so the same file shape has several
## on-disk encodings and each has its own read and write path. A suite whose
## fixtures all hold single-digit counts exercises exactly one of them. Same
## for a state track, whose key index widens once there are more terms than a
## byte can number.
##
## Deep coverage of WGBS counts matters because that is what real data holds:
## a pooled methylome reaches six and seven figures per site.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
cd "$d"

n=24
## one fixture per magnitude: single digits, past a byte, past two bytes,
## past four, and a mixture that forces the widest of them for the record
mags="9 300 70000 5000000"
units=""
for mag in $mags; do
  awk -v n=$n -v m=$mag 'BEGIN {
    for (i = 0; i < n; i++) {
      if (i % 8 == 5) { print 0 "\t" 0; continue }     # an uncovered site
      a = int(m * (i + 1) / n); b = int(m * (n - i) / n)
      print a "\t" b } }' > w$mag.txt
  "$YAME" pack -f m w$mag.txt > w$mag.cg
  u=$("$YAME" info w$mag.cg 2>/dev/null | tail -1 | cut -f6)
  units="$units $u"

  ## the round trip must hold at every width
  "$YAME" unpack -f -1 w$mag.cg 2>/dev/null > w$mag.back
  diff w$mag.txt w$mag.back || { echo "M/U up to $mag did not round-trip (unit $u)"; exit 1; }

  ## and beta must agree with the counts, which is where a wrong width shows
  "$YAME" unpack -f 1 w$mag.cg 2>/dev/null > w$mag.beta
  awk -F'\t' '{ c = $1 + $2; print (c < 1) ? "NA" : sprintf("%1.3f", $1 / c) }' w$mag.txt > w$mag.bwant
  diff w$mag.bwant w$mag.beta || { echo "beta disagrees with the counts at magnitude $mag"; exit 1; }
done

## the widths must actually differ, or this file is testing one path four times
[ "$(echo $units | tr ' ' '\n' | sort -u | wc -l)" -ge 3 ] ||
  { echo "the four magnitudes gave only these units:$units"; exit 1; }

## ---- every command over every width ---------------------------------------
for mag in $mags; do
  "$YAME" rowsub -B 4_12 w$mag.cg 2>/dev/null | "$YAME" unpack -f -1 - 2>/dev/null > s.txt
  sed -n '5,12p' w$mag.txt > s.want
  diff s.want s.txt || { echo "rowsub at magnitude $mag returned the wrong rows"; exit 1; }

  "$YAME" binarize w$mag.cg >/dev/null 2>&1 || { echo "binarize failed at $mag"; exit 1; }
  "$YAME" summary w$mag.cg >/dev/null 2>&1 || { echo "summary failed at $mag"; exit 1; }
  "$YAME" hprint -c -g w$mag.cg >/dev/null 2>&1 || { echo "hprint failed at $mag"; exit 1; }
  "$YAME" dsample -s 1 -N 8 w$mag.cg >/dev/null 2>&1 || { echo "dsample failed at $mag"; exit 1; }

  ## musum across two records of the same width, checked against awk
  cat w$mag.cg w$mag.cg > pair$mag.cg
  "$YAME" rowop -o musum pair$mag.cg 2>/dev/null | "$YAME" unpack -f -1 - 2>/dev/null > ms.txt
  awk -F'\t' '{ print $1 * 2 "\t" $2 * 2 }' w$mag.txt > ms.want
  diff ms.want ms.txt || { echo "musum at magnitude $mag is not twice the input"; exit 1; }

  ## and a summary of one against another of a DIFFERENT width
  "$YAME" pairwise -S w9.cg w$mag.cg >/dev/null 2>&1 ||
    { echo "pairwise across widths failed at $mag"; exit 1; }
done

## a store mixing the widths, read through one cdata_t
cat w9.cg w300.cg w70000.cg w5000000.cg > mixed.cg
printf 'a\nb\nc\ne\n' > nm.txt
"$YAME" index -s nm.txt mixed.cg
for pair in "a w9" "b w300" "c w70000" "e w5000000"; do
  set -- $pair
  "$YAME" subset mixed.cg "$1" 2>/dev/null | "$YAME" unpack -f -1 - 2>/dev/null > got.txt
  diff "$2.txt" got.txt || { echo "record $1 changed inside a mixed-width store"; exit 1; }
done
"$YAME" unpack -a -f 1 mixed.cg 2>/dev/null | awk -F'\t' 'NF != 4 { print "row " NR " has " NF " columns"; exit 1 }'

## ---- unpack -u: the inflated-width override -------------------------------
## Valid widths are 1, 2, 4, 6 and 8; anything else is refused rather than
## used to stride through the data at the wrong pitch.
for u in 1 2 4 6 8; do
  "$YAME" unpack -u $u w9.cg >/dev/null 2>&1 || true    # may or may not suit the record
done
if "$YAME" unpack -u 3 w9.cg >/dev/null 2>&1; then
  echo "unpack -u 3 was accepted; 3 is not an allowed unit"; exit 1
fi
if "$YAME" unpack -u 0 w9.cg >/dev/null 2>&1; then :; fi   # 0 means auto

## ---- format 4: NA runs, the endpoints, and the commands over them ---------
## A beta record stores NA as a negative float and run-compresses stretches of
## it, so a fixture without long NA runs never reaches the run encoder. Real
## arrays are mostly NA outside their probe set, so this is the common shape.
{ for i in $(seq 1 130); do echo NA; done
  printf '0.000
1.000
0.500
'
  for i in $(seq 1 70); do echo NA; done
  printf '0.001
0.999
'; } > na.txt
"$YAME" pack -f n na.txt > na.cg
"$YAME" unpack na.cg 2>/dev/null > na.back
diff na.txt na.back || { echo "a beta record with long NA runs did not round-trip"; exit 1; }
[ "$("$YAME" info na.cg 2>/dev/null | tail -1 | cut -f4)" -eq 205 ] ||
  { echo "the NA fixture lost rows"; exit 1; }

## all NA, and a single value: the degenerate ends of the run encoder
for i in $(seq 1 40); do echo NA; done > allna.txt
"$YAME" pack -f n allna.txt > allna.cg
[ "$("$YAME" unpack allna.cg 2>/dev/null | sort -u)" = "NA" ] ||
  { echo "an all-NA beta record did not decode as NA"; exit 1; }
printf '0.250
' > one4.txt
"$YAME" pack -f n one4.txt > one4.cg
[ "$("$YAME" unpack one4.cg 2>/dev/null)" = "0.250" ] || { echo "a single beta did not round-trip"; exit 1; }

## the commands, over a record that is mostly NA
"$YAME" rowsub -B 128_136 na.cg 2>/dev/null | "$YAME" unpack - 2>/dev/null > na.sub
sed -n '129,136p' na.txt > na.subwant
diff na.subwant na.sub || { echo "slicing across an NA run returned the wrong rows"; exit 1; }
for cmd in "summary na.cg" "hprint -c na.cg" "hprint -c -g na.cg"; do
  "$YAME" $cmd >/dev/null 2>&1 || { echo "yame $cmd failed on a mostly-NA beta record"; exit 1; }
done
## a mask over it, and a pairing against itself: agreement must be perfect
awk 'BEGIN { for (i = 0; i < 205; i++) print (i % 3 == 0) ? 1 : 0 }' |
  "$YAME" pack -f b - > na.mask
"$YAME" summary -m na.mask na.cg >/dev/null 2>&1 || { echo "summary -m failed on the NA record"; exit 1; }
"$YAME" pairwise -S na.cg na.cg > self.txt 2>/dev/null
awk -F'	' 'NR == 2 && !($4 == "0.000000" && $7 == "1.000000") {
  print "a beta record scored against itself is not a perfect match"; print; exit 1 }' self.txt
## only the non-NA sites are scored
nn=$(grep -vc NA na.txt)
[ "$(awk -F'	' 'NR==2{print $3}' self.txt)" = "$nn" ] ||
  { echo "pairwise scored $(awk -F'	' 'NR==2{print $3}' self.txt) sites, want the $nn non-NA ones"; exit 1; }

## ---- a state track wider than one byte of keys ----------------------------
## The key index is a byte while there are fewer than 256 terms and widens
## past that, so a track with 300 distinct states reads through the other path.
awk 'BEGIN { for (i = 0; i < 600; i++) printf "state%03d\n", i % 300 }' > many.txt
"$YAME" pack -f s many.txt > many.cg
"$YAME" unpack many.cg 2>/dev/null > many.back
diff many.txt many.back || { echo "a 300-key state track did not round-trip"; exit 1; }
nkeys=$("$YAME" info many.cg 2>/dev/null | tail -1 | cut -f7 | sed 's/N=\([0-9]*\).*/\1/')
[ "$nkeys" = "300" ] || { echo "a 300-key track reports $nkeys keys"; exit 1; }
## and it still slices correctly, by every selector
printf '1\n300\n600\n' > idx.txt
"$YAME" rowsub -l idx.txt many.cg 2>/dev/null | "$YAME" unpack - 2>/dev/null > many.sub
{ sed -n '1p;300p;600p' many.txt; } > many.subwant
diff many.subwant many.sub || { echo "slicing a wide-key state track returned the wrong rows"; exit 1; }
"$YAME" summary -m many.cg many.cg >/dev/null 2>&1 || { echo "summary over a wide-key mask failed"; exit 1; }
