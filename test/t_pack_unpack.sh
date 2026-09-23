#!/bin/bash
## Layer 1: the codec round-trip, one case per format.
##
##   pack -f X in.txt | unpack  ==  in.txt
##
## This is the property most of the format bugs broke without any caller
## noticing: fmt4 indexed as float_t (16 bytes under one FP mode, 4 under
## another), the fmt1 run-length overread, the unaligned uint16 stores, the
## FMT6 macro dropping its index. A format grid reaches branches no command
## reaches.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
cd "$d"

## ---- fixtures, one per format, each with the awkward values inline --------
## fmt0 binary
printf '1\n0\n0\n1\n1\n1\n0\n1\n' > f0.txt                    # 8 rows: byte-exact boundary
## fmt1 one character per line -- 'a' is the "first character only" trap
printf '0\n9\na\n2\n2\n1\nA\nC\n' > f1.txt
## fmt2 states, repeated so the key table is exercised
printf 'Promoter\nPromoter\nExon\nIntron\nIntergenic\nIntergenic\nPromoter\nExon\n' > f2.txt
## fmt3 M/U, with an uncovered site and a long zero run (the RLE flush path)
{ printf '3\t1\n0\t0\n5\t5\n'; for i in $(seq 1 40); do printf '0\t0\n'; done
  printf '7\t0\n1\t9\n'; } > f3.txt
## fmt4 betas with NA
printf '0.750\nNA\n0.500\n0.000\n1.000\nNA\n0.125\n' > f4.txt
## fmt6 set + universe: TWO columns, S and U. Every site in-universe here, so
## the text round-trips exactly; the third state (outside the universe) prints
## as NA rather than 0 and is checked on its own below.
printf '1\t1\n0\t1\n1\t1\n0\t1\n0\t1\n1\t1\n1\t1\n0\t1\n' > f6.txt
## fmt7 reference coordinates. pack reads chrom and START only; the default
## print mode renders each site as the CpG DINUCLEOTIDE it stands for, so the
## end is begin + 2. The fixture is written in that convention, which is what
## a .cr round-trips through.
printf 'chr1\t100\t102\nchr1\t200\t202\nchr2\t50\t52\nchr2\t70\t72\n' > f7.txt

## ---- round-trip each ------------------------------------------------------
## Diagnostics go to stderr and the format code is appended to a file, so a
## failure inside the helper is visible rather than captured by a redirect.
: > fmts.got
rt() {                          # <letter> <fixture> [unpack args...]
  local letter=$1 src=$2; shift 2
  "$YAME" pack -f "$letter" "$src" > rt.cg 2>pack.err ||
    { echo "pack -f $letter failed" >&2; cat pack.err >&2; exit 1; }
  "$YAME" unpack "$@" rt.cg > rt.txt 2>/dev/null
  diff "$src" rt.txt > rt.diff ||
    { echo "round-trip -f $letter changed the data" >&2; head -6 rt.diff >&2; exit 1; }
  ## and the format the header claims is the one the table documents
  printf '%s %s\n' "$letter" \
    "$("$YAME" info rt.cg 2>/dev/null | tail -1 | cut -f5)" >> fmts.got
}
rt b f0.txt
rt c f1.txt
rt s f2.txt
rt m f3.txt -f -1            # M<TAB>U, the two-column form f3.txt is written in
rt n f4.txt
rt d f6.txt -f -1            # value<TAB>universe, the two-column form above
rt r f7.txt
printf 'b 0\nc 1\ns 2\nm 3\nn 4\nd 6\nr 7\n' | sort > fmts.want
sort fmts.got > fmts.sorted
diff fmts.want fmts.sorted ||
  { echo "a pack letter produced the wrong format code"; exit 1; }

## ---- fmt6's third state ---------------------------------------------------
## A site outside the universe is neither set nor unset. It prints as '2' in
## the default mode and as NA<TAB>0 with -f -1; 0<TAB>0 would be a lie, since
## "not set" and "not asked" are different answers.
printf '1\t1\n0\t1\n0\t0\n1\t1\n' > f6u.txt
"$YAME" pack -f d f6u.txt > f6u.cg
[ "$("$YAME" unpack f6u.cg 2>/dev/null | paste -sd, -)" = "1,0,2,1" ] ||
  { echo "fmt6 default mode: out-of-universe is not '2'"; exit 1; }
[ "$("$YAME" unpack -f -1 f6u.cg 2>/dev/null | sed -n 3p)" = "$(printf 'NA\t0')" ] ||
  { echo "fmt6 -f -1: out-of-universe is not NA<TAB>0"; exit 1; }

## ---- fmt7 coordinate print modes ------------------------------------------
## -r 0 (default) chrm/beg0/end1 spanning the dinucleotide; -r 1 the allc form;
## -r 2 the single-string chrm_pos1 form. All three describe the same site, so
## the begin column must agree and only the rendering may differ.
"$YAME" pack -f r f7.txt > f7.cg
[ "$("$YAME" unpack -r 1 f7.cg 2>/dev/null | head -1)" = "$(printf 'chr1\t100\t101')" ] ||
  { echo "fmt7 -r 1 is not the allc rendering"; "$YAME" unpack -r 1 f7.cg | head -1; exit 1; }
[ "$("$YAME" unpack -r 2 f7.cg 2>/dev/null | head -1)" = "chr1_101" ] ||
  { echo "fmt7 -r 2 is not chrm_pos1"; "$YAME" unpack -r 2 f7.cg | head -1; exit 1; }
## every mode reports the same number of sites
for r in 0 1 2; do
  [ "$("$YAME" unpack -r $r f7.cg 2>/dev/null | wc -l)" -eq 4 ] ||
    { echo "fmt7 -r $r lost sites"; exit 1; }
done

## ---- fmt3 print modes agree with each other -------------------------------
## -f -1 gives M and U; -f 1 gives beta; they must describe the same sites.
"$YAME" pack -f m f3.txt > f3.cg
"$YAME" unpack -f -1 f3.cg 2>/dev/null |
  awk -F'\t' '{ print ($1 + $2 < 1) ? "NA" : sprintf("%1.3f", $1 / ($1 + $2)) }' > beta.want
"$YAME" unpack -f 1 f3.cg 2>/dev/null > beta.got
diff beta.want beta.got || { echo "-f 1 disagrees with -f -1"; exit 1; }

## ---- pack -f c refuses a multi-character line ------------------------------
## It used to keep the first character silently, so 255 became '2'.
printf '5\n10\n255\n' > bad1.txt
if "$YAME" pack -f c bad1.txt > /dev/null 2>pack.err; then
  echo "pack -f c accepted a multi-character line"; exit 1
fi
[ -s pack.err ] || { echo "pack -f c failed silently"; exit 1; }

## ---- chunked unpack equals unchunked, including at an exact multiple ------
## `unpack -c` ran one iteration too many when the row count divided by the
## chunk size: every row printed, then "Slicing negative span" and rc 1.
awk 'BEGIN { for (i = 0; i < 40; i++) print (i % 3 == 0) ? 1 : 0 }' > f0_40.txt
"$YAME" pack -f b f0_40.txt > f0_40.cg
"$YAME" unpack f0_40.cg 2>/dev/null > plain.txt
for s in 1 7 10 13 40 100; do
  "$YAME" unpack -c -s $s f0_40.cg > chunk.txt 2>/dev/null ||
    { echo "unpack -c -s $s exited non-zero"; exit 1; }
  diff plain.txt chunk.txt >/dev/null ||
    { echo "unpack -c -s $s differs from unchunked"; exit 1; }
done

## ---- a multi-record store keeps each record's own format ------------------
## read_cdata2() walks a whole file through ONE cdata_t; a stale unit or key
## table from the previous record would show up here and nowhere else.
"$YAME" pack -f n f4.txt > f4.cg
"$YAME" pack -f s f2.txt > f2.cg
cat f3.cg f4.cg f2.cg > mixed.cg
[ "$("$YAME" info mixed.cg 2>/dev/null | tail -n +2 | cut -f5 | paste -sd, -)" = "3,4,2" ] ||
  { echo "mixed-format store: formats not reported per record"; "$YAME" info mixed.cg; exit 1; }
printf 'a\nb\nc\n' > mixed.names
"$YAME" index -s mixed.names mixed.cg
for pair in "a f3.txt -f -1" "b f4.txt" "c f2.txt"; do
  set -- $pair
  name=$1; want=$2; shift 2
  "$YAME" subset mixed.cg "$name" > one.cg 2>/dev/null
  "$YAME" unpack "$@" one.cg > one.txt 2>/dev/null
  diff "$want" one.txt >/dev/null ||
    { echo "record $name did not survive a mixed-format store"; exit 1; }
done

## ---- an output that cannot be opened is an error, not a warning ------------
## pack used to print "Error opening file for writing" and exit 0, so a script
## under set -e learned of a missing output directory three steps later.
if "$YAME" pack -f m f3.txt nodir/out.cg 2>perr.txt; then
  echo "pack exited 0 with its output unopenable"; exit 1
fi
grep -q "Error opening file for writing: nodir/out.cg" perr.txt ||
  { echo "pack did not say which file it could not open"; cat perr.txt; exit 1; }
[ ! -e nodir ] || { echo "pack created the missing directory"; exit 1; }
