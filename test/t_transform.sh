#!/bin/bash
## The transforms that rewrite values in place: binarize, mask, dsample,
## perturb. Two of them take a seed, so "reproducible for a fixed seed" is
## part of the contract; dsample also wrote one byte past its bitset whenever
## the row count divided by 8.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
cd "$d"

## 48 rows (divides by 8 -- the dsample bitset boundary), betas spread across
## the threshold, with uncovered sites mixed in
awk 'BEGIN {
  for (i = 0; i < 48; i++) {
    if (i % 8 == 3) print 0 "\t" 0          # uncovered
    else print (i % 10) "\t" (10 - (i % 10))
  } }' > mu.txt
"$YAME" pack -f m mu.txt > mu.cg

## ---- 1. binarize: beta >= t is set; only covered sites are in-universe ----
"$YAME" binarize -t 0.5 -c 1 mu.cg > bin.cg 2>/dev/null
"$YAME" unpack bin.cg 2>/dev/null > bin.got
awk -F'\t' '{ cov = $1 + $2
              if (cov < 1) print 2                       # outside the universe
              else print ($1 / cov >= 0.5) ? 1 : 0 }' mu.txt > bin.want
diff bin.want bin.got || { echo "binarize -t 0.5 disagrees with the beta rule"; exit 1; }

## ---- 1b. binarize -o on an INDEXED input: every record, and an index -----
## The documented example is `binarize -t 0.5 -c 3 -o calls.cg <indexed>`, and
## until v1.45 it wrote the file SHORT. binarize read its own output back to
## learn each record's BGZF offset, and did that while the output was still
## open with the last record in the buffer. The read hit a truncated record and
## exited, so the buffer never reached disk. The stdout form was fine, and
## every test here used it, which is why nothing caught this.
"$YAME" pack -f m mu.txt > one.cg
cat one.cg one.cg one.cg > three.cg
printf 'sA\nsB\nsC\n' > three.names
"$YAME" index -s three.names three.cg >/dev/null 2>&1
[ -s three.cg.idx ] || { echo "fixture: index -s wrote nothing"; exit 1; }

"$YAME" binarize -t 0.5 -c 1 -o bo.cg three.cg 2>bo.err || {
  echo "binarize -o failed on an indexed input"; cat bo.err; exit 1; }

n=$("$YAME" info bo.cg 2>/dev/null | tail -n +2 | grep -c .)
[ "$n" -eq 3 ] || { echo "binarize -o wrote $n records, want 3"; exit 1; }

## the index must name all three, and each name must resolve
[ "$(grep -c . bo.cg.idx)" -eq 3 ] ||
  { echo "binarize -o wrote $(grep -c . bo.cg.idx) index rows, want 3"; cat bo.cg.idx; exit 1; }
for nm in sA sB sC; do
  "$YAME" subset bo.cg "$nm" > sub.cg 2>/dev/null
  [ -s sub.cg ] || { echo "the index binarize wrote does not resolve $nm"; exit 1; }
done

## -o and the stdout form must produce the same bytes
"$YAME" binarize -t 0.5 -c 1 three.cg > bstdout.cg 2>/dev/null
cmp -s bo.cg bstdout.cg ||
  { echo "binarize -o differs from the stdout form"; ls -la bo.cg bstdout.cg; exit 1; }

## the threshold actually moves the calls
"$YAME" binarize -t 0.9 mu.cg 2>/dev/null | "$YAME" unpack - 2>/dev/null > bin9.got
[ "$(grep -c '^1$' bin9.got)" -lt "$(grep -c '^1$' bin.got)" ] ||
  { echo "binarize: raising -t did not reduce the set"; exit 1; }

## a higher coverage floor can only shrink the universe
"$YAME" binarize -c 20 mu.cg 2>/dev/null | "$YAME" unpack - 2>/dev/null > binc.got
[ "$(grep -c '^2$' binc.got)" -ge "$(grep -c '^2$' bin.got)" ] ||
  { echo "binarize: raising -c did not grow the out-of-universe set"; exit 1; }

## ---- 2. mask: blank the sites the mask covers, -v the complement ----------
## With a format-3 mask, "covered" means M+U > 0. This is how truth is
## restricted to the sites a model was NOT shown.
awk 'BEGIN { for (i = 0; i < 48; i++) print (i % 3 == 0) ? "1\t1" : "0\t0" }' > msk.txt
"$YAME" pack -f m msk.txt > msk.cg
"$YAME" mask mu.cg msk.cg 2>/dev/null | "$YAME" unpack -f -1 - 2>/dev/null > masked.txt
paste msk.txt mu.txt |
  awk -F'\t' '{ if ($1 + $2 > 0) print "0\t0"; else print $3 "\t" $4 }' > masked.want
diff masked.want masked.txt || { echo "mask did not blank exactly the covered sites"; exit 1; }

"$YAME" mask -v mu.cg msk.cg 2>/dev/null | "$YAME" unpack -f -1 - 2>/dev/null > maskedv.txt
paste msk.txt mu.txt |
  awk -F'\t' '{ if ($1 + $2 > 0) print $3 "\t" $4; else print "0\t0" }' > maskedv.want
diff maskedv.want maskedv.txt || { echo "mask -v is not the complement"; exit 1; }

## every site is blanked by exactly one of the two senses
paste masked.txt maskedv.txt |
  awk -F'\t' '$1 + $2 > 0 && $3 + $4 > 0 { bad++ } END { exit bad > 0 }' ||
  { echo "mask and mask -v both kept the same site"; exit 1; }

## ---- 2b. the same sense for a format-6 input, and a format-6 mask --------
## The fmt6 path once kept only the covered sites, the opposite of fmt3 and
## fmt0, so masking a binarized methylome with a blacklist kept the blacklist.
"$YAME" binarize -t 0.5 -c 1 mu.cg > bin.cg 2>/dev/null
"$YAME" mask bin.cg msk.cg 2>/dev/null | "$YAME" unpack -f -1 - 2>/dev/null > masked6.txt
"$YAME" unpack -f -1 bin.cg 2>/dev/null > bin.txt
paste msk.txt bin.txt |
  awk -F'\t' '{ if ($1 + $2 > 0) print "NA\t0"; else print $3 "\t" $4 }' > masked6.want
diff masked6.want masked6.txt || { echo "fmt6 mask did not blank exactly the covered sites"; exit 1; }
## a methylome as the mask: its universe is what it covers, so `mask truth
## query` blanks the sites the query was shown
"$YAME" mask mu.cg bin.cg 2>/dev/null | "$YAME" unpack -f -1 - 2>/dev/null > maskedq.txt
paste bin.txt mu.txt |
  awk -F'\t' '{ if ($1 != "NA") print "0\t0"; else print $3 "\t" $4 }' > maskedq.want
diff maskedq.want maskedq.txt || { echo "a fmt6 mask did not reduce to its universe"; exit 1; }
## a name that is not a file resolves in the store for this row space, and
## the refusal says so rather than "Error opening file"
"$YAME" mask mu.cg Blacklist > /dev/null 2> nm.err && { echo "an unresolvable mask name was accepted"; exit 1; }
grep -q "Blacklist" nm.err || { echo "the refusal does not name the mask"; cat nm.err; exit 1; }
grep -q "Error opening file" nm.err && { echo "a bare name still goes straight to fopen"; cat nm.err; exit 1; }
## and a mask of the wrong format is named, with its format
"$YAME" pack -f i <(awk 'BEGIN{for(i=0;i<48;i++) print i}') > f2.cg 2>/dev/null || true
if [ -s f2.cg ]; then
  "$YAME" mask mu.cg f2.cg > /dev/null 2> wf.err && { echo "a format-2 mask was accepted"; exit 1; }
  grep -q "f2.cg is format" wf.err || { echo "the format refusal does not name the file and format"; cat wf.err; exit 1; }
fi

## ---- 3. dsample: reproducible for a fixed seed, and a no-op at the edge ---
"$YAME" dsample -s 42 -N 20 mu.cg > ds1.cg 2>/dev/null
"$YAME" dsample -s 42 -N 20 mu.cg > ds2.cg 2>/dev/null
cmp ds1.cg ds2.cg || { echo "dsample -s 42 is not reproducible"; exit 1; }
"$YAME" dsample -s 7 -N 20 mu.cg > ds3.cg 2>/dev/null
cmp -s ds1.cg ds3.cg && { echo "dsample ignored the seed"; exit 1; }

## the kept sites are a SUBSET of the covered ones, and there are N of them
"$YAME" unpack -f -1 ds1.cg 2>/dev/null > ds1.txt
[ "$(wc -l < ds1.txt)" -eq 48 ] || { echo "dsample changed the row count"; exit 1; }
[ "$(awk -F'\t' '$1 + $2 > 0' ds1.txt | wc -l)" -eq 20 ] ||
  { echo "dsample -N 20 kept $(awk -F'\t' '$1+$2>0' ds1.txt | wc -l) covered sites"; exit 1; }
paste ds1.txt mu.txt | awk -F'\t' '$1 + $2 > 0 && $3 + $4 == 0 { bad++ } END { exit bad > 0 }' ||
  { echo "dsample kept a site the input did not cover"; exit 1; }

## ---- 4. perturb: reproducible, and -p 0 changes nothing -------------------
"$YAME" binarize mu.cg > b.cg 2>/dev/null
"$YAME" perturb -s 42 -p 0.25 b.cg > p1.cg 2>/dev/null
"$YAME" perturb -s 42 -p 0.25 b.cg > p2.cg 2>/dev/null
cmp p1.cg p2.cg || { echo "perturb -s 42 is not reproducible"; exit 1; }
"$YAME" perturb -s 42 -p 0 b.cg 2>/dev/null | "$YAME" unpack - 2>/dev/null > p0.txt
"$YAME" unpack b.cg 2>/dev/null > b.txt
diff b.txt p0.txt || { echo "perturb -p 0 changed the data"; exit 1; }
## and a flip rate above zero changes something, but not the row count
"$YAME" unpack p1.cg 2>/dev/null > p1.txt
[ "$(wc -l < p1.txt)" -eq "$(wc -l < b.txt)" ] || { echo "perturb changed the row count"; exit 1; }
cmp -s b.txt p1.txt && { echo "perturb -p 0.25 changed nothing"; exit 1; }
exit 0
