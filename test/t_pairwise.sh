#!/bin/bash
## yame pairwise: the set output, and the -S summary added for scoring a
## prediction against truth (see the 20260911 entry in the YAME org).
##
## Every expected value is recomputed here in awk from the same text the
## fixtures are packed from, in the same row order and with the same
## formulas, so the check is against an independent implementation rather
## than a pasted constant. Nothing is downloaded.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
cd "$d"

## ---- fixtures: 20 rows; "0 0" is an uncovered site, NA a missing beta ----
cat > a.txt <<'T'
3	1
0	0
5	5
2	8
7	0
0	0
1	3
4	4
9	1
0	6
2	2
0	0
6	1
1	9
3	3
5	0
0	0
8	2
2	5
4	1
T
cat > b.txt <<'T'
1	3
2	2
5	5
0	0
6	1
4	0
1	3
0	0
8	2
1	5
3	1
2	2
0	0
2	8
3	3
0	5
1	1
7	3
0	0
1	4
T
cat > c.txt <<'T'
0.750
NA
0.500
0.200
1.000
0.100
NA
0.500
0.900
0.000
0.500
0.300
NA
0.100
0.500
1.000
0.400
0.800
NA
0.800
T
"$YAME" pack -f m a.txt > a.cg
"$YAME" pack -f m b.txt > b.cg
"$YAME" pack -f n c.txt > c.cg
cat a.cg b.cg c.cg > abc.cg
printf 'sA\nsB\nsC\n' > names.txt
"$YAME" index -s names.txt abc.cg

## beta from an M/U line, or "NA" when uncovered (M+U < 1)
beta_mu='function beta(m, u) { return (m + u < 1) ? "NA" : m / (m + u) }'

## ---- 1. set output (-H 3, any difference): 1 differs, 0 same, 2 outside ----
paste a.txt b.txt |
  awk -F'\t' "$beta_mu"'{
    b1 = beta($1, $2); b2 = beta($3, $4)
    if (b1 == "NA" || b2 == "NA") print 2
    else print (b1 != b2) ? 1 : 0 }' > set.expected
"$YAME" pairwise -H 3 a.cg b.cg | "$YAME" unpack - > set.got
diff set.expected set.got || { echo "set output (-H 3) differs"; exit 1; }

## ---- 2. set output (-H 1, hyper in sample 1, delta 0.2) ----
paste a.txt b.txt |
  awk -F'\t' "$beta_mu"'{
    b1 = beta($1, $2); b2 = beta($3, $4)
    if (b1 == "NA" || b2 == "NA") print 2
    else print (b1 > b2 && b1 - b2 > 0.2) ? 1 : 0 }' > set1.expected
"$YAME" pairwise -H 1 -d 0.2 a.cg b.cg | "$YAME" unpack - > set1.got
diff set1.expected set1.got || { echo "set output (-H 1 -d 0.2) differs"; exit 1; }

## summary statistics over rows where both sides have a value; raw-moment
## Pearson, accumulated in row order exactly as the C code does
stats='function flush(l1, l2,   den, r) {
    den = (n * sxx - sx * sx) * (n * syy - sy * sy)
    r = (n >= 2 && den > 0) ? (n * sxy - sx * sy) / sqrt(den) : "nan"
    printf "%s\t%s\t%d\t%.6f\t%.6f\t%.6f\t%.6f\n", l1, l2, n,
           sae / n, sqrt(sse / n), agree / n, r }
  function add(x, y,   dd) {
    dd = x - y; n++; sae += (dd < 0 ? -dd : dd); sse += dd * dd
    agree += ((x >= t) == (y >= t)); sx += x; sy += y
    sxx += x * x; syy += y * y; sxy += x * y }'

## ---- 3. -S, format 3 vs format 3 ----
paste a.txt b.txt |
  awk -F'\t' -v t=0.5 "$beta_mu $stats"'
    { b1 = beta($1, $2); b2 = beta($3, $4)
      if (b1 != "NA" && b2 != "NA") add(b1, b2) }
    END { flush("1", "1") }' > s33.expected
"$YAME" pairwise -S a.cg b.cg | tail -n +2 > s33.got
diff s33.expected s33.got || { echo "-S fmt3 vs fmt3 differs"; exit 1; }

## ---- 4. -S, format 3 vs format 4, and -t moves only acc ----
for t in 0.5 0.3; do
  paste a.txt c.txt |
    awk -F'\t' -v t=$t "$beta_mu $stats"'
      { b1 = beta($1, $2); b2 = $3
        if (b1 != "NA" && b2 != "NA") add(b1, b2) }
      END { flush("1", "1") }' > s34.expected
  "$YAME" pairwise -S -t $t a.cg c.cg | tail -n +2 > s34.got
  diff s34.expected s34.got || { echo "-S fmt3 vs fmt4 (-t $t) differs"; exit 1; }
done

## ---- 5. names through the index, and the header line ----
printf 'sample1\tsample2\tn\tmae\trmse\tacc\tpearson\n' > named.expected
## A literal tab, not \t: whether sed understands the escape is a dialect
## question, and BSD sed's answer is "a letter t".
tab=$(printf '\t')
sed "s/^1${tab}1${tab}/sA${tab}sC${tab}/" s34.expected >> named.expected  # last s34 run was t=0.3
"$YAME" pairwise -S -t 0.3 -1 sA -2 sC abc.cg abc.cg > named.got
diff named.expected named.got || { echo "-1/-2 by name differs"; exit 1; }

## ---- 6. broadcast: one named record against every record of file 2 ----
"$YAME" pairwise -S -1 sA abc.cg abc.cg > bc.got
[ "$(wc -l < bc.got)" -eq 4 ] || { echo "broadcast: expected 3 lines + header"; cat bc.got; exit 1; }
[ "$(cut -f2 bc.got | tail -n +2 | paste -sd, -)" == "sA,sB,sC" ] ||
  { echo "broadcast: labels not in file order"; cat bc.got; exit 1; }
## sA against itself: perfect agreement
awk -F'\t' 'NR == 2 && !($4 == "0.000000" && $6 == "1.000000" && $7 == "1.000000") {
  print "broadcast: sA vs sA is not a perfect match"; print; exit 1 }' bc.got

## ---- 7. refusals: each exits 1, says why, and leaves stdout empty ----
refuse() {                      # <expected message fragment> <args...>
  local want=$1; shift
  local out err rc
  out=$("$YAME" pairwise "$@" 2>err.txt) && rc=0 || rc=$?
  err=$(cat err.txt)
  [ "$rc" -eq 1 ] || { echo "expected exit 1 for: $*"; echo "$err"; exit 1; }
  [ -z "$out" ] || { echo "expected empty stdout for: $*"; echo "$out"; exit 1; }
  case "$err" in *"$want"*) ;; *) echo "wrong message for: $*"; echo "$err"; exit 1;; esac
}
printf '1\t1\n2\t2\n' | "$YAME" pack -f m - > short.cg
refuse "more than one record"   -S abc.cg abc.cg
refuse "differ in length"       -S a.cg short.cg
refuse "not in the index"       -S -1 nope abc.cg b.cg
refuse "need an index"          -S -1 sA a.cg b.cg
"$YAME" binarize a.cg > a6.cg
refuse "need format 3"          -S a6.cg b.cg
