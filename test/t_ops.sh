#!/bin/bash
## The remaining rowop operations, and the option surfaces of subset, summary
## and dsample that the first pass did not reach.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
cd "$d"

rows=24
for s in 1 2 3; do
  awk -v s=$s -v n=$rows 'BEGIN {
    for (i = 0; i < n; i++) {
      if (i == 5) { print 0 "\t" 0; continue }       # uncovered everywhere
      m = ((i * s) % 11); print m "\t" (10 - m) } }' > s$s.txt
  "$YAME" pack -f m s$s.txt > s$s.cg
done
cat s1.cg s2.cg s3.cg > three.cg
printf 'a\nb\nc\n' > names.txt
"$YAME" index -s names.txt three.cg

## ---- 1. rowop binstring: one binary string per row ------------------------
## Not a plain threshold: an ambiguous site is FILLED from the majority of the
## confident calls in its row, and when that majority is not convincing the
## whole line is emitted as the all-'2' sentinel. So the exact-agreement check
## uses a fixture with no ambiguity at all, and the fill rules are checked on
## their own below.
for s in 1 2 3; do
  awk -v s=$s -v n=$rows 'BEGIN {
    for (i = 0; i < n; i++) {
      if (i == 5) { print 0 "\t" 0; continue }          # uncovered everywhere
      hi = ((i + s) % 2)                                 # 0.9 or 0.1, never near 0.5
      print (hi ? 9 : 1) "\t" (hi ? 1 : 9) } }' > c$s.txt
  "$YAME" pack -f m c$s.txt > c$s.cg
done
cat c1.cg c2.cg c3.cg > clear.cg
"$YAME" rowop -o binstring -b 0.5 -c 1 clear.cg > bs.txt 2>/dev/null
[ "$(wc -l < bs.txt)" -eq $rows ] || { echo "binstring printed $(wc -l < bs.txt) rows, want $rows"; exit 1; }
awk -v n=3 'length($0) != n { print "row " NR " is " length($0) " chars, want " n; exit 1 }
            /[^012]/ { print "row " NR " has a glyph outside 012: " $0; exit 1 }' bs.txt
paste c1.txt c2.txt c3.txt |
  awk -F'\t' '{ out = ""
    for (k = 1; k <= 5; k += 2) { cov = $k + $(k+1)
      out = out ((cov < 1) ? "2" : (($k / cov > 0.5) ? "1" : "0")) }
    print out }' > bs.want
diff bs.want bs.txt || { echo "binstring disagrees on unambiguous betas"; exit 1; }

## A row whose confident calls are evenly split cannot be filled, so the line
## becomes the sentinel rather than a guess.
printf '9\t1\n' > t1.txt; printf '5\t5\n' > t2.txt; printf '1\t9\n' > t3.txt
for i in 1 2 3; do "$YAME" pack -f m t$i.txt > t$i.cg; done
cat t1.cg t2.cg t3.cg > tie.cg
[ "$("$YAME" rowop -o binstring -b 0.5 tie.cg 2>/dev/null)" = "222" ] ||
  { echo "a tied row was not emitted as the sentinel"; "$YAME" rowop -o binstring tie.cg; exit 1; }

## With a clear majority the ambiguous sample is filled rather than sentinelled.
printf '9\t1\n' > u1.txt; printf '5\t5\n' > u2.txt; printf '8\t2\n' > u3.txt
for i in 1 2 3; do "$YAME" pack -f m u$i.txt > u$i.cg; done
cat u1.cg u2.cg u3.cg > maj.cg
[ "$("$YAME" rowop -o binstring -b 0.5 maj.cg 2>/dev/null)" = "111" ] ||
  { echo "a clear majority did not fill the ambiguous call"; "$YAME" rowop -o binstring maj.cg; exit 1; }

## -m: tighten the allowed ambiguous fraction and the sentinel appears
"$YAME" rowop -o binstring -c 99 -m 0.1 three.cg 2>/dev/null > bs2.txt
grep -qx '222' bs2.txt || { echo "-m did not emit an all-2 sentinel line"; head -3 bs2.txt; exit 1; }

## ---- 2. rowop cometh on a real format-3 input -----------------------------
"$YAME" rowop -o cometh -w 2 s1.cg > cm.txt 2>/dev/null
[ -s cm.txt ] || { echo "cometh printed nothing"; exit 1; }
## -v spells each window's packed uint64 out as UU-UM-MU-MM instead, so every
## field after the row index is four dash-separated counts.
"$YAME" rowop -o cometh -w 2 -v s1.cg > cmv.txt 2>/dev/null
awk -F'\t' 'NR == 1 { for (i = 2; i <= NF; i++)
                  if ($i !~ /^[0-9]+-[0-9]+-[0-9]+-[0-9]+$/) {
                    print "field " i " is not UU-UM-MU-MM: " $i; exit 1 }
                if (NF < 2) { print "no window columns"; exit 1 } }' cmv.txt
## and the two forms describe the same number of windows
[ "$(head -1 cm.txt | awk -F'\t' '{print NF}')" -eq "$(head -1 cmv.txt | awk -F'\t' '{print NF}')" ] ||
  { echo "cometh -v has a different column count from the packed form"; exit 1; }

## ---- 3. rowop on the other input formats ----------------------------------
## binasum takes fmt0 and fmt1 as well as fmt3.
awk -v n=$rows 'BEGIN { for (i = 0; i < n; i++) print (i % 2) }' > b.txt
"$YAME" pack -f b b.txt > b1.cg
awk -v n=$rows 'BEGIN { for (i = 0; i < n; i++) print ((i + 1) % 2) }' |
  "$YAME" pack -f b - > b2.cg
cat b1.cg b2.cg > bb.cg
"$YAME" rowop -o binasum bb.cg 2>/dev/null | "$YAME" unpack -f -1 - 2>/dev/null > ba0.txt
[ "$(wc -l < ba0.txt)" -eq $rows ] || { echo "binasum on fmt0 lost rows"; exit 1; }
## every row has exactly one methylated and one unmethylated call
awk -F'\t' '$1 + $2 != 2 { print "row " NR ": " $0; exit 1 }' ba0.txt

## ---- 4. rowop -t threads must agree with the serial answer ----------------
"$YAME" rowop -o stat three.cg > st1.txt 2>/dev/null
"$YAME" rowop -o stat -t 3 three.cg > st4.txt 2>/dev/null
diff st1.txt st4.txt || { echo "rowop -t 3 disagrees with the serial run"; exit 1; }

## ---- 4b. rowop -t without an index: the streaming record queue ------------
## With an index the parallel path seeks to each record; without one it reads
## the stream once and hands records to the workers through a queue. Both must
## give the serial answer, and the queue is only reached by the second.
cp three.cg noidx.cg                       # deliberately no .idx beside it
for op in binasum musum stat; do
  "$YAME" rowop -o $op three.cg > par_idx.txt 2>/dev/null
  "$YAME" rowop -o $op -t 3 noidx.cg > par_str.txt 2>/dev/null
  diff par_idx.txt par_str.txt >/dev/null ||
    { echo "rowop -o $op -t 3 without an index disagrees with the serial answer"; exit 1; }
done
## and through a pipe, where there is no file to seek in at all
"$YAME" rowop -o stat -t 2 - < three.cg > par_pipe.txt 2>/dev/null
"$YAME" rowop -o stat three.cg > ser.txt 2>/dev/null
diff ser.txt par_pipe.txt >/dev/null ||
  { echo "rowop -t over a pipe disagrees with the serial answer"; exit 1; }

## ---- 5. subset -l: names from a list file ---------------------------------
printf 'c\na\n' > want.txt
"$YAME" subset -l want.txt three.cg > sl.cg 2>/dev/null
"$YAME" unpack -f -1 -a sl.cg 2>/dev/null | cut -f1 > sl.first
"$YAME" unpack -f -1 s3.cg 2>/dev/null | cut -f1 > c.first
diff c.first sl.first || { echo "subset -l did not put c first"; exit 1; }
[ "$("$YAME" info sl.cg 2>/dev/null | tail -n +2 | wc -l)" -eq 2 ] ||
  { echo "subset -l returned the wrong record count"; exit 1; }

## ---- 6. summary options ---------------------------------------------------
awk -v n=$rows 'BEGIN { for (i = 0; i < n; i++) print (i < 8) ? 1 : 0 }' | "$YAME" pack -f b - > m1.cg
awk -v n=$rows 'BEGIN { for (i = 0; i < n; i++) print (i % 2 == 0) ? 1 : 0 }' | "$YAME" pack -f b - > m2.cg
cat m1.cg m2.cg > masks.cg
printf 'first8\neven\n' > mnames.txt
"$YAME" index -s mnames.txt masks.cg
"$YAME" summary -m masks.cg three.cg > sm.txt 2>/dev/null
[ "$(tail -n +2 sm.txt | wc -l)" -ge 6 ] ||
  { echo "summary over 3 queries x 2 masks gave $(tail -n +2 sm.txt | wc -l) rows"; cat sm.txt; exit 1; }
grep -q 'first8' sm.txt || { echo "summary does not name the masks"; head -2 sm.txt; exit 1; }
## -H suppresses the header, -M preloads, -F prints full paths
[ "$("$YAME" summary -H -m masks.cg three.cg 2>/dev/null | wc -l)" \
  -eq "$(( $(wc -l < sm.txt) - 1 ))" ] || { echo "-H did not remove exactly the header"; exit 1; }
"$YAME" summary -M -m masks.cg three.cg 2>/dev/null | tail -n +2 > smM.txt
diff <(tail -n +2 sm.txt) smM.txt || { echo "-M changed the answer"; exit 1; }
## -F prints the path as given rather than its basename, so it can only
## differ when the path HAS directories -- pass one that does.
mkdir -p sub && cp three.cg sub/ && cp masks.cg sub/ && cp masks.cg.idx sub/ 2>/dev/null || true
"$YAME" summary -m sub/masks.cg sub/three.cg 2>/dev/null | tail -1 | cut -f1 > base.txt
"$YAME" summary -F -m sub/masks.cg sub/three.cg 2>/dev/null | tail -1 | cut -f1 > full.txt
[ "$(cat base.txt)" = "three.cg" ] || { echo "default QFile is not the basename: $(cat base.txt)"; exit 1; }
[ "$(cat full.txt)" = "sub/three.cg" ] || { echo "-F QFile is not the given path: $(cat full.txt)"; exit 1; }
## a state (format 2) query summarises by term
printf 'Promoter\nExon\nExon\nIntron\n' > st.txt
for i in $(seq 5 $rows); do echo Intergenic; done >> st.txt
"$YAME" pack -f s st.txt > states.cg
"$YAME" summary -m masks.cg states.cg > sms.txt 2>/dev/null
grep -qE 'Promoter|Exon' sms.txt || { echo "summary of a fmt2 query does not name terms"; head -3 sms.txt; exit 1; }

## ---- 7. dsample options ---------------------------------------------------
## -r replicates, -p prefix, -b binarize, and fmt6 input
"$YAME" dsample -s 3 -N 8 -r 3 -p rep three.cg > ds.cg 2>/dev/null
[ "$("$YAME" info ds.cg 2>/dev/null | tail -n +2 | wc -l)" -eq 9 ] ||
  { echo "dsample -r 3 on 3 samples gave $("$YAME" info ds.cg | tail -n +2 | wc -l) records, want 9"; exit 1; }
"$YAME" dsample -s 3 -N 8 -b s1.cg > dsb.cg 2>/dev/null
[ "$("$YAME" info dsb.cg 2>/dev/null | tail -1 | cut -f5)" = "3" ] ||
  { echo "dsample -b changed the format"; exit 1; }
"$YAME" binarize s1.cg > bin.cg 2>/dev/null
"$YAME" dsample -s 3 -N 8 bin.cg > ds6.cg 2>/dev/null
[ "$("$YAME" info ds6.cg 2>/dev/null | tail -1 | cut -f5)" = "6" ] ||
  { echo "dsample on fmt6 did not stay fmt6"; exit 1; }
