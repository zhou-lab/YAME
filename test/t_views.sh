#!/bin/bash
## The format-6 views, the state-mask path, and subset's own option surface.
##
## Format 6 packs two bits per row and three different meanings have been put
## on them: a set within a universe, a methylation call within a covered site,
## and four raw quaternary states. `summary -V` picks which one is being
## asked about, and each has its own counting code -- so a file summarised
## under the wrong view returns confident numbers about a question nobody
## asked. Each view is checked against an independent recount here.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
cd "$d"

n=40
## S and U chosen so all four quaternary states occur
awk -v n=$n 'BEGIN {
  for (i = 0; i < n; i++) {
    r = i % 4
    if (r == 0) print "1\t1"        # in universe, set
    else if (r == 1) print "0\t1"   # in universe, unset
    else if (r == 2) print "0\t0"   # outside the universe
    else print "1\t1"
  } }' > q.txt
"$YAME" pack -f d q.txt > q.cg
awk -v n=$n 'BEGIN { for (i = 0; i < n; i++) print (i < 20) ? 1 : 0 }' > m.txt
"$YAME" pack -f b m.txt > m.cg

## ---- 1. -V set: universe is the background, set is feature membership ----
"$YAME" summary -V set -m m.cg q.cg > vset.txt 2>/dev/null
read -r nu nq nm no <<<"$(tail -1 vset.txt | cut -f5-8 | tr '\t' ' ')"
paste q.txt m.txt | awk -F'\t' '
  { inuni = ($1 + $2 > 0); set = ($1 == 1 && $2 == 1); inmask = ($3 == 1)
    if (inuni) u++
    if (inuni && set) q++
    if (inuni && inmask) mm++
    if (inuni && set && inmask) o++ }
  END { print u, q, mm, o }' > vset.want
read -r wu wq wm wo < vset.want
[ "$nu" = "$wu" ] || { echo "-V set N_univ is $nu, want $wu"; exit 1; }
[ "$nq" = "$wq" ] || { echo "-V set N_query is $nq, want $wq"; exit 1; }
[ "$nm" = "$wm" ] || { echo "-V set N_mask is $nm, want $wm"; exit 1; }
[ "$no" = "$wo" ] || { echo "-V set N_overlap is $no, want $wo"; exit 1; }

## ---- 2. -V meth: universe is coverage, set is methylated ------------------
"$YAME" summary -V meth -m m.cg q.cg > vmeth.txt 2>/dev/null
head -1 vmeth.txt | grep -q 'N_covered' ||
  { echo "-V meth did not switch the column names"; head -1 vmeth.txt; exit 1; }
read -r nc nmeth ncm nmm <<<"$(tail -1 vmeth.txt | cut -f5-8 | tr '\t' ' ')"
[ "$nc" = "$wu" ] || { echo "-V meth N_covered is $nc, want $wu"; exit 1; }
[ "$nmeth" = "$wq" ] || { echo "-V meth N_meth is $nmeth, want $wq"; exit 1; }

## ---- 3. -V 2bit: the four states counted separately -----------------------
"$YAME" summary -V 2bit -m m.cg q.cg > v2.txt 2>/dev/null
[ -s v2.txt ] || { echo "-V 2bit produced nothing"; exit 1; }
[ "$(tail -n +2 v2.txt | wc -l)" -ge 1 ] || { echo "-V 2bit printed no rows"; exit 1; }
## -6 is the deprecated alias and must agree
"$YAME" summary -6 -m m.cg q.cg 2>/dev/null | tail -n +2 > v6.txt
diff <(tail -n +2 v2.txt) v6.txt || { echo "-6 disagrees with -V 2bit"; exit 1; }

## the three views are not the same answer
diff -q <(tail -1 vset.txt) <(tail -1 vmeth.txt) >/dev/null &&
  { echo "-V set and -V meth returned identical rows"; exit 1; }

## an unknown view is refused
if "$YAME" summary -V nonsense -m m.cg q.cg >/dev/null 2>verr.txt; then
  echo "an unknown -V view was accepted"; exit 1
fi
[ -s verr.txt ] || { echo "an unknown -V failed silently"; exit 1; }

## -V meth on a non-fmt6 query is refused
awk -v n=$n 'BEGIN { for (i = 0; i < n; i++) print (i % 7) "\t" (7 - (i % 7)) }' |
  "$YAME" pack -f m - > mu.cg
if "$YAME" summary -V meth -m m.cg mu.cg >/dev/null 2>merr.txt; then
  echo "-V meth accepted a format-3 query"; exit 1
fi

## ---- 4. no mask at all: the universe plays the mask role ------------------
"$YAME" summary q.cg > nomask.txt 2>/dev/null
grep -q 'global' nomask.txt || { echo "a maskless summary does not report Mask as global"; head -2 nomask.txt; exit 1; }

## ---- 5. a state mask summarises per key ----------------------------------
awk -v n=$n 'BEGIN { split("Alpha Beta Gamma", t, " ")
                     for (i = 0; i < n; i++) print t[(i % 3) + 1] }' > st.txt
"$YAME" pack -f s st.txt > st.cg
"$YAME" summary -m st.cg q.cg > bykey.txt 2>/dev/null
for k in Alpha Beta Gamma; do
  grep -q "$k" bykey.txt || { echo "the state mask did not report $k"; cat bykey.txt; exit 1; }
done
[ "$(tail -n +2 bykey.txt | wc -l)" -ge 3 ] ||
  { echo "a 3-key state mask gave $(tail -n +2 bykey.txt | wc -l) rows"; exit 1; }
## -T keeps the state name in the label even when it is not needed
"$YAME" summary -T -m st.cg q.cg >/dev/null 2>&1 || { echo "summary -T failed"; exit 1; }

## ---- 6. subset of a state track by TERM (-s) ------------------------------
## For format 2, subset selects terms rather than samples.
"$YAME" subset -s st.cg Alpha > alpha.cg 2>/dev/null
[ -s alpha.cg ] || { echo "subset -s produced nothing"; exit 1; }
"$YAME" info alpha.cg >/dev/null 2>&1 || { echo "subset -s output is not readable"; exit 1; }
"$YAME" unpack alpha.cg >/dev/null 2>&1 || { echo "subset -s output cannot be unpacked"; exit 1; }
## a term the track does not carry
if "$YAME" subset -s st.cg NoSuchTerm >/dev/null 2>serr.txt; then
  echo "subset -s accepted an unknown term"; exit 1
fi

## ---- 7. subset -v says which path ran ------------------------------------
cat q.cg q.cg > two.cg
printf 'one\ntwo\n' > nm.txt
"$YAME" index -s nm.txt two.cg
"$YAME" subset -v two.cg one >/dev/null 2>v.txt
grep -qiE 'raw|copy|re-encode|encode' v.txt ||
  { echo "subset -v did not say which path ran"; cat v.txt; exit 1; }
## and -z forces the re-encode path, with the same bytes back
"$YAME" subset two.cg one 2>/dev/null | "$YAME" unpack -f -1 - 2>/dev/null > raw.txt
"$YAME" subset -z1 two.cg one 2>/dev/null | "$YAME" unpack -f -1 - 2>/dev/null > enc.txt
diff raw.txt enc.txt || { echo "the raw and re-encoded subsets differ"; exit 1; }
