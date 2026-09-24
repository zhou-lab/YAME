#!/bin/bash
## summary with a format 0 or format 4 QUERY, against every mask kind.
##
## Format 3 and 6 queries are what the suite and the knowledgebase regression
## exercise; a binary (fmt0) or beta (fmt4) query went through its own
## summarize1_queryfmt0/4 with no test at all -- no mask, a fmt6 mask and a
## state mask were never run. Every expected number below was worked out by
## hand from the 12 rows written here, and both summary paths (the one-pass
## kernel and YAME_SUMMARY_KERNEL=0, the per-mask oracle) must print it.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
cd "$d"

## ---- fixtures: 12 rows -----------------------------------------------------
## -v: the readers report what they loaded -- the only reader output there is
printf '%s\n' 1 0 1 1 0 0 1 0 1 1 0 1 > q0.txt            # in the set: 1 3 4 7 9 10 12
"$YAME" pack -v -f b q0.txt q0.cg 2>v.txt
grep -q 'Vector of length 12 loaded' v.txt || { echo "pack -v -f b: no load report"; exit 1; }
printf '%s\n' 0.1 NA 0.5 0.9 NA 0.3 0.7 0.2 NA 1.0 0.4 0.6 > q4.txt
"$YAME" pack -v -f n q4.txt q4.cg 2>v.txt                 # NA at 2 5 9
grep -q 'Vector of length 12 loaded' v.txt || { echo "pack -v -f n: no load report"; exit 1; }
printf '%s\n' 1 1 0 1 0 1 0 0 1 0 1 0 > m0.txt            # 1 2 4 6 9 11
"$YAME" pack -f b m0.txt m0.cm
## set, universe: universe = 1 2 4 6 7 8 10 11 12; set within it = 1 4 6 8 10
## (row 3 is set but outside the universe, so it must not count)
printf '1\t1\n0\t1\n1\t0\n1\t1\n0\t0\n1\t1\n0\t1\n1\t1\n0\t0\n1\t1\n0\t1\n0\t1\n' > m6.txt
"$YAME" pack -f d m6.txt m6.cm
printf '%s\n' A A B B C C A B C A B C > m2.txt            # A 1 2 7 10, B 3 4 8 11, C 5 6 9 12
"$YAME" pack -f s m2.txt m2.cm

## run summary both ways; the two must agree before either is checked
both() {
  "$YAME" summary "$@" 2>/dev/null > k1.txt
  YAME_SUMMARY_KERNEL=0 "$YAME" summary "$@" 2>/dev/null > k0.txt
  diff k1.txt k0.txt >/dev/null ||
    { echo "summary $*: the kernel and the per-mask path disagree"; diff k1.txt k0.txt; exit 1; }
  ## Mask, N_univ, N_query, N_mask, N_overlap, Beta
  tail -n +2 k1.txt | cut -f4-8,10
}
want() { printf '%s\n' "$@" | tr ' ' '\t'; }

## ---- 1. a binary query is a set ----------------------------------------------
## A set has no methylation level, so Beta is NA. No mask: the universe plays
## the mask role, for every query format, so N_mask is the universe.
diff <(want "global 12 7 12 7 NA") <(both q0.cg) ||
  { echo "fmt0 query, no mask"; exit 1; }
diff <(want "1 12 7 6 3 NA") <(both -m m0.cm q0.cg) ||
  { echo "fmt0 query, fmt0 mask"; exit 1; }
## a fmt6 mask's universe bounds EVERYTHING for a set query: 9 rows, of which
## the query holds 1 4 7 10 12 and the mask 1 4 6 8 10
diff <(want "1 9 5 5 3 NA") <(both -m m6.cm q0.cg) ||
  { echo "fmt0 query, fmt6 mask"; exit 1; }
diff <(want "A 12 7 4 3 NA" "B 12 7 4 2 NA" "C 12 7 4 2 NA") \
     <(both -m m2.cm q0.cg) || { echo "fmt0 query, state mask"; exit 1; }

## ---- 2. a beta query: counts of non-NA rows, and their mean --------------------
## nine non-NA betas summing to 4.7
diff <(want "global 12 9 12 9 0.522") <(both q4.cg) ||
  { echo "fmt4 query, no mask"; exit 1; }
## mask rows 1 2 4 6 9 11: non-NA 0.1 0.9 0.3 0.4
diff <(want "1 12 9 6 4 0.425") <(both -m m0.cm q4.cg) ||
  { echo "fmt4 query, fmt0 mask"; exit 1; }
## as for fmt3, the universe stays every row; the mask is set AND in universe:
## rows 1 4 6 8 10, betas 0.1 0.9 0.3 0.2 1.0
diff <(want "1 12 9 5 5 0.500") <(both -m m6.cm q4.cg) ||
  { echo "fmt4 query, fmt6 mask"; exit 1; }
## A: 0.1 0.7 1.0; B: 0.5 0.9 0.2 0.4; C: 0.3 0.6
diff <(want "A 12 9 4 3 0.600" "B 12 9 4 4 0.500" "C 12 9 4 2 0.450") \
     <(both -m m2.cm q4.cg) || { echo "fmt4 query, state mask"; exit 1; }

## -T prefixes each state with the mask's own name (here the default, 1)
for q in q0 q4; do
  [ "$(both -T -m m2.cm $q.cg | cut -f1 | paste -sd' ' -)" = "1-A 1-B 1-C" ] ||
    { echo "$q.cg: -T did not label the states 1-A 1-B 1-C"; exit 1; }
done

## ---- 2b. no mask means the universe is the mask, for every query format ------
## (formats 2, 3, 4 and 7 used to report N_mask 0 here while 0, 1 and 6
## reported the universe; the global row now reads the same way everywhere)
printf '3\t1\n0\t0\n2\t2\n' | "$YAME" pack -f m - q3.cg          # 2 of 3 covered
diff <(want "global 3 2 3 2") <(both q3.cg | cut -f1-5) || { echo "fmt3 query, no mask"; exit 1; }
printf '1\t1\n0\t1\n0\t0\n' | "$YAME" pack -f d - q6.cg          # universe 2, set 1
diff <(want "global 2 1 2 1") <(both q6.cg | cut -f1-5) || { echo "fmt6 query, no mask"; exit 1; }
## a state query: one row per state, the state's count is query and overlap
diff <(want "12 4 12 4" "12 4 12 4" "12 4 12 4") <(both m2.cm | cut -f2-5) ||
  { echo "fmt2 query, no mask"; exit 1; }
## coordinates: one row per chromosome
printf 'chr1\t10\t12\nchr1\t20\t22\nchr2\t5\t7\n' | "$YAME" pack -f r - q7.cg
diff <(want "3 2 3 2" "3 1 3 1") <(both q7.cg | cut -f2-5) || { echo "fmt7 query, no mask"; exit 1; }

## ---- 3. what must be refused -------------------------------------------------
## a mask one row short, of each kind, for each query -- on both paths
head -11 m0.txt | "$YAME" pack -f b - m0s.cm
head -11 m6.txt | "$YAME" pack -f d - m6s.cm
head -11 m2.txt | "$YAME" pack -f s - m2s.cm
for q in q0 q4; do
  for m in m0s m6s m2s; do
    for k in 1 0; do
      if YAME_SUMMARY_KERNEL=$k "$YAME" summary -m $m.cm $q.cg >/dev/null 2>err.txt; then
        echo "$q.cg against the 11-row $m.cm was accepted (kernel=$k)"; exit 1
      fi
      grep -q 'different lengths' err.txt ||
        { echo "$q.cg vs $m.cm (kernel=$k) failed without saying why"; cat err.txt; exit 1; }
    done
  done
  ## a beta track is not a mask
  if YAME_SUMMARY_KERNEL=0 "$YAME" summary -m q4.cg $q.cg >/dev/null 2>err.txt; then
    echo "$q.cg accepted a format 4 mask"; exit 1
  fi
  grep -q 'Mask format 4 unsupported' err.txt ||
    { echo "$q.cg vs a format 4 mask failed without saying why"; cat err.txt; exit 1; }
done
echo "ok: t_summary_formats"
