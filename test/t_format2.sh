#!/bin/bash
## Format 2, categorical states: the key table, its value widths, and a state
## track as a summary QUERY.
##
## The suite packed small state tracks as masks; nothing packed enough states
## to need 2- or 3-byte values, grew the key table past its first allocation,
## sent an empty line through, or summarized a state query against a fmt6 or
## state mask. Counts below are worked by hand from the seven rows written
## here, and both summary paths must print them.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
cd "$d"

## ---- 1. value width follows the number of states --------------------------
## 300 states need 2 bytes a row and 70000 need 3; both grow the key table
## past its first 1024 slots. (8 bytes starts at 2^24 states -- too many for
## a test.) Keys are scattered, not sorted, so order is what is checked.
for n in 300 70000; do
  awk -v n=$n 'BEGIN { for (i = 0; i < n; i++) print "k" (i * 7919) % n }' > k$n.txt
  "$YAME" pack -f s k$n.txt k$n.cg 2>/dev/null
  cmp -s k$n.txt <("$YAME" unpack k$n.cg 2>/dev/null) ||
    { echo "$n states did not round-trip"; exit 1; }
done
[ "$("$YAME" info k300.cg 2>/dev/null | tail -1 | cut -f6)" -eq 2 ] ||
  { echo "300 states are not stored 2 bytes a row"; exit 1; }
[ "$("$YAME" info k70000.cg 2>/dev/null | tail -1 | cut -f6)" -eq 3 ] ||
  { echo "70000 states are not stored 3 bytes a row"; exit 1; }

## an empty line is the state NA -- an empty key would confuse the separator
printf 'A\n\nB\nA\n' > e.txt
"$YAME" pack -f s e.txt e.cg 2>/dev/null
[ "$("$YAME" unpack e.cg 2>/dev/null | paste -sd' ' -)" = "A NA B A" ] ||
  { echo "an empty line did not become NA"; exit 1; }
## -v says what it stored
"$YAME" pack -v -f s e.txt ev.cg 2>v.txt
grep -q 'Format 2' v.txt || { echo "pack -v -f s said nothing of the format"; exit 1; }

## ---- 2. a state track as the summary query --------------------------------
printf '%s\n' A B C A C B A > st.txt                 # A 1 4 7 . B 2 6 . C 3 5
"$YAME" pack -f s st.txt q2.cg
printf '%s\n' 1 1 0 1 0 1 0 | "$YAME" pack -f b - m0.cm          # rows 1 2 4 6
## set, universe: set AND universe at rows 1 4 6 (row 3 is set, not in it)
printf '1\t1\n0\t1\n1\t0\n1\t1\n0\t0\n1\t1\n0\t1\n' | "$YAME" pack -f d - m6.cm
printf '%s\n' x y x y x y x | "$YAME" pack -f s - m2.cm          # x 1 3 5 7 . y 2 4 6

both() {                          # Query, Mask, N_univ, N_query, N_mask, N_overlap
  "$YAME" summary "$@" 2>/dev/null > k1.txt
  YAME_SUMMARY_KERNEL=0 "$YAME" summary "$@" 2>/dev/null > k0.txt
  diff k1.txt k0.txt >/dev/null ||
    { echo "summary $*: the kernel and the per-mask path disagree"; diff k1.txt k0.txt; exit 1; }
  tail -n +2 k1.txt | cut -f2,4-8
}
want() { printf '%s\n' "$@" | tr ' ' '\t'; }

## no mask: one row per state; the universe is the mask, so overlap = query
diff <(want "A global 7 3 7 3" "B global 7 2 7 2" "C global 7 2 7 2") <(both q2.cg) ||
  { echo "state query, no mask"; exit 1; }
## a binary mask: A 1 4, B 2 6, C none
diff <(want "A 1 7 3 4 2" "B 1 7 2 4 2" "C 1 7 2 4 0") <(both -m m0.cm q2.cg) ||
  { echo "state query, fmt0 mask"; exit 1; }
## a fmt6 mask: the universe stays every row (as for fmt3); A 1 4, B 6, C none
diff <(want "A 1 7 3 3 2" "B 1 7 2 3 1" "C 1 7 2 3 0") <(both -m m6.cm q2.cg) ||
  { echo "state query, fmt6 mask"; exit 1; }
## a state mask: every query state against every mask state
diff <(want "A x 7 3 4 2" "B x 7 2 4 0" "C x 7 2 4 2" "A y 7 3 3 1" "B y 7 2 3 2" "C y 7 2 3 0") \
     <(both -m m2.cm q2.cg) || { echo "state query, state mask"; exit 1; }
## -T puts the record's name in front of each state, query and mask alike,
## whatever the mask (a binary one used to prefix even without -T)
[ "$(both -T q2.cg | cut -f1 | paste -sd' ' -)" = "1-A 1-B 1-C" ] ||
  { echo "-T did not label the query states"; exit 1; }
[ "$(both -T -m m0.cm q2.cg | cut -f1 | paste -sd' ' -)" = "1-A 1-B 1-C" ] ||
  { echo "-T did not label the query states against a binary mask"; exit 1; }
[ "$(both -T -m m2.cm q2.cg | cut -f2 | sort -u | paste -sd' ' -)" = "1-x 1-y" ] ||
  { echo "-T did not label the mask states"; exit 1; }

## a state mask one row short is refused
printf '%s\n' x y x y x y | "$YAME" pack -f s - m2s.cm
if YAME_SUMMARY_KERNEL=0 "$YAME" summary -m m2s.cm q2.cg >/dev/null 2>err.txt; then
  echo "a 6-row state mask was accepted for a 7-row query"; exit 1
fi
grep -q 'different lengths' err.txt || { echo "the short mask failed without saying why"; exit 1; }
echo "ok: t_format2"
