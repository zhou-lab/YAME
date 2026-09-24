#!/bin/bash
## Format 3, M/U counts: the lossy fit, the zero-run coding, and the per-mask
## summary path a format 6 mask takes.
##
## fitMU -- halve M and U together until both fit, so the ratio survives and
## the absolute depth does not -- is documented and was never exercised: no
## fixture held a count past the width it was packed at. Nor did any test run
## the per-mask oracle (YAME_SUMMARY_KERNEL=0) for a fmt3 query against a
## fmt6 mask, the arithmetic the one-pass kernel is checked against.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
cd "$d"

## ---- 1. fitMU: a count past its width is halved, with its partner -----------
## -u 1 is one byte a site, 4 bits each for M and U, so both must end below
## 15: 100/50 halves three times to 12/6. 14/1 already fits and stays.
printf '100\t50\n3\t4\n14\t1\n' > u.txt
"$YAME" pack -u 1 -f m u.txt u1.cg 2>/dev/null
[ "$("$YAME" unpack -f -1 u1.cg 2>/dev/null | paste -sd' ' -)" = "$(printf '12\t6 3\t4 14\t1')" ] ||
  { echo "-u 1 did not fit 100/50 to 12/6 and leave the rest"; "$YAME" unpack -f -1 u1.cg; exit 1; }
## past 31 bits the compressed form halves too: 5e9/2.5e9 -> 1.25e9/6.25e8
printf '5000000000\t2500000000\n7\t9\n' > big.txt
"$YAME" pack -f m big.txt big.cg 2>/dev/null
[ "$("$YAME" unpack -f -1 big.cg 2>/dev/null | head -1)" = "$(printf '1250000000\t625000000')" ] ||
  { echo "a count past 2^31 was not halved to fit"; exit 1; }
## the ratio is what survives
[ "$("$YAME" unpack -f 1 big.cg 2>/dev/null | head -1)" = "0.667" ] ||
  { echo "fitting a large count changed its beta"; exit 1; }

## ---- 2. runs of uncovered sites, in every position -------------------------
printf '0\t0\n0\t0\n1\t1\n0\t0\n5\t5\n0\t0\n0\t0\n0\t0\n2\t0\n0\t0\n' > z.txt
"$YAME" pack -f m z.txt z.cg 2>/dev/null
diff z.txt <("$YAME" unpack -f -1 z.cg 2>/dev/null) || { echo "zero runs did not round-trip"; exit 1; }
"$YAME" pack -v -f m z.txt zv.cg 2>v.txt
grep -q 'Format 3' v.txt || { echo "pack -v -f m said nothing of the format"; exit 1; }

## ---- 3. what pack refuses ----------------------------------------------------
if printf '1\n' | "$YAME" pack -f m - x.cg 2>err.txt; then echo "pack -f m took one field"; exit 1; fi
grep -q 'fewer than 2 fields' err.txt || { echo "a one-field line failed without saying why"; exit 1; }
if printf 'a\t1\n' | "$YAME" pack -f m - x.cg 2>err.txt; then echo "pack -f m took a non-number"; exit 1; fi
grep -q 'not a nonnegative integer' err.txt || { echo "a non-number failed without saying why"; exit 1; }

## ---- 4. a fmt3 query against a fmt6 mask, on both summary paths ------------
## covered: 1 3 4 5 6 7 (M/U 3/1 5/5 2/8 0/4 6/0 1/1); the mask is set AND in
## its universe at 1 4 6: betas 0.75 0.2 1.0 (mean 0.650), depths 4 10 6 (6.667)
printf '3\t1\n0\t0\n5\t5\n2\t8\n0\t4\n6\t0\n1\t1\n' | "$YAME" pack -f m - q3.cg
printf '1\t1\n0\t1\n1\t0\n1\t1\n0\t0\n1\t1\n0\t1\n' | "$YAME" pack -f d - m6.cm
for k in 1 0; do
  [ "$(YAME_SUMMARY_KERNEL=$k "$YAME" summary -m m6.cm q3.cg 2>/dev/null | tail -1 | cut -f5-8,10,11)" = \
    "$(printf '7\t6\t3\t3\t0.650\t6.667')" ] || { echo "fmt3 query, fmt6 mask (kernel=$k)"; exit 1; }
done
## a mask one row short, of both kinds, is refused by the per-mask path
printf '1\t1\n0\t1\n1\t0\n' | "$YAME" pack -f d - m6s.cm
printf '%s\n' a b a | "$YAME" pack -f s - m2s.cm
for m in m6s m2s; do
  if YAME_SUMMARY_KERNEL=0 "$YAME" summary -m $m.cm q3.cg >/dev/null 2>err.txt; then
    echo "a short $m mask was accepted"; exit 1
  fi
  grep -q 'different lengths' err.txt || { echo "the short $m mask failed without saying why"; exit 1; }
done
echo "ok: t_format3"
