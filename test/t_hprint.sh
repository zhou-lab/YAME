#!/bin/bash
## hprint's views for binary data: format 6 per CpG and genome-wide, format 0
## windowed and genome-wide, and the refusals. t_coords.sh covers the same
## views over M/U counts; these formats take their own paths through them.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
cd "$d"

## two chromosomes of 10 CpGs, 100 bp apart
awk 'BEGIN { split("chr1 chr2", c, " ")
  for (k = 1; k <= 2; k++) for (i = 0; i < 10; i++)
    printf "%s\t%d\t%d\n", c[k], 100 + i * 100, 102 + i * 100 }' > ref.bed
"$YAME" pack -f r ref.bed ref.cr
## fmt6 over chr1: set/universe 11 01 10 00 11 01 11 11 01 00 -> meth, unmeth,
## NA (set outside the universe is not a call), NA, then █░██░. ; chr2 all meth
printf '1\t1\n0\t1\n1\t0\n0\t0\n1\t1\n0\t1\n1\t1\n1\t1\n0\t1\n0\t0\n' > a.txt
for i in $(seq 10); do printf '1\t1\n'; done >> a.txt
"$YAME" pack -f d a.txt q6.cg
## fmt0: every third row set
awk 'BEGIN { for (i = 0; i < 20; i++) print (i % 3 == 0) }' | "$YAME" pack -f b - q0.cg

## ---- 1. fmt6, one column per CpG, in a region and genome-wide -------------
[ "$("$YAME" hprint -c -R ref.cr -r chr1 q6.cg 2>/dev/null | tail -1 | tr -d ' ')" = "█░..█░██░." ] ||
  { echo "fmt6 region view"; "$YAME" hprint -c -R ref.cr -r chr1 q6.cg; exit 1; }
[ "$("$YAME" hprint -c -R ref.cr q6.cg 2>/dev/null | tail -1 | tr -d ' ')" = "█░..█░██░.██████████" ] ||
  { echo "fmt6 whole-genome view"; exit 1; }

## ---- 2. fmt0, window-averaged ----------------------------------------------
## chr1's rows 0-9 are 1001001001; in windows of 4,4,2 that is 1001 0010 01
## -> 0.5 0.25 0.5 -> M L M (0.33-0.67 is M). Genome-wide, -w 4 does the same
## over each chromosome.
[ "$("$YAME" hprint -c -R ref.cr -r chr1 -w 3 q0.cg 2>/dev/null | tail -1 | tr -d ' ')" = "MLM" ] ||
  { echo "fmt0 windowed region view"; "$YAME" hprint -c -R ref.cr -r chr1 -w 3 q0.cg; exit 1; }
[ "$("$YAME" hprint -c -R ref.cr -w 4 q0.cg 2>/dev/null | tail -1 | tr -d ' ')" = "MMLM" ] ||
  { echo "fmt0 whole-genome view"; exit 1; }

## ---- 3. what hprint refuses -------------------------------------------------
refuse() {      # <what the message must say> <args...>
  local msg=$1; shift
  if "$YAME" hprint -c "$@" >/dev/null 2>err.txt; then echo "hprint $* was accepted"; exit 1; fi
  grep -q "$msg" err.txt || { echo "hprint $* failed without saying '$msg'"; cat err.txt; exit 1; }
}
refuse 'Cannot parse region' -R ref.cr -r 'chr1:abc' q6.cg
refuse 'must be format 7' -R q0.cg -r chr1 q6.cg
refuse 'must be format 7' -R q0.cg q6.cg
awk 'BEGIN { for (i = 0; i < 19; i++) print 1 }' | "$YAME" pack -f b - q19.cg
refuse 'Dimension mismatch' -R ref.cr -r chr1 q19.cg
refuse 'Dimension mismatch' -R ref.cr q19.cg
## -s with a name count the records do not have: the view streams, so the
## rows are already out, but the run still fails and says why
printf 'a\nb\n' > two.txt
refuse 'gave 2 names but the input holds 1 record' -s two.txt q6.cg
echo "ok: t_hprint"
