#!/bin/bash
## Format 5, the obsolete one.
##
## yame still DECODES it -- "still decodable by unpack, no longer packable",
## as the docs put it -- and nothing can write one any more, so its decoder
## had no coverage at all and no way to get any. test/make_fmt5.py writes the
## bytes by hand, exactly as format5.c documents them, which is enough to
## exercise the decoder and to prove the claim in the docs is still true.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
command -v python3 >/dev/null || { echo "skip: no python3" >&2; exit 0; }
here=$(cd "$(dirname "$0")" && pwd)
d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
cd "$d"

## ---- 1. a hand-written record decodes to what was encoded -----------------
python3 "$here/make_fmt5.py" f5.cg > want.line
tr ' ' '\n' < want.line | sed '/^$/d' > want.txt
[ "$("$YAME" info f5.cg 2>/dev/null | tail -1 | cut -f5)" = "5" ] ||
  { echo "the crafted record does not report format 5"; "$YAME" info f5.cg; exit 1; }
"$YAME" unpack f5.cg 2>/dev/null > got.txt
diff want.txt got.txt || { echo "format 5 did not decode to what was encoded"; exit 1; }

## ---- 2. the run encodings at their edges ----------------------------------
## a long NA run (the length field is 7 bits, so 127 is its ceiling), a vector
## of only 0/1, one of only NA, and a single value
python3 "$here/make_fmt5.py" long.cg $(python3 -c 'print(" ".join(["2"]*127 + ["1","0"]))') > /dev/null
[ "$("$YAME" unpack long.cg 2>/dev/null | wc -l)" -eq 129 ] ||
  { echo "a 127-long NA run did not decode to 129 values"; exit 1; }
[ "$("$YAME" unpack long.cg 2>/dev/null | grep -c NA)" -eq 127 ] ||
  { echo "the NA run lost values"; exit 1; }

python3 "$here/make_fmt5.py" bits.cg 1 0 1 1 0 0 1 0 1 > /dev/null
[ "$("$YAME" unpack bits.cg 2>/dev/null | paste -sd, -)" = "1,0,1,1,0,0,1,0,1" ] ||
  { echo "a pure 0/1 vector did not round-trip"; "$YAME" unpack bits.cg | paste -sd, -; exit 1; }

python3 "$here/make_fmt5.py" nas.cg 2 2 2 2 > /dev/null
[ "$("$YAME" unpack nas.cg 2>/dev/null | sort -u)" = "NA" ] ||
  { echo "an all-NA vector did not decode as NA"; exit 1; }

python3 "$here/make_fmt5.py" one.cg 1 > /dev/null
[ "$("$YAME" unpack one.cg 2>/dev/null)" = "1" ] || { echo "a single value did not decode"; exit 1; }

## ---- 3. the commands that walk any record still work on it ----------------
"$YAME" info f5.cg >/dev/null 2>&1 || { echo "info failed on format 5"; exit 1; }
"$YAME" rowsub -B 3_8 f5.cg > sub.cg 2>/dev/null
"$YAME" unpack sub.cg 2>/dev/null > sub.txt
sed -n '4,8p' want.txt > sub.want
diff sub.want sub.txt || { echo "rowsub on format 5 returned the wrong rows"; exit 1; }

## it can be carried through a store beside other formats
printf '1\n0\n1\n0\n1\n0\n1\n0\n1\n0\n1\n0\n1\n0\n1\n0\n' | "$YAME" pack -f b - > f0.cg
cat f5.cg f0.cg > mixed.cg
printf 'five\nzero\n' > nm.txt
"$YAME" index -s nm.txt mixed.cg
[ "$("$YAME" info mixed.cg 2>/dev/null | tail -n +2 | cut -f5 | paste -sd, -)" = "5,0" ] ||
  { echo "a store holding format 5 beside format 0 misreports"; "$YAME" info mixed.cg; exit 1; }
"$YAME" subset mixed.cg five 2>/dev/null | "$YAME" unpack - 2>/dev/null > back.txt
diff want.txt back.txt || { echo "format 5 did not survive subset out of a mixed store"; exit 1; }

## ---- 4. and it is still not packable, which is the other half of the claim -
if "$YAME" pack -f 5 want.txt >/dev/null 2>&1; then
  echo "pack accepted -f 5; format 5 is supposed to be unwritable"; exit 1
fi
