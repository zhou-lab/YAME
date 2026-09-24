#!/bin/bash
## subset: pick samples out of a store by name, in the order asked, or split a
## state track into one bitset per state.
##
## Other tests reach subset only as a way to get one record out; nothing ran
## -o (so the output index had no test), -H/-T, -l, the -s state split, or the
## decode/encode fallback that unaligned stores need. Each value below encodes
## its sample, so a record in the wrong place cannot pass for the right one.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
command -v python3 >/dev/null || { echo "skip: no python3" >&2; exit 0; }
here=$(cd "$(dirname "$0")" && pwd)
d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
cd "$d"

## four samples of 20 rows: sample s holds M = 10*s + row
for s in 1 2 3 4; do
  awk -v s=$s 'BEGIN { for (i = 0; i < 20; i++) print s * 10 + i "\t" i }' > s$s.txt
  "$YAME" pack -f m s$s.txt s$s.cg
done
printf 's1\ns2\ns3\ns4\n' > nm.txt
cat s1.cg s2.cg s3.cg s4.cg > al.cg               # aligned: each record its own block
"$YAME" index -s nm.txt al.cg
python3 "$here/make_unaligned.py" un.cg s1.cg s2.cg s3.cg s4.cg   # one block for all
"$YAME" index -s nm.txt un.cg

## the M column of every record in a store, side by side, one row per CpG
ms() { "$YAME" unpack -a -f -1 "$1" 2>/dev/null |
       awk -F'\t' '{ o = $1; for (k = 3; k <= NF; k += 2) o = o "\t" $k; print o }'; }
## what a subset of the named samples should hold, from the text they came from
want() { for s in "$@"; do cut -f1 $s.txt > $s.col; done
         paste $(for s in "$@"; do printf '%s.col ' $s; done); }

## ---- 1. names, in the order asked, on both paths ---------------------------
for src in al un; do
  "$YAME" subset -v -o o_$src.cg $src.cg s3 s1 2>v.txt
  [ "$(cut -f1 o_$src.cg.idx | paste -sd' ' -)" = "s3 s1" ] ||
    { echo "$src: -o index names $(cut -f1 o_$src.cg.idx | paste -sd' ' -), want s3 s1"; exit 1; }
  diff <(want s3 s1) <(ms o_$src.cg) || { echo "$src: subset s3 s1 did not return s3 then s1"; exit 1; }
  ## the index points at the records it names
  "$YAME" subset o_$src.cg s1 2>/dev/null > back.cg
  diff <(want s1) <(ms back.cg) || { echo "$src: the output index is wrong"; exit 1; }
done
grep -q 'raw block passthrough' <("$YAME" subset -v al.cg s2 2>&1 >/dev/null) ||
  { echo "an aligned store did not take the raw copy"; exit 1; }
grep -q 're-encoding' <("$YAME" subset -v un.cg s2 2>&1 >/dev/null) ||
  { echo "an unaligned store did not fall back to re-encoding"; exit 1; }
## -z forces the re-encode on an aligned store, and reads back the same
"$YAME" subset -z 1 -v -o z.cg al.cg s3 s1 2>v.txt
grep -q 're-encoding' v.txt || { echo "-z did not turn the raw copy off"; exit 1; }
diff <(ms o_al.cg) <(ms z.cg) || { echo "-z 1 changed the data"; exit 1; }

## stdout carries the same records (and, having no file, no index)
"$YAME" subset al.cg s4 2>/dev/null > so.cg
diff <(want s4) <(ms so.cg) || { echo "subset to stdout returned the wrong record"; exit 1; }

## ---- 2. -l, -H, -T ----------------------------------------------------------
printf 's4\ns2\n' > l.txt
"$YAME" subset -l l.txt -o l.cg al.cg 2>/dev/null
[ "$(cut -f1 l.cg.idx | paste -sd' ' -)" = "s4 s2" ] || { echo "-l l.txt: wrong samples"; exit 1; }
"$YAME" subset -l l.txt -o l2.cg al.cg s1 2>/dev/null      # names given: -l ignored
[ "$(cut -f1 l2.cg.idx | paste -sd' ' -)" = "s1" ] || { echo "names did not override -l"; exit 1; }
"$YAME" subset -H 2 -o h.cg al.cg 2>/dev/null
[ "$(cut -f1 h.cg.idx | paste -sd' ' -)" = "s1 s2" ] || { echo "-H 2 is not the first two"; exit 1; }
"$YAME" subset -T 2 -o t.cg un.cg 2>/dev/null
[ "$(cut -f1 t.cg.idx | paste -sd' ' -)" = "s3 s4" ] || { echo "-T 2 is not the last two"; exit 1; }
diff <(want s3 s4) <(ms t.cg) || { echo "-T 2 carried the wrong records"; exit 1; }
"$YAME" subset -H 99 -o h99.cg al.cg 2>/dev/null           # clamped to what there is
[ "$(wc -l < h99.cg.idx)" -eq 4 ] || { echo "-H 99 on 4 samples did not give 4"; exit 1; }
"$YAME" subset -o h0.cg al.cg 2>/dev/null                   # no names at all: the first
[ "$(cut -f1 h0.cg.idx)" = "s1" ] || { echo "no names did not default to the first sample"; exit 1; }

## ---- 3. a name the index does not hold is refused, on both paths -----------
for src in al un; do
  if "$YAME" subset -o bad.cg $src.cg s1 nope >/dev/null 2>err.txt; then
    echo "$src: an unknown sample name was accepted"; exit 1
  fi
  grep -q 'not in the index' err.txt || { echo "$src: an unknown name failed without saying why"; exit 1; }
done

## ---- 4. -s: a state track becomes one bitset per state ---------------------
printf '%s\n' A B C A C B A > st.txt
"$YAME" pack -f s st.txt st.cg
"$YAME" subset -s -o so.cg st.cg A C 2>/dev/null
[ "$(cut -f1 so.cg.idx | paste -sd' ' -)" = "A C" ] || { echo "-s: index names are not the states"; exit 1; }
[ "$("$YAME" info so.cg 2>/dev/null | tail -n +2 | cut -f5 | sort -u)" = "0" ] ||
  { echo "-s did not write format 0"; exit 1; }
awk '{ print ($1 == "A") "\t" ($1 == "C") }' st.txt > so.want
diff so.want <("$YAME" unpack -a so.cg 2>/dev/null) || { echo "-s: a bitset does not follow its state"; exit 1; }
if "$YAME" subset -s s1.cg A >/dev/null 2>err.txt; then
  echo "-s accepted a format 3 input"; exit 1
fi
grep -q 'format 2 state track' err.txt || { echo "-s on format 3 failed without saying why"; exit 1; }
echo "ok: t_subset"
