#!/bin/bash
## index / subset / split: the name-to-offset path every multi-sample store
## is read through. Two shipped bugs lived here -- `index -1` writing one
## constant offset in every conda build (the seek was the argument of an
## assert, which -DNDEBUG deletes), and `subset` leaving a zero-record file
## behind when a name was missing.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
cd "$d"

n=6
for i in $(seq 1 $n); do
  awk -v v=$i 'BEGIN { for (r = 0; r < 30; r++) print v "\t" (10 - v) }' > s$i.txt
  "$YAME" pack -f m s$i.txt > s$i.cg
  echo "sample$i" >> names.txt
done
cat s1.cg s2.cg s3.cg s4.cg s5.cg s6.cg > store.cg

## ---- 1. index -s: one entry per record, all offsets distinct ---------------
"$YAME" index -s names.txt store.cg
[ "$(wc -l < store.cg.idx)" -eq $n ] || { echo "index -s: wrong entry count"; exit 1; }
[ "$(cut -f1 store.cg.idx | sort | diff - <(sort names.txt) >/dev/null && echo ok)" = ok ] ||
  { echo "index -s: names do not match"; exit 1; }
[ "$(cut -f2 store.cg.idx | sort -u | wc -l)" -eq $n ] ||
  { echo "index -s: offsets are not distinct"; cat store.cg.idx; exit 1; }

## ---- 2. index -1: appending one sample at a time gives the same offsets ----
## This is what mergeCG2 calls. Under -DNDEBUG the seek vanished and every
## entry after the second carried the second's address -- 644 of 645 on a real
## store, and `subset` by name then returned the wrong sample.
cp store.cg one.cg
while read -r nm; do "$YAME" index -1 "$nm" one.cg; done < names.txt
[ "$(cut -f2 one.cg.idx | sort -u | wc -l)" -eq $n ] ||
  { echo "index -1: offsets are not distinct"; cat one.cg.idx; exit 1; }
diff <(sort store.cg.idx) <(sort one.cg.idx) ||
  { echo "index -1 disagrees with index -s"; exit 1; }

## ---- 3. subset by name returns that record, on both code paths -------------
## The raw path copies compressed bytes; -z re-encodes. They must agree.
for nm in sample1 sample4 sample6; do
  i=${nm#sample}
  "$YAME" subset store.cg "$nm" > got.cg 2>/dev/null
  "$YAME" unpack -f -1 got.cg 2>/dev/null > got.txt
  diff "s$i.txt" got.txt >/dev/null || { echo "subset $nm returned the wrong record"; exit 1; }
  "$YAME" subset -z6 store.cg "$nm" > gotz.cg 2>/dev/null
  "$YAME" unpack -f -1 gotz.cg 2>/dev/null > gotz.txt
  diff got.txt gotz.txt >/dev/null || { echo "subset -z6 $nm differs from the raw path"; exit 1; }
done

## several names at once, in the order asked for rather than store order
"$YAME" subset store.cg sample5 sample2 > multi.cg 2>/dev/null
"$YAME" unpack -f -1 -a multi.cg 2>/dev/null | cut -f1 > multi.first
[ "$(head -1 multi.first)" = "5" ] && [ "$(head -1 multi.first | wc -l)" -eq 1 ] ||
  { echo "subset with several names lost the requested order"; exit 1; }

## ---- 4. a missing name leaves no output file behind ------------------------
## The check used to sit inside the write loop, so a bad name aborted with the
## output already open, leaving a .cg that reads as a valid EMPTY store --
## indistinguishable downstream from a real one.
rm -f out.cg
if "$YAME" subset -o out.cg store.cg sample2 nosuchsample 2>err.txt; then
  echo "subset accepted a missing name"; exit 1
fi
[ -s err.txt ] || { echo "subset failed silently on a missing name"; exit 1; }
[ ! -e out.cg ] || { echo "subset left an output file behind: $(wc -c < out.cg) bytes"; exit 1; }

## ---- 5. split writes one file per record, named from the index -------------
## Output names are <prefix><name>.cx with -s, and <prefix>_split_<i>.cx
## without. Both forms must hold the right record.
mkdir sp && cd sp
"$YAME" split -s ../names.txt ../store.cg out >/dev/null 2>&1
for i in $(seq 1 $n); do
  f="outsample$i.cx"
  [ -f "$f" ] || { echo "split -s did not produce $f"; /bin/ls; exit 1; }
  "$YAME" unpack -f -1 "$f" 2>/dev/null > s.txt
  diff "../s$i.txt" s.txt >/dev/null || { echo "split -s: $f holds the wrong record"; exit 1; }
done
"$YAME" split ../store.cg pre >/dev/null 2>&1
for i in $(seq 1 $n); do
  f="pre_split_$i.cx"
  [ -f "$f" ] || { echo "split did not produce $f"; /bin/ls; exit 1; }
  "$YAME" unpack -f -1 "$f" 2>/dev/null > s.txt
  diff "../s$i.txt" s.txt >/dev/null || { echo "split: $f holds the wrong record"; exit 1; }
done
cd ..

## ---- 6. a duplicate name is refused ---------------------------------------
printf 'dup\ndup\n' > dupnames.txt
cat s1.cg s2.cg > two.cg
if "$YAME" index -s dupnames.txt two.cg 2>err.txt; then
  echo "index -s accepted a duplicate sample name"; exit 1
fi
