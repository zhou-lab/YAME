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

## ---- 6. unpack's four ways of choosing records ----------------------------
## By name, the first N, the last N (needs the index), and everything. Each
## has its own reader, and they must agree about which records they name.
"$YAME" unpack -a -f -1 store.cg 2>/dev/null | awk -F'\t' '{print NF}' | sort -u > ncols.txt
[ "$(cat ncols.txt)" = "12" ] || { echo "-a on 6 M/U samples gave $(cat ncols.txt) columns"; exit 1; }

## named samples, in the order asked for
"$YAME" unpack -f -1 store.cg sample4 sample1 2>/dev/null | cut -f1,2 > byname.txt
paste <(cut -f1 s4.txt) <(cut -f2 s4.txt) > byname.want
diff byname.want byname.txt || { echo "unpack by name did not put sample4 first"; exit 1; }

## -l takes the same names from a file
printf 'sample4\nsample1\n' > pick.txt
"$YAME" unpack -f -1 -l pick.txt store.cg 2>/dev/null | cut -f1,2 > bylist.txt
diff byname.txt bylist.txt || { echo "unpack -l disagrees with naming the samples"; exit 1; }

## -H N is the first N records, -T N the last N (the tail reader needs the index)
"$YAME" unpack -f -1 -H 2 store.cg 2>/dev/null | awk -F'\t' '{print NF}' | sort -u > h.txt
[ "$(cat h.txt)" = "4" ] || { echo "-H 2 gave $(cat h.txt) columns, want 4"; exit 1; }
"$YAME" unpack -f -1 -H 2 store.cg 2>/dev/null | cut -f1,2 > head2.txt
paste <(cut -f1 s1.txt) <(cut -f2 s1.txt) > head2.want
diff head2.want head2.txt || { echo "-H 2 did not start at the first record"; exit 1; }

"$YAME" unpack -f -1 -T 2 store.cg 2>/dev/null | awk -F'\t' '{print NF}' | sort -u > t.txt
[ "$(cat t.txt)" = "4" ] || { echo "-T 2 gave $(cat t.txt) columns, want 4"; exit 1; }
"$YAME" unpack -f -1 -T 2 store.cg 2>/dev/null | cut -f3,4 > tail2.txt
paste <(cut -f1 s6.txt) <(cut -f2 s6.txt) > tail2.want
diff tail2.want tail2.txt || { echo "-T 2 did not end at the last record"; exit 1; }

## -H beyond the record count is the whole file, not an error
"$YAME" unpack -f -1 -H 99 store.cg 2>/dev/null | awk -F'\t' '{print NF}' | sort -u > hbig.txt
[ "$(cat hbig.txt)" = "12" ] || { echo "-H 99 on 6 records gave $(cat hbig.txt) columns"; exit 1; }

## a name the store does not have
if "$YAME" unpack store.cg nosuchsample >/dev/null 2>uerr.txt; then
  echo "unpack accepted a name the store does not have"; exit 1
fi
[ -s uerr.txt ] || { echo "unpack failed silently on an unknown name"; exit 1; }

## ---- 7. a duplicate name is refused ---------------------------------------
printf 'dup\ndup\n' > dupnames.txt
cat s1.cg s2.cg > two.cg
if "$YAME" index -s dupnames.txt two.cg 2>err.txt; then
  echo "index -s accepted a duplicate sample name"; exit 1
fi
