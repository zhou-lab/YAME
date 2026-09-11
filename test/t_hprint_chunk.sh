#!/bin/bash
## hprint (the horizontal/ANSI view) and chunk/chunkchar (splitting a stream).
## hprint is the largest untested file in the tree and the one a person looks
## at, so a crash or a shifted alphabet there is visible damage.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
cd "$d"

rows=40
awk -v n=$rows 'BEGIN { for (i = 0; i < n; i++) print (i % 10) "\t" (10 - (i % 10)) }' > mu.txt
"$YAME" pack -f m mu.txt > f3.cg
"$YAME" unpack -f 1 f3.cg 2>/dev/null |
  awk '{ print ($1 == "NA") ? "NA" : $1 }' > beta.txt
"$YAME" pack -f n beta.txt > f4.cg
"$YAME" binarize f3.cg > f6.cg 2>/dev/null
awk -v n=$rows 'BEGIN { for (i = 0; i < n; i++) print (i % 3 == 0) ? 1 : 0 }' |
  "$YAME" pack -f b - > f0.cg

## ---- 1. full-dataset mode accepts fmt0/3/4/6 ------------------------------
## It took format 6 alone until v1.39, so `cat truth.cg reconstructed.cg |
## yame hprint -c -` died on exactly the comparison people wanted.
for f in f0 f3 f4 f6; do
  out=$("$YAME" hprint -c $f.cg 2>err.txt) ||
    { echo "hprint on $f.cg exited non-zero"; cat err.txt; exit 1; }
  [ -n "$out" ] || { echo "hprint on $f.cg printed nothing"; exit 1; }
done

## ---- 2. the decile alphabet is stable -------------------------------------
## -g prints one digit per site. Its glyphs are parsed elsewhere, so the
## alphabet is a contract: digits, '.' for no coverage, nothing else.
"$YAME" hprint -c -g f3.cg 2>/dev/null | tr -d '\n' > g.txt
[ -s g.txt ] || { echo "hprint -g printed nothing"; exit 1; }
LC_ALL=C grep -q '^[0-9.]*$' g.txt ||
  { echo "hprint -g emitted a glyph outside [0-9.]"; cat g.txt; exit 1; }
[ "$(wc -c < g.txt)" -eq $rows ] ||
  { echo "hprint -g printed $(wc -c < g.txt) glyphs, want $rows"; exit 1; }

## the digit must follow the beta: row i has beta (i%10)/10
"$YAME" unpack -f 1 f3.cg 2>/dev/null |
  awk '{ printf "%d", ($1 == "NA") ? 0 : int($1 * 10) }' > g.want
diff <(cat g.want) <(cat g.txt) >/dev/null ||
  { echo "hprint -g digits do not follow the betas"; paste <(cat g.want) <(cat g.txt); exit 1; }

## ---- 3. fmt6 keeps its own alphabet ---------------------------------------
## 1/0/2 is old enough to be parsed elsewhere, so it may not drift to H/M/L.
"$YAME" hprint -c f6.cg 2>/dev/null | tr -d '\n' > six.txt
LC_ALL=C grep -q '^[012]*$' six.txt ||
  { echo "hprint on fmt6 left the 1/0/2 alphabet"; cat six.txt; exit 1; }

## ---- 4. -w clips the view without changing what precedes it ---------------
full=$("$YAME" hprint -c -g f3.cg 2>/dev/null | head -1)
clip=$("$YAME" hprint -c -g -w 10 f3.cg 2>/dev/null | head -1)
[ "${#clip}" -le "${#full}" ] || { echo "-w made the line longer"; exit 1; }
case "$full" in "$clip"*) ;; *) echo "-w 10 is not a prefix of the full view"; exit 1 ;; esac

## ---- 5. chunk: the pieces reassemble into the original --------------------
## Pieces land in <input>_chunks/ as 0.cx, 1.cx, ... in row order, so
## concatenating them in numeric order must give the input back. 40 rows at
## -s 40 is the exact-multiple case that broke the chunked reader.
for s in 7 10 40; do
  rm -rf f3.cg_chunks
  "$YAME" chunk -s $s f3.cg >/dev/null 2>&1 || { echo "chunk -s $s failed"; exit 1; }
  [ -d f3.cg_chunks ] || { echo "chunk -s $s made no output directory"; /bin/ls; exit 1; }
  parts=$(/bin/ls f3.cg_chunks | wc -l)
  [ "$parts" -gt 0 ] || { echo "chunk -s $s produced no pieces"; exit 1; }
  want=$(( (rows + s - 1) / s ))
  [ "$parts" -eq "$want" ] ||
    { echo "chunk -s $s made $parts pieces, want $want"; exit 1; }
  ## Each piece is its own RECORD, so concatenating the files gives a
  ## multi-record store rather than one long record. The property that holds
  ## is on the text: the pieces, unpacked in order, are the original rows.
  : > rejoined.txt
  for part in $(/bin/ls f3.cg_chunks | sort -n); do
    "$YAME" unpack -f -1 "f3.cg_chunks/$part" 2>/dev/null >> rejoined.txt
  done
  diff mu.txt rejoined.txt >/dev/null ||
    { echo "chunk -s $s did not reassemble to the original"; exit 1; }
done
rm -rf f3.cg_chunks

## ---- 6. chunkchar splits text on line boundaries --------------------------
## No line may be cut in half, and the pieces must reassemble exactly.
awk 'BEGIN { for (i = 0; i < 25; i++) print "line" i }' > t.txt
"$YAME" chunkchar -s 10 t.txt >/dev/null 2>&1 || { echo "chunkchar failed"; exit 1; }
[ -d t.txt_chunks ] || { echo "chunkchar made no output directory"; /bin/ls; exit 1; }
cat $(/bin/ls t.txt_chunks | sort -n | sed 's|^|t.txt_chunks/|') > t.rejoined
diff t.txt t.rejoined || { echo "chunkchar did not reassemble"; exit 1; }
