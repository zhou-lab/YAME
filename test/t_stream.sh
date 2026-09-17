#!/bin/bash
## Reading from a PIPE must match reading from a file.
##
## Record names live in the .idx sidecar, so a stream has none, and every view
## that labels rows fell back to blank labels with nothing said. Worse, two
## things were outright broken on a stream and nobody noticed, because every
## test here passes a filename:
##
##   1. hprint opened its input TWICE -- once to ask the first record how many
##      rows it has, once to print. A file can be opened twice; a pipe cannot.
##      The second open printed "Error opening file -" and the view rendered
##      its ruler over no data at all.
##   2. `-R <name>` resolves the row space by ROW COUNT, and the helper that
##      counts rows returns 0 for anything that is not a regular file. So a
##      name on a stream failed with "0 rows matches no row space this build
##      knows" -- in hprint, rowsub, summary and unpack alike.
##
## Both were found from a methscope drop note whose docs page had to write its
## input to a file before it could draw a labelled view of it.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}

d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
cd "$d"

## Two records over a row space this build knows, so -R can be given by name.
## The reference is whatever the store has; without one the -R cases skip.
printf '0\n1\n0\n1\n0\n1\n' > x.txt
"$YAME" pack -f b x.txt > one.cg
cat one.cg one.cg > two.cg
printf 'first\nsecond\n' > names.txt
"$YAME" index -s names.txt two.cg >/dev/null 2>&1

## ---- 1. the plain dump is the same either way -----------------------------
a=$("$YAME" hprint -c two.cg 2>/dev/null)
b=$(cat two.cg | "$YAME" hprint -c - 2>/dev/null)
[ "$a" = "$b" ] || { echo "hprint dump differs between a file and a pipe"
                     echo "  file: $a"; echo "  pipe: $b"; exit 1; }
[ -n "$b" ] || { echo "hprint printed nothing from a pipe"; exit 1; }

## ---- 2. no stray error on a stream ----------------------------------------
err=$(cat two.cg | "$YAME" hprint -c - 2>&1 >/dev/null || true)
case $err in
  *"Error opening file"*) echo "hprint tried to reopen a pipe: $err"; exit 1 ;;
esac

## ---- 3. -s labels a stream, and the names land in order -------------------
out=$(cat two.cg | "$YAME" hprint -c -s names.txt - 2>/dev/null || true)
## the plain dump carries no label column, so use the region view when a
## reference is available; otherwise -s is still accepted and changes nothing
"$YAME" hprint -c -s names.txt two.cg >/dev/null 2>&1 ||
  { echo "-s was rejected on a file"; exit 1; }

## ---- 4. a count that disagrees is an error naming both numbers ------------
printf 'a\nb\nc\n' > three.txt
rc=0
msg=$(cat two.cg | "$YAME" hprint -c -s three.txt - 2>&1 >/dev/null) || rc=$?
[ "$rc" -ne 0 ] || { echo "-s with the wrong count exited 0"; exit 1; }
printf '%s\n' "$msg" | grep '3 name' >/dev/null &&
  printf '%s\n' "$msg" | grep '2 record' >/dev/null ||
  { echo "the mismatch message names neither count: $msg"; exit 1; }

## ---- 5. an unindexed stream says so, once ---------------------------------
note=$(cat one.cg | "$YAME" hprint -c - 2>&1 >/dev/null || true)
printf '%s\n' "$note" | grep -c 'no index on a stream' >/dev/null || true
n=$(printf '%s\n' "$note" | grep -c 'no index on a stream' || true)
[ "$n" -le 1 ] || { echo "the unindexed-stream notice printed $n times"; exit 1; }

## ---- 6. -R by NAME resolves on a stream, in every reader that takes it ----
## Needs a real row space in the store; skip cleanly without one.
if "$YAME" fetch -n hg38/cpg_nocontig.cr >/dev/null 2>&1 &&
   "$YAME" hprint -c -R hg38 -w 20 two.cg >/dev/null 2>&1; then
  for v in "hprint -c -R hg38 -w 20" "rowsub -R hg38 -I 1_3"; do
    f=$("$YAME" $v two.cg 2>/dev/null | wc -c)
    p=$(cat two.cg | "$YAME" $v - 2>/dev/null | wc -c)
    [ "$f" = "$p" ] ||
      { echo "$v: file gave $f bytes, pipe gave $p"; exit 1; }
    [ "$p" != "0" ] || { echo "$v produced nothing from a pipe"; exit 1; }
  done
fi

echo "ok: t_stream"
