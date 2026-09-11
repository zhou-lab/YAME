#!/bin/bash
## Every command, over every format it accepts.
##
## The format codecs are where the bugs that nothing reached lived: fmt4
## indexed as float_t, the fmt1 run-length overread, unaligned uint16 stores,
## the FMT6 macro dropping its index, cdata_nbytes ignoring unit. Each of
## those sat in a branch no COMMAND exercised, because the suite only ever ran
## one or two formats through each verb. This runs the grid.
##
## For each (command, format) the rule is the same as everywhere else: the
## command either produces a readable result of the right shape, or exits
## non-zero saying why. Exit 0 with something unreadable is the failure.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
cd "$d"

n=32                            # divides by 8 (fmt0/fmt6 byte boundaries)

## ---- one record per format, all n rows ------------------------------------
awk -v n=$n 'BEGIN { for (i = 0; i < n; i++) print (i % 3 == 0) ? 1 : 0 }' > f0.txt
awk -v n=$n 'BEGIN { for (i = 0; i < n; i++) printf "%c\n", 97 + (i % 5) }' > f1.txt
awk -v n=$n 'BEGIN { split("Promoter Exon Intron Intergenic", t, " ")
                     for (i = 0; i < n; i++) print t[(i % 4) + 1] }' > f2.txt
awk -v n=$n 'BEGIN { for (i = 0; i < n; i++) {
                       if (i % 7 == 3) { print 0 "\t" 0; continue }
                       m = i % 11; print m "\t" (10 - m) } }' > f3.txt
awk -v n=$n 'BEGIN { for (i = 0; i < n; i++)
                       print (i % 6 == 1) ? "NA" : sprintf("%.3f", (i % 10) / 10) }' > f4.txt
awk -v n=$n 'BEGIN { for (i = 0; i < n; i++)
                       print ((i % 5 == 2) ? "0\t0" : (((i % 3) ? "1" : "0") "\t1")) }' > f6.txt
awk -v n=$n 'BEGIN { for (i = 0; i < n; i++) {
                       c = (i < n / 2) ? "chr1" : "chr2"
                       b = 100 + (i % (n / 2)) * 50; printf "%s\t%d\t%d\n", c, b, b + 2 } }' > f7.txt

"$YAME" pack -f b f0.txt > f0.cg
"$YAME" pack -f c f1.txt > f1.cg
"$YAME" pack -f s f2.txt > f2.cg
"$YAME" pack -f m f3.txt > f3.cg
"$YAME" pack -f n f4.txt > f4.cg
"$YAME" pack -f d f6.txt > f6.cg
"$YAME" pack -f r f7.txt > f7.cg
FMTS="0 1 2 3 4 6 7"

for f in $FMTS; do
  got=$("$YAME" info f$f.cg 2>/dev/null | tail -1 | cut -f5)
  [ "$got" = "$f" ] || { echo "f$f.cg reports format $got"; exit 1; }
  rows=$("$YAME" info f$f.cg 2>/dev/null | tail -1 | cut -f4)
  [ "$rows" -eq $n ] || { echo "f$f.cg has $rows rows, want $n"; exit 1; }
done

## A command run on one format: it must either succeed and leave a readable
## store of `want` rows, or fail with a message. Never exit 0 with garbage.
## `want` of -1 means "same as the input".
check() {                       # <label> <expect-rows> <cmd...>
  local label=$1 want=$2; shift 2
  local rc=0
  "$@" > out.cg 2> out.err || rc=$?
  if [ "$rc" -ne 0 ]; then
    [ -s out.err ] || { echo "$label: exit $rc with an empty stderr"; exit 1; }
    return 0                    # a clean refusal is an acceptable answer
  fi
  "$YAME" info out.cg > out.info 2>/dev/null ||
    { echo "$label: exit 0 but the output is not a readable store"; exit 1; }
  local got; got=$(tail -n +2 out.info | head -1 | cut -f4)
  if [ "$want" -ge 0 ] && [ "$got" != "$want" ]; then
    echo "$label: produced $got rows, want $want"; exit 1
  fi
  ## and it can be printed back without dying
  "$YAME" unpack out.cg >/dev/null 2>&1 ||
    { echo "$label: output cannot be unpacked"; exit 1; }
}

## ---- row selection, every format ------------------------------------------
printf '1\n9\n32\n' > idx3.txt
for f in $FMTS; do
  check "rowsub -B f$f"  8 "$YAME" rowsub -B 8_16 f$f.cg
  check "rowsub -I f$f"  8 "$YAME" rowsub -I 1_8  f$f.cg
  check "rowsub -l f$f"  3 "$YAME" rowsub -l idx3.txt f$f.cg
done

## All four selectors must agree on the same rows, for every format. This is
## what caught -l on format 2: a state track carries its key table in front of
## the rows, and only -l indexed straight into it -- `info` accepted the
## result and `unpack` called it corrupted, exit 0 either way.
printf '9\n10\n11\n12\n13\n14\n15\n16\n' > idx8.txt
awk -v n=$n 'BEGIN { for (i = 0; i < n; i++) print (i >= 8 && i < 16) ? 1 : 0 }' |
  "$YAME" pack -f b - > mid8.cg
for f in $FMTS; do
  "$YAME" rowsub -B 8_16    f$f.cg 2>/dev/null | "$YAME" unpack - 2>/dev/null > by_b.txt
  "$YAME" rowsub -I 1_8     f$f.cg 2>/dev/null | "$YAME" unpack - 2>/dev/null > by_i.txt
  "$YAME" rowsub -l idx8.txt f$f.cg 2>/dev/null | "$YAME" unpack - 2>/dev/null > by_l.txt
  "$YAME" rowsub -m mid8.cg f$f.cg 2>/dev/null | "$YAME" unpack - 2>/dev/null > by_m.txt
  for other in by_i by_l by_m; do
    diff by_b.txt $other.txt >/dev/null ||
      { echo "f$f: rowsub selectors disagree on rows 8-15 (-B vs ${other#by_})"
        diff by_b.txt $other.txt | head -4; exit 1; }
  done
  [ "$(wc -l < by_b.txt)" -eq 8 ] || { echo "f$f: the common slice is not 8 rows"; exit 1; }
done

## a mask keeps the same rows whatever the payload format
awk -v n=$n 'BEGIN { for (i = 0; i < n; i++) print (i % 4 == 0) ? 1 : 0 }' |
  "$YAME" pack -f b - > msk.cg
for f in $FMTS; do check "rowsub -m f$f" 8 "$YAME" rowsub -m msk.cg f$f.cg; done

## ---- unpack, every format and every print mode ----------------------------
for f in $FMTS; do
  for mode in "" "-f 0" "-f 1" "-f -1" "-f 5"; do
    "$YAME" unpack $mode f$f.cg >/dev/null 2>&1 ||
      { echo "unpack $mode on f$f exited non-zero"; exit 1; }
  done
  ## Chunked output equals unchunked, at a size that divides the row count and
  ## one that does not -- except for formats 2 and 7, which are not flat row
  ## vectors (a state track carries a key table, a .cr is a delta stream) and
  ## are refused by name rather than sliced positionally.
  "$YAME" unpack f$f.cg 2>/dev/null > plain.txt
  for s in 8 7; do
    if [ "$f" = 2 ] || [ "$f" = 7 ]; then
      if "$YAME" unpack -c -s $s f$f.cg >/dev/null 2>ch.err; then
        echo "unpack -c on f$f succeeded; it is not a flat row vector"; exit 1
      fi
      grep -qiE 'rowsub|does not support' ch.err ||
        { echo "unpack -c on f$f refused without saying why"; cat ch.err; exit 1; }
    else
      "$YAME" unpack -c -s $s f$f.cg 2>/dev/null > chunk.txt
      diff plain.txt chunk.txt >/dev/null ||
        { echo "unpack -c -s $s on f$f differs from unchunked"; exit 1; }
    fi
  done
done

## ---- a store holding every format at once ---------------------------------
## read_cdata2 walks it through ONE cdata_t, so a stale unit, key table or aux
## from the previous record shows up here and nowhere else.
cat f0.cg f1.cg f2.cg f3.cg f4.cg f6.cg f7.cg > all.cg
printf 'z0\nz1\nz2\nz3\nz4\nz6\nz7\n' > allnames.txt
"$YAME" index -s allnames.txt all.cg
[ "$("$YAME" info all.cg 2>/dev/null | tail -n +2 | cut -f5 | paste -sd,)" = "0,1,2,3,4,6,7" ] ||
  { echo "the mixed store does not report one format per record"; "$YAME" info all.cg; exit 1; }
## each record still round-trips out of the middle of it
for f in $FMTS; do
  "$YAME" subset all.cg "z$f" > one.cg 2>/dev/null
  [ "$("$YAME" info one.cg 2>/dev/null | tail -1 | cut -f5)" = "$f" ] ||
    { echo "subset z$f out of the mixed store changed its format"; exit 1; }
  case $f in
    3) "$YAME" unpack -f -1 one.cg 2>/dev/null > r.txt; diff f3.txt r.txt >/dev/null || { echo "z3 changed"; exit 1; } ;;
    ## fmt6 prints an out-of-universe site as NA<TAB>0 rather than 0<TAB>0 --
    ## "not set" and "not asked" are different answers -- so compare against
    ## the file's own first unpack rather than the packed text.
    6) "$YAME" unpack -f -1 one.cg 2>/dev/null > r.txt
       "$YAME" unpack -f -1 f6.cg 2>/dev/null > w.txt
       diff w.txt r.txt >/dev/null || { echo "z6 changed"; exit 1; } ;;
    *) "$YAME" unpack one.cg 2>/dev/null > r.txt; diff f$f.txt r.txt >/dev/null || { echo "z$f changed"; exit 1; } ;;
  esac
done
## and split takes them all apart
mkdir sp && ( cd sp && "$YAME" split -s ../allnames.txt ../all.cg p >/dev/null 2>&1 )
[ "$(/bin/ls sp | wc -l)" -eq 7 ] || { echo "split of the mixed store made $(/bin/ls sp | wc -l) files"; exit 1; }

## ---- the transforms, on the formats they accept ---------------------------
for f in $FMTS; do
  check "binarize f$f" -1 "$YAME" binarize f$f.cg
  check "mask f$f"     -1 "$YAME" mask f$f.cg msk.cg
  check "mask -v f$f"  -1 "$YAME" mask -v f$f.cg msk.cg
  check "dsample f$f"  -1 "$YAME" dsample -s 5 -N 8 f$f.cg
  check "perturb f$f"  -1 "$YAME" perturb -s 5 -p 0.2 f$f.cg
done

## ---- summary and hprint over the formats ----------------------------------
for f in $FMTS; do
  "$YAME" summary -m msk.cg f$f.cg >/dev/null 2>sum.err || {
    [ -s sum.err ] || { echo "summary on f$f exited non-zero silently"; exit 1; }
  }
  "$YAME" hprint -c f$f.cg >/dev/null 2>hp.err || {
    [ -s hp.err ] || { echo "hprint on f$f exited non-zero silently"; exit 1; }
  }
  "$YAME" hprint -c -g f$f.cg >/dev/null 2>&1 || true
done

## ---- rowop's inputs -------------------------------------------------------
for f in $FMTS; do
  cat f$f.cg f$f.cg > pair.cg
  for op in binasum musum stat binstring; do
    "$YAME" rowop -o $op pair.cg >/dev/null 2>op.err || {
      [ -s op.err ] || { echo "rowop $op on f$f exited non-zero silently"; exit 1; }
    }
  done
done
