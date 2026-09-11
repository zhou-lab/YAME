#!/bin/bash
## Layer 5: the commands a reader runs first, on the real fixtures the docs
## fetch, through the registry. Skips cleanly when the store is not set --
## every other test is self-contained, and this one must not turn a laptop
## without the shared store into a failing suite.
##
## Meant for release tags rather than every push: it needs the network or a
## populated YAME_DATA_HOME, and the fixtures it wants are 100 KB and 27.7 MB.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
if [ -z "${YAME_DATA_HOME:-}" ]; then
  echo "skip: YAME_DATA_HOME unset (layer 5 needs the shared store)"; exit 0
fi
d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
cd "$d"

## ---- 1. a real single-cell store: info and summary agree with the file ----
"$YAME" fetch -y -c hg38/data/human_hg38_test.cg >/dev/null 2>&1 ||
  { echo "skip: fetch failed (offline?)"; exit 0; }
[ -s human_hg38_test.cg ] || { echo "fetch reported success but wrote nothing"; exit 1; }

## one record, FORMAT 6 (set + universe: the sparse query is binary calls,
## not M/U), over the genome's row count. This assertion was first written
## for format 3 and failed -- which is the point of pinning it.
"$YAME" info human_hg38_test.cg > info.txt 2>/dev/null
awk -F'\t' 'NR==2 { if ($2 != "1" || $5 != "6" || $4 != "29401795") exit 1 }' info.txt ||
  { echo "info: expected 1 record, fmt6, 29,401,795 rows; got:"; cat info.txt; exit 1; }

## summary on a fmt6 store: N_univ is the universe (sites observed), N_query
## the set (sites called 1). Both must match an independent recount from
## unpack -f -1, which prints value<tab>universe per row.
"$YAME" summary human_hg38_test.cg > sum.txt 2>/dev/null
n_univ=$(awk -F'\t' 'NR==2 {print $5}' sum.txt)
n_query=$(awk -F'\t' 'NR==2 {print $6}' sum.txt)
read -r r_univ r_query < <("$YAME" unpack -f -1 human_hg38_test.cg 2>/dev/null |
  awk -F'\t' '$2 == 1 { u++; if ($1 == 1) q++ } END { print u+0, q+0 }')
[ "$n_univ" = "$r_univ" ] && [ "$n_query" = "$r_query" ] ||
  { echo "summary N_univ=$n_univ N_query=$n_query; unpack recounts $r_univ $r_query"; exit 1; }

## ---- 2. a bundle: info must refuse, NOT read its prefix and exit 0 -------
## A .updecx is a CX record followed by a non-CX container. Before v1.40 the
## reader stopped silently at the prefix; since then it is a hard error, and
## that is the documented contract (the bound is opt-in, via cx_read_record).
## A bundle reading as a valid one-record store with exit 0 would mean the
## reader has gone lenient again -- the regression this test exists for.
"$YAME" fetch -y -c hg38/models/hg38_10k1.updecx >/dev/null 2>&1 ||
  { echo "skip: model fetch failed (offline?)"; exit 0; }
if "$YAME" info hg38_10k1.updecx >/dev/null 2>&1; then
  echo "info on a bundle exited 0: the reader accepted a non-CX tail"; exit 1
fi

## and the prefix IS a valid mask, so slicing it must work -- the recipe
## methscope's lite-decoder build uses. Section 0 of hg38_10k1.updecx is
## 8,331,342 bytes (the whole-genome model's is 14,746,955; this test was
## first written with that one and failed -- the constant is per bundle).
head -c 8331342 hg38_10k1.updecx > prefix.cm
"$YAME" info prefix.cm > pinfo.txt 2>/dev/null
awk -F'\t' 'NR==2 { if ($5 != "2" || $4 != "29401795") exit 1 }' pinfo.txt ||
  { echo "the bundle prefix did not read as a fmt2 mask over the genome:"; cat pinfo.txt; exit 1; }
echo "ok: store fixtures read as documented; bundle refused; prefix readable"
