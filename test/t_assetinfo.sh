#!/bin/bash
## The descriptions and the catalogue must agree.
##
## data/assets.tsv is what the browser says a file IS; the registry is what
## fetch can download. They are edited separately -- a pin bump moves the
## registry, a content request moves the TSV -- and nothing tied them: v9's
## four new models showed a bare "-" for a title for three days, one row
## still described a model replaced three tags earlier, and one row named a
## file the catalogue no longer carried. Both directions are checked here.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/.." && pwd)
tsv="$root/data/assets.tsv"
[ -f "$tsv" ] || { echo "skip: no data/assets.tsv in this tree" >&2; exit 0; }
d=$(mktemp -d); trap 'rm -rf "$d"' EXIT

## ---- 1. every file the registry lists has a title ---------------------------
"$YAME" fetch -l </dev/null 2>/dev/null > "$d/l.tsv"
/usr/bin/awk -F'\t' 'NR > 1 && $10 == "" { print $1 "/" $5 }' "$d/l.tsv" > "$d/untitled.txt"
if [ -s "$d/untitled.txt" ]; then
  echo "files in the catalogue with no description in data/assets.tsv:"
  sed 's/^/  /' "$d/untitled.txt"; exit 1
fi

## ---- 2. every MODEL row in the TSV names a file the catalogue carries -------
## Keys are file stems (hg38_celltype -> hg38_celltype.clfx). A row for a
## withdrawn file describes something nobody can fetch. Only the model rows
## are held to this: annotation keys map to files by other rules.
/usr/bin/awk -F'\t' 'NR > 1 && $1 ~ /models$/ { f = $5; sub(/\.[a-z]+$/, "", f); print f }' "$d/l.tsv" | sort -u > "$d/model_files.txt"
/usr/bin/awk -F'\t' '!/^#/ && NF == 8 && $1 != "key" && ($2 ~ /^(hg38|mm10)$/) && $4 ~ /classifier|decoder|Deconvolution/ { print $1 }' "$tsv" | sort -u > "$d/model_rows.txt"
comm -23 "$d/model_rows.txt" "$d/model_files.txt" > "$d/orphans.txt"
if [ -s "$d/orphans.txt" ]; then
  echo "model rows in data/assets.tsv for files the catalogue does not carry:"
  sed 's/^/  /' "$d/orphans.txt"; exit 1
fi

## ---- 3. the compiled header is what the TSV says ---------------------------
## assetinfo.h is generated from the TSV and committed; an edit to one without
## the other is the drift this file exists to catch.
[ -x "$root/tools/make_assetinfo.sh" ] || exit 0
"$root/tools/make_assetinfo.sh" -o "$d/assetinfo.h" >/dev/null 2>&1 ||
  { echo "make_assetinfo.sh failed"; exit 1; }
diff -q "$root/src/assetinfo.h" "$d/assetinfo.h" >/dev/null ||
  { echo "src/assetinfo.h is stale: regenerate with tools/make_assetinfo.sh -o src/assetinfo.h"; exit 1; }

## ---- 4. the rows themselves are well formed -------------------------------
/usr/bin/awk -F'\t' '!/^#/ && NF && NF != 8 { print "row " NR " has " NF " fields: " $1; bad = 1 } END { exit bad }' "$tsv"
