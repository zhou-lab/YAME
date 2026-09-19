#!/bin/bash
## The registry generator: yame's own projection regenerates byte-identically
## from the cached manifests (no network), and the methscope projection is the
## same shape over only the sources methscope consumes -- and compiles into a
## working fetch. This is what lets a downstream tool ship fetch over the YAME
## code it bundles, with a registry that pins what its own binary pins.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/.." && pwd)
[ -x "$root/tools/make_registry.sh" ] || { echo "skip: no generator in this tree" >&2; exit 0; }
d=$(mktemp -d); trap 'rm -rf "$d"' EXIT

## ---- 1. yame's registry is a pure function of the catalogue -----------------
## the two generations are independent and each takes ~2 s; run them together
"$root/tools/make_registry.sh" --tool=methscope -o "$d/ms.h" >/dev/null 2>&1 &
gen_ms=$!
"$root/tools/make_registry.sh" --tool=yame -o "$d/yame.h" >/dev/null 2>&1 ||
  { echo "--tool=yame failed"; exit 1; }
diff -q "$root/src/registry.h" "$d/yame.h" >/dev/null ||
  { echo "regenerating --tool=yame does not reproduce the committed src/registry.h"
    diff "$root/src/registry.h" "$d/yame.h" | head -6; exit 1; }

## ---- 2. the methscope projection: same shape, its sources only --------------
wait $gen_ms || { echo "--tool=methscope failed"; exit 1; }
head -1 "$d/ms.h" | grep -- '--tool=methscope' >/dev/null ||
  { echo "the header does not carry its own regenerate line"; head -1 "$d/ms.h"; exit 1; }
for want in '"methscope", "hg38/data"' '"methscope", "hg38/models"' '"methscope", "mm10/models"'; do
  grep "$want" "$d/ms.h" >/dev/null || { echo "methscope registry lacks $want"; exit 1; }
done

## The COORDINATE stream comes too, and only that. methscope reads its own
## upscale output against CpG positions, so hprint -r and rowsub -R need this
## one file; it has no use for a knowledgebase. Requested 2026-09-17.
##
## It comes from zhou-lab/genomes, which publishes it as of tag v4. Before
## that the only publisher was the per-genome knowledgebase repo, so this
## registry had to carry a KYCGKB row and LIFT the file out of hg38/KYCG.
for want in '"genomes", "hg38"' '"genomes", "mm10"' 'cpg_nocontig.cr'; do
  grep "$want" "$d/ms.h" >/dev/null ||
    { echo "methscope registry lacks the coordinate stream: $want"; exit 1; }
done
## and NOTHING else from that source: no seqinfo/gaps/cytoband/gene models
for absent in seqinfo gaps cytoband genes.bed; do
  grep "$absent" "$d/ms.h" >/dev/null &&
    { echo "methscope registry carries $absent, which it does not consume"; exit 1; }
done
## nor a knowledgebase, which is the source it used to come from
n_kb=$(grep -c '\.cm"' "$d/ms.h" || true)
[ "$n_kb" -eq 0 ] ||
  { echo "methscope registry carries $n_kb knowledgebase .cm files"; exit 1; }

## It must land where yame puts it, or the two tools fetch 30 MB twice into
## different places and neither sees the other's copy. Both store it by the
## ROW, under <genome>/, which is now also the directory that publishes it --
## so neither side needs a lift and there is nothing left to keep in step.
grep '"genomes", "hg38", "[^"]*", "[^"]*", "hg38", "hg38"' "$d/ms.h" >/dev/null ||
  { echo "methscope does not store the coordinate stream under hg38/"
    grep '"genomes", "hg38"' "$d/ms.h"; exit 1; }
grep 'cpg_nocontig.cr", "[0-9a-f]*", [0-9]*, NULL }' "$root/src/registry.h" >/dev/null ||
  { echo "yame no longer publishes the coordinate stream at hg38/"; exit 1; }
## and the knowledgebase no longer offers a second copy of it
grep 'YAME_FILES_KYCGKB_hg38' -A 40 "$root/src/registry.h" |
  sed -n '/^};/q;p' | grep 'cpg_nocontig.cr' >/dev/null &&
  { echo "the knowledgebase still lists the coordinate stream"; exit 1; }

for absent in InfiniumAnnotation KYCGKB YAME_REF_ROWS; do
  grep "$absent" "$d/ms.h" >/dev/null && { echo "methscope registry carries $absent, which it does not consume"; exit 1; }
done
## prior anchors come along: that is what lets a store at v8 upgrade without -f
grep -q 'YAME_PRIOR_methscope_models' "$d/ms.h" || { echo "methscope registry has no prior anchors"; exit 1; }
## and the pin is the one TAGS names
tag=$(grep -o 'methscope_models *v[0-9]*' "$root/tools/registry/TAGS" | /usr/bin/awk '{print $2}')
grep -q "\"hg38/models\", \"[^\"]*\", \"$tag\"" "$d/ms.h" ||
  { echo "methscope registry does not pin hg38/models at $tag"; exit 1; }

## ---- 3. it compiles into a fetch that lists exactly its catalogue ----------
[ -f "$root/libyame.a" ] && [ -x "$root/yame-config" ] || { echo "skip: no libyame.a for the link check" >&2; exit 0; }
cat > "$d/ms.c" <<'EOF'
#include "assets.h"
#include "ms.h"
int main(int argc, char **argv) {
  yame_fetch_cfg_t cfg = { YAME_ASSETS, YAME_ASSETS_N, "methscope", "METHSCOPE_DATA_HOME" };
  return yame_fetch_main(&cfg, argc, argv);
}
EOF
${CC:-cc} -O1 -std=gnu99 -I"$d" $("$root/yame-config" --cflags) -o "$d/ms" "$d/ms.c" \
  $("$root/yame-config" --libs) 2>"$d/cc.err" || { echo "the methscope registry does not compile"; cat "$d/cc.err"; exit 1; }
export METHSCOPE_DATA_HOME="$d/store"; mkdir -p "$METHSCOPE_DATA_HOME"
"$d/ms" -l </dev/null > "$d/l.txt" 2>/dev/null
n=$(tail -n +2 "$d/l.txt" | wc -l)
## methscope's catalogue is its own source PLUS the three coordinate streams,
## and nothing else. Counting only the methscope source would miss the
## coordinates; counting everything would let a knowledgebase creep in.
want=$("$YAME" fetch -l </dev/null 2>/dev/null |
       /usr/bin/awk -F'\t' '$2 == "methscope" || $5 == "cpg_nocontig.cr"' | wc -l)
[ "$n" -eq "$want" ] || { echo "methscope fetch -l lists $n files; expected $want"
                          tail -n +2 "$d/l.txt" | cut -f1,5 | sort | head -30; exit 1; }
## the coordinate streams are there, and land in the shared place
tail -n +2 "$d/l.txt" | /usr/bin/awk -F'\t' '$5 == "cpg_nocontig.cr" {print $4}' |
  sort -u > "$d/coord_paths.txt"
"$YAME" fetch -l </dev/null 2>/dev/null |
  /usr/bin/awk -F'\t' '$5 == "cpg_nocontig.cr" {print $4}' | sort -u > "$d/coord_yame.txt"
cmp -s "$d/coord_paths.txt" "$d/coord_yame.txt" ||
  { echo "the coordinate stream lands in different store paths for the two tools"
    paste "$d/coord_paths.txt" "$d/coord_yame.txt"; exit 1; }
## Two sources and no more: its own, and genomes for the coordinate stream.
tail -n +2 "$d/l.txt" | cut -f2 | sort -u > "$d/srcs.txt"
printf 'genomes\nmethscope\n' > "$d/srcs.want"
cmp -s "$d/srcs.txt" "$d/srcs.want" ||
  { echo "methscope fetch -l lists sources it should not:"; cat "$d/srcs.txt"; exit 1; }
## and genomes contributes ONLY the coordinate stream
tail -n +2 "$d/l.txt" | /usr/bin/awk -F'\t' '$2 == "genomes" && $5 != "cpg_nocontig.cr"' |
  head -3 > "$d/extra.txt"
[ ! -s "$d/extra.txt" ] ||
  { echo "genomes brought more than the coordinate stream:"; cut -f5 "$d/extra.txt"; exit 1; }
help=$("$d/ms" -h </dev/null 2>&1) || true          # -h exits 1 by convention
printf '%s\n' "$help" | grep '^  methscope fetch ' >/dev/null || { echo "the usage is not in methscope's voice"; exit 1; }
