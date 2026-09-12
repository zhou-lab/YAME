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
"$root/tools/make_registry.sh" --tool=yame -o "$d/yame.h" >/dev/null 2>&1 ||
  { echo "--tool=yame failed"; exit 1; }
diff -q "$root/src/registry.h" "$d/yame.h" >/dev/null ||
  { echo "regenerating --tool=yame does not reproduce the committed src/registry.h"
    diff "$root/src/registry.h" "$d/yame.h" | head -6; exit 1; }

## ---- 2. the methscope projection: same shape, its sources only --------------
"$root/tools/make_registry.sh" --tool=methscope -o "$d/ms.h" >/dev/null 2>&1 ||
  { echo "--tool=methscope failed"; exit 1; }
head -1 "$d/ms.h" | grep -q -- '--tool=methscope' ||
  { echo "the header does not carry its own regenerate line"; head -1 "$d/ms.h"; exit 1; }
for want in '"methscope", "hg38/data"' '"methscope", "hg38/models"' '"methscope", "mm10/models"'; do
  grep -q "$want" "$d/ms.h" || { echo "methscope registry lacks $want"; exit 1; }
done
for absent in InfiniumAnnotation KYCGKB '"genomes"' YAME_REF_ROWS; do
  grep -q "$absent" "$d/ms.h" && { echo "methscope registry carries $absent, which it does not consume"; exit 1; }
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
want=$("$YAME" fetch -l </dev/null 2>/dev/null | /usr/bin/awk -F'\t' '$2 == "methscope"' | wc -l)
[ "$n" -eq "$want" ] || { echo "methscope fetch -l lists $n files; yame lists $want for that source"; exit 1; }
tail -n +2 "$d/l.txt" | cut -f2 | sort -u | grep -qvx methscope && { echo "methscope fetch -l lists another source"; exit 1; }
help=$("$d/ms" -h </dev/null 2>&1) || true          # -h exits 1 by convention
printf '%s\n' "$help" | grep -q '^  methscope fetch ' || { echo "the usage is not in methscope's voice"; exit 1; }
