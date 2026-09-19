#!/bin/bash
## The registry generator: yame's own registry is a pure function of the cached
## catalogue, regenerating byte-identically with no network; --check is the CI
## form of that; and the generator writes ONE file and cannot be aimed anywhere
## else.
##
## It used to also drive --tool=kycg|sesame|methscope and assert on those
## projections. Each of those tools now owns its emitter and sources
## tools/registry/lib.sh, so their registries are their own repos' tests. What
## stays here is what this repo still owns: its own registry, the lookups it
## publishes, and the boundary around where output may land.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/.." && pwd)
gen=$root/tools/make_registry.sh
[ -x "$gen" ] || { echo "skip: no generator in this tree" >&2; exit 0; }
d=$(mktemp -d); trap 'rm -rf "$d"' EXIT

## ---- 1. the committed registry matches the catalogue -----------------------
"$gen" --check >/dev/null 2>&1 ||
  { echo "src/registry.h does not match the catalogue"; "$gen" --check 2>&1 | head -8; exit 1; }

## ---- 2. --check actually compares, rather than always passing --------------
cp "$root/src/registry.h" "$d/keep.h"
printf '/* drift */\n' >> "$root/src/registry.h"
if "$gen" --check >/dev/null 2>&1; then
  cp "$d/keep.h" "$root/src/registry.h"
  echo "--check passed a registry that does not match the catalogue"; exit 1
fi
cp "$d/keep.h" "$root/src/registry.h"
cmp -s "$d/keep.h" "$root/src/registry.h" ||
  { echo "the drift probe did not restore src/registry.h"; exit 1; }

## ---- 3. there is one destination, and no way to name another ---------------
## A generator that could be pointed at a path is how this repo came to write
## other tools' headers. Both spellings that allowed it must now fail, and
## neither may leave a file behind.
for bad in "--tool=kycg" "-o $d/elsewhere.h" "--out=$d/elsewhere.h"; do
  if "$gen" $bad >/dev/null 2>&1; then
    echo "the generator accepted '$bad'"; exit 1
  fi
done
[ -e "$d/elsewhere.h" ] && { echo "a rejected run still wrote a file"; exit 1; }

## ---- 4. the catalogue lookups are sourceable on their own ------------------
## This is the surface a downstream tool depends on; it must work from a
## caller in another directory, with its own $0, having sourced nothing else.
cat > "$d/consumer.sh" <<'EOF'
set -euo pipefail
. "$1/tools/registry/lib.sh"
printf '%s\t%s\t%s\n' "$(tag_of genomes)" "$(base_of genomes)" "$cat_dir"
printf '%s\n' "$(anchor_of genomes "$(tag_of genomes)" hg38)"
printf '%s\n' "$(nsets_of KYCGKB "$(tag_of KYCGKB)" hg38)"
EOF
( cd "$d" && bash "$d/consumer.sh" "$root" ) > "$d/lib.out" 2>"$d/lib.err" ||
  { echo "lib.sh is not sourceable from another directory"; cat "$d/lib.err"; exit 1; }

tag=$(cut -f1 < <(head -1 "$d/lib.out"))
[ -n "$tag" ] || { echo "tag_of returned nothing through lib.sh"; exit 1; }
grep -q "^genomes[[:space:]]\+$tag[[:space:]]" "$root/tools/registry/TAGS" ||
  { echo "tag_of disagrees with TAGS: got '$tag'"; exit 1; }
## the anchor is a sha256 computed from bytes on disk, not a placeholder
sed -n '2p' "$d/lib.out" | grep -Eq '^[0-9a-f]{64}$' ||
  { echo "anchor_of did not return a sha256"; sed -n '2p' "$d/lib.out"; exit 1; }
sed -n '3p' "$d/lib.out" | grep -Eq '^[0-9]+$' ||
  { echo "nsets_of did not return a count"; sed -n '3p' "$d/lib.out"; exit 1; }
## and it found the catalogue beside itself, not beside the caller
[ "$(cut -f3 < <(head -1 "$d/lib.out"))" = "$root/tools/registry/catalog" ] ||
  { echo "lib.sh resolved cat_dir against the caller, not itself"
    head -1 "$d/lib.out"; exit 1; }

## ---- 5. this repo's own registry still says what it should -----------------
## The coordinate stream is published at <genome>/ by the genomes unit, with no
## lift, and the knowledgebase no longer offers a second copy of it.
grep -q 'cpg_nocontig.cr", "[0-9a-f]*", [0-9]* }' "$root/src/registry.h" ||
  { echo "yame no longer publishes the coordinate stream at <genome>/"; exit 1; }
grep 'YAME_FILES_KYCGKB_hg38' -A 40 "$root/src/registry.h" |
  sed -n '/^};/q;p' | grep -q 'cpg_nocontig.cr' &&
  { echo "the knowledgebase still lists the coordinate stream"; exit 1; }

exit 0
