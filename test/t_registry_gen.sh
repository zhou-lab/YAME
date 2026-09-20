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
## line 1: the table is found beside lib.sh, not beside the caller
printf '%s\n' "$files"
## line 2: a glob selects rows, field reads one by name, url_of applies the
## host rule -- including the hf: source, whose colon is why keys split on
## the LAST one
r=$(files_of 'hg38/cpg_nocontig.cr' | head -1)
printf '%s\t%s\t%s\n' "$(field "$r" store_path)" "$(field "$r" kind)" "$(url_of "$(field "$r" key)")"
## lines 3-5: `*` stops at a slash and `**` crosses it, so a directory glob
## is that directory's own files -- the same set `yame fetch <dir>` means
printf '%s\n' "$(files_of 'EPICv2/*' | wc -l)"
printf '%s\n' "$(files_of 'EPICv2/**' | wc -l)"
printf '%s\n' "$(files_of '*/KYCG/*' | wc -l)"
## lines 6-8: the key parser, on the source that carries a colon
k='hf:zhou-lab/methscope@v11:hg38_wg.updecx'
printf '%s\n' "$(key_source "$k")" "$(key_tag "$k")" "$(key_path "$k")"
printf '%s\n' "$(url_of "$k")"
EOF
( cd "$d" && bash "$d/consumer.sh" "$root" ) > "$d/lib.out" 2>"$d/lib.err" ||
  { echo "lib.sh is not sourceable from another directory"; cat "$d/lib.err"; exit 1; }

[ "$(sed -n '1p' "$d/lib.out")" = "$root/tools/registry/files.tsv" ] ||
  { echo "lib.sh resolved the table against the caller, not itself"; sed -n '1p' "$d/lib.out"; exit 1; }
[ "$(sed -n '2p' "$d/lib.out")" = "$(printf 'hg38/cpg_nocontig.cr\tgenome\thttps://raw.githubusercontent.com/zhou-lab/genomes/v4/hg38/cpg_nocontig.cr')" ] ||
  { echo "files_of / field / url_of did not agree on the coordinate stream"; sed -n '2p' "$d/lib.out"; exit 1; }
[ "$(sed -n '3p' "$d/lib.out")" -eq 6 ] ||   # -eq: BSD wc pads its count
  { echo "'EPICv2/*' should be that directory's 6 own files, got $(sed -n '3p' "$d/lib.out")"; exit 1; }
[ "$(sed -n '4p' "$d/lib.out")" -gt 20 ] ||
  { echo "'EPICv2/**' should include EPICv2/KYCG/, got $(sed -n '4p' "$d/lib.out")"; exit 1; }
[ "$(sed -n '5p' "$d/lib.out")" -gt 100 ] ||
  { echo "'*/KYCG/*' selected only $(sed -n '5p' "$d/lib.out") rows"; exit 1; }
[ "$(sed -n '6,8p' "$d/lib.out" | paste -sd' ' -)" = "hf:zhou-lab/methscope v11 hg38_wg.updecx" ] ||
  { echo "key_source/key_tag/key_path mis-split an hf: key"; sed -n '6,8p' "$d/lib.out"; exit 1; }
[ "$(sed -n '9p' "$d/lib.out")" = "https://huggingface.co/zhou-lab/methscope/resolve/v11/hg38_wg.updecx" ] ||
  { echo "url_of mishandled an hf: key"; sed -n '9p' "$d/lib.out"; exit 1; }

## ---- 5. this repo's own registry still says what it should -----------------
## The coordinate stream is published at <genome>/ by zhou-lab/genomes, and the
## knowledge base directory no longer offers a second copy of it.
grep -q '"zhou-lab/genomes@v[0-9]*:hg38/cpg_nocontig.cr", "hg38/cpg_nocontig.cr"' "$root/src/registry.h" ||
  { echo "the coordinate stream is not published at hg38/ from genomes"; exit 1; }
grep -q '"hg38/KYCG/cpg_nocontig.cr"' "$root/src/registry.h" &&
  { echo "the knowledge base still lists the coordinate stream"; exit 1; }
## and the registry is exactly the table: one C row per data row
n_tsv=$(grep -vc '^#' "$root/tools/registry/files.tsv")
n_h=$(grep -c '^    { "' "$root/src/registry.h")
n_ref=$(grep -cE '^    \{ "[^"]*", "(genome|array)", ' "$root/src/registry.h")
[ "$((n_h - n_ref))" = "$n_tsv" ] ||
  { echo "registry.h has $((n_h - n_ref)) file rows; files.tsv has $n_tsv"; exit 1; }

exit 0
