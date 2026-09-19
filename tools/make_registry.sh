#!/usr/bin/env bash
## Generate a tool's compiled asset registry from the shared catalog.
##
##   tools/make_registry.sh --tool=yame   [-o src/registry.h]
##   tools/make_registry.sh --tool=kycg   [-o src/registry.h]
##   tools/make_registry.sh --tool=methscope [-o src/registry.h]   yame shape, methscope's sources
##   tools/make_registry.sh --tool=sesame [-o src/registry.h]
##   tools/make_registry.sh --refresh [--tag=vN]
##
## kycg and sesame-cli each grew their own generator, pinned their own tags, and
## so could not share a download. This is the one generator; tools/registry/ is
## the one catalog. A tool's registry.h is a projection of it.
##
## EMISSION IS OFFLINE. Every anchor is sha256 of a manifest cached under
## tools/registry/sums/<source>/<tag>/, hashed from the bytes on disk, so a pin
## can never disagree with the content it claims to describe and CI regenerates
## byte-identically with no network. --refresh is the only mode that downloads:
## it re-caches those manifests (and the KYCGKB file sizes, the one fact no
## manifest carries) so a tag bump is a reviewable diff of cached bytes rather
## than a silent change of digests.
set -euo pipefail

## The catalog lookups live beside the catalog, in a file any tool can source.
## Sourcing sets reg, cat_dir and sums_dir as well.
here=$(cd "$(dirname "$0")" && pwd)
. "$here/registry/lib.sh"

tool=""
out=""
refresh=0
new_tag=""

for arg in "$@"; do
  case $arg in
    --tool=*)  tool=${arg#--tool=} ;;
    --tag=*)   new_tag=${arg#--tag=} ;;
    --refresh) refresh=1 ;;
    -o)        die_next_is_o=1 ;;   ## handled below
    *)         ;;
  esac
done

## -o takes a value, which the loop above cannot see; parse it separately.
prev=""
for arg in "$@"; do
  if [ "$prev" = "-o" ]; then out=$arg; fi
  prev=$arg
done

## The tags this build supersedes for one directory: every cached tag STRICTLY
## EARLIER than the pinned one, as "<anchor>\t<tag>".
##
## "Earlier" is enforced here rather than assumed. The rule used to be every
## cached tag except the current one, which is correct only while the procedure
## is followed -- cache a manifest and bump TAGS in the same commit. Pre-caching
## a future tag's manifest (to make a later bump a one-line change, which is a
## reasonable thing to want) then listed it as an ancestor, and
## YAME_PIN_ANCESTOR grants overwrite-without--f on the promise that this build
## holds the LATER tag. A store already at the newer tag would be walked back to
## the pinned one and the message would call it catching up. Comparing versions
## makes that structurally impossible instead of merely discouraged, so a cached
## manifest is now safe to land ahead of its bump.
priors_of() {   ## source tag subpath
  local d t p
  for d in "$sums_dir/$1"/*/; do
    [ -d "$d" ] || continue
    t=$(basename "$d")
    [ "$t" = "$2" ] && continue
    ## sort -V puts the pair in version order; keep t only when it sorts first.
    [ "$(printf '%s\n%s\n' "$t" "$2" | sort -V | head -1)" = "$t" ] || continue
    p=$(sums_path "$1" "$t" "$3")
    [ -s "$p" ] || continue
    printf '%s\t%s\n' "$(sha256_of "$p")" "$t"
  done
}

## The ancestry table for one directory, emitted only when there is one.
emit_prior_table() {   ## slug source tag subpath
  local rows; rows=$(priors_of "$2" "$3" "$4")
  [ -n "$rows" ] || return 0
  printf 'static const yame_pin_prior_t YAME_PRIOR_%s[] = {\n' "$1"
  printf '%s\n' "$rows" | while IFS=$'\t' read -r sha t; do
    printf '    { "%s", "%s" },\n' "$t" "$sha"
  done
  printf '};\n\n'
}

## The two struct fields naming that table, or the empty ancestry.
prior_ref() {   ## slug source tag subpath
  if [ -n "$(priors_of "$2" "$3" "$4")" ]; then
    printf 'YAME_PRIOR_%s, YAME_NPRIOR(YAME_PRIOR_%s)' "$1" "$1"
  else
    printf 'NULL, 0'
  fi
}

# --------------------------------------------------------------- refresh mode

if [ "$refresh" = 1 ]; then
  ia_tag=${new_tag:-$(tag_of InfiniumAnnotation)}
  ia_base=$(base_of InfiniumAnnotation)
  echo "refreshing InfiniumAnnotation @ $ia_tag" >&2
  for p in $(rows_of "$cat_dir/InfiniumAnnotation.tsv" | cut -f1); do
    for sub in "$p" "$p/KYCG"; do
      d=$sums_dir/InfiniumAnnotation/$ia_tag/$sub
      mkdir -p "$d"
      curl -sfL -o "$d/SHA256SUMS" "$ia_base/$ia_tag/$sub/SHA256SUMS" \
        || { echo "  MISS $sub" >&2; rm -f "$d/SHA256SUMS"; }
    done
  done

  g_tag=$(tag_of genomes); g_base=$(base_of genomes)
  echo "refreshing genomes @ $g_tag" >&2
  for g in $(rows_of "$cat_dir/genomes.tsv" | cut -f1); do
    d=$sums_dir/genomes/$g_tag/$g; mkdir -p "$d"
    curl -sfL -o "$d/SHA256SUMS" "$g_base/$g_tag/$g/SHA256SUMS" \
      || { echo "  MISS $g" >&2; rm -f "$d/SHA256SUMS"; }
  done

  kb_tag=$(tag_of KYCGKB); kb_base=$(base_of KYCGKB)
  echo "refreshing KYCGKB @ $kb_tag" >&2
  rows_of "$cat_dir/KYCGKB.tsv" | while IFS=$'\t' read -r g _tools repo _rest; do
    d=$sums_dir/KYCGKB/$kb_tag/$g; mkdir -p "$d"
    curl -sfL -o "$d/SHA256SUMS" "$kb_base/$repo/raw/$kb_tag/SHA256SUMS" \
      || { echo "  MISS $g" >&2; rm -f "$d/SHA256SUMS"; }
  done

  ms_tag=$(tag_of methscope); ms_base=$(base_of methscope)
  echo "refreshing methscope @ $ms_tag" >&2
  rows_of "$cat_dir/methscope.tsv" | while IFS=$'\t' read -r _g sub; do
    d=$sums_dir/methscope/$ms_tag/$sub; mkdir -p "$d"
    curl -sfL -o "$d/SHA256SUMS" "$ms_base/$ms_tag/$sub/SHA256SUMS"
  done

  msm_tag=$(tag_of methscope_models); msm_base=$(base_of methscope_models)
  echo "refreshing methscope_models @ $msm_tag" >&2
  ## One manifest at the repo root serves every row, so fetch it ONCE rather
  ## than per row, and fail loudly. This used to be a bare curl inside the row
  ## loop with no status check: HuggingFace answered 429 during the v10 bump,
  ## the line above still said "refreshing", and the failure only surfaced
  ## later as "No such file or directory" from the emit step.
  msm_d=$sums_dir/methscope_models/$msm_tag
  mkdir -p "$msm_d"
  msm_url="$msm_base/$msm_tag/SHA256SUMS"
  msm_ok=0
  for try in 1 2 3; do
    code=$(curl -sL -o "$msm_d/SHA256SUMS.part" -w '%{http_code}' "$msm_url" || echo 000)
    if [ "$code" = 200 ]; then
      mv "$msm_d/SHA256SUMS.part" "$msm_d/SHA256SUMS"; msm_ok=1; break
    fi
    rm -f "$msm_d/SHA256SUMS.part"
    ## 429 and 503 are "come back later" and are worth a wait; anything else
    ## will not improve by asking again.
    case $code in
      429|503) echo "  HTTP $code from $msm_url; waiting 45 s (try $try of 3)" >&2
               sleep 45 ;;
      *)       break ;;
    esac
  done
  if [ "$msm_ok" != 1 ]; then
    echo "make_registry.sh: cannot fetch $msm_url (HTTP $code)." >&2
    echo "  The pin in TAGS is $msm_tag. Check the tag exists, then re-run;" >&2
    echo "  a 429 is HuggingFace rate-limiting and clears on its own." >&2
    exit 1
  fi

  ## File sizes: the GitHub contents API is the only source, so this is the one
  ## part of the catalog that cannot be rebuilt from a manifest. Every source,
  ## because a browser that shows a size for some files and not others reads as
  ## if the others were free.
  echo "refreshing file sizes" >&2

  ## The API rate limit is 60/hour unauthenticated and this walks ~20
  ## directories, so use gh (5000/hour with a token) when it is there.
  api() {
    if command -v gh >/dev/null 2>&1; then gh api "$1" 2>/dev/null
    else curl -sfL "https://api.github.com/$1"; fi
  }

  sizes_of() {                     ## scope repo path ref
    api "repos/zhou-lab/$2/contents/$3?ref=$4" | python3 -c "
import json, sys
try:
    d = json.load(sys.stdin)
except Exception:
    sys.exit(0)
if not isinstance(d, list): sys.exit(0)
for f in sorted(d, key=lambda x: x['name']):
    if f.get('type') != 'file': continue
    if f['name'] in ('README.md', 'SHA256SUMS'): continue
    print('$1\t%s\t%d' % (f['name'], f['size']))
"
  }

  tmp=$(mktemp)
  {
    echo "# Published file sizes, for display only -- nothing depends on them,"
    echo "# so a file added upstream needs no rebuild. Refreshed from the GitHub"
    echo "# contents API (the one fact in this catalog that no manifest carries);"
    echo "# ASCII-sorted by name within a scope, the order the emitted tables"
    echo "# must keep."
    echo "#"
    echo "# scope is <source>/<subpath>: the directory the file lives in."
    echo "#"
    printf "# scope\tname\tsize\n"

    rows_of "$cat_dir/InfiniumAnnotation.tsv" | cut -f1 | while read -r p; do
      sizes_of "InfiniumAnnotation/$p"      InfiniumAnnotation "$p"      "$ia_tag"
      sizes_of "InfiniumAnnotation/$p/KYCG" InfiniumAnnotation "$p/KYCG" "$ia_tag"
    done
    rows_of "$cat_dir/KYCGKB.tsv" | while IFS=$'\t' read -r g _tools repo _rest; do
      sizes_of "KYCGKB/$g" "$repo" "" "$kb_tag"
    done
    rows_of "$cat_dir/genomes.tsv" | cut -f1 | while read -r g; do
      sizes_of "genomes/$g" genomes "$g" "$g_tag"
    done
    rows_of "$cat_dir/methscope.tsv" | while IFS=$'\t' read -r g sub; do
      sizes_of "methscope/$g" methscope_data "$sub" "$ms_tag"
    done
  } > "$tmp"
  mv "$tmp" "$cat_dir/file_sizes.tsv"
  echo "done; review the diff under tools/registry/ before committing" >&2
  exit 0
fi

[ -n "$tool" ] || { sed -n '2,9p' "$0" | sed 's/^## \{0,1\}//' >&2; exit 1; }


## Per-directory file list, straight out of the cached manifest: the names and
## digests are already there, so compiling them in costs nothing and lets a
## browser offer individual files instead of whole directories. Sizes come from
## catalog/file_sizes.tsv (the contents API); 0 means "not published", and
## nothing depends on it.
emit_file_table() {   ## slug source tag subpath scope [prefix] [name:slot] [skip]
  local slug=$1 src=$2 tg=$3 sub=$4 scope=${5:-} pfx=${6:-} lift=${7:-} skip=${8:-}
  local p; p=$(sums_path "$src" "$tg" "$sub")
  printf 'static const yame_asset_file_t YAME_FILES_%s[] = {\n' "$slug"
  while read -r sha name; do
    [ -n "${name:-}" ] || continue
    ## One upstream directory can back several entries; take only this one's.
    case "$name" in "$pfx"*) ;; *) continue ;; esac
    ## A file this directory publishes but another source owns. It stays in
    ## the upstream manifest -- we do not control that repo -- but listing it
    ## here too would fetch the same bytes twice into the same store slot.
    [ -n "$skip" ] && [ "$name" = "$skip" ] && continue
    local size=0
    if [ -n "$scope" ]; then
      size=$(rows_of "$cat_dir/file_sizes.tsv" |
             /usr/bin/awk -F'\t' -v g="$scope" -v n="$name" \
                          '$1==g && $2==n {print $3}')
      [ -n "$size" ] || size=0
    fi
    ## A file that belongs elsewhere in the store than the directory
    ## publishing it, given as name:slot.
    local slot=NULL
    if [ -n "$lift" ] && [ "$name" = "${lift%%:*}" ]; then
      slot="\"${lift#*:}\""
    fi
    printf '    { "%s", "%s", %s, %s },\n' "$name" "$sha" "$size" "$slot"
  done < "$p"
  printf '    { NULL, NULL, 0, NULL }\n};\n\n'
}

## A C identifier for "<source>/<target>".
slug_of() { echo "$1/$2" | tr -c 'A-Za-z0-9\n' '_'; }

# ------------------------------------------------------------------ emitters

## Emit a registry in the yame shape (yame_asset_reg_t rows) for TOOL, over
## the SOURCES named -- every source by default, which is yame's own registry;
## a downstream tool that ships fetch over the YAME code it bundles asks for
## only the sources it consumes. One function, so the two cannot drift.
emit_yame() {
  local tool=${1:-yame}
  local sources=${2:-"InfiniumAnnotation genomes KYCGKB methscope methscope_models"}
  wants() { case " $sources " in *" $1 "*) return 0 ;; esac; return 1; }
  local ia_tag ia_base g_tag g_base kb_tag kb_base ms_tag ms_base msm_tag msm_base
  ia_tag=$(tag_of InfiniumAnnotation); ia_base=$(base_of InfiniumAnnotation)
  g_tag=$(tag_of genomes);             g_base=$(base_of genomes)
  kb_tag=$(tag_of KYCGKB);             kb_base=$(base_of KYCGKB)
  ms_tag=$(tag_of methscope);          ms_base=$(base_of methscope)
  msm_tag=$(tag_of methscope_models);  msm_base=$(base_of methscope_models)

  cat <<EOF
/* registry.h -- GENERATED by tools/make_registry.sh --tool=$tool. Do not edit.
 *
 * What \`$tool fetch\` knows how to download, and the digest each directory's
 * manifest must have. One row per fetchable directory.
 *
 * store_sub is the browser path: what the tree shows under a name is what
 * lands under that name in the store, so \`hg38/data\` in the browser is
 * <root>/hg38/data on disk. Which upstream repo filled a directory is in
 * \`source\`, not in the path.
 *
 * The anchor is sha256(SHA256SUMS) at the pinned tag, hashed from the copy
 * cached under tools/registry/sums/. Verifying a fetched manifest against it
 * is what makes every per-file digest trustworthy, and comparing a STORED
 * manifest against it is how a directory says which tag filled it.
 *
 * \`prior\` is the same hash for every EARLIER tag still cached in the repo.
 * Without it a store one tag behind is indistinguishable from a store some
 * other tool filled, and both need -f; with it, catching up is just a fetch.
 *
 * SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2021-present Wanding Zhou
 */
#ifndef _YAME_REGISTRY_H
#define _YAME_REGISTRY_H

#include <stddef.h>
#include <stdint.h>

#include "assets.h"   /* yame_pin_prior_t, yame_asset_file_t, yame_asset_reg_t */

/* yame_asset_file_t and yame_asset_reg_t -- the row types this table is made
 * of -- are declared in assets.h, not here: library code takes ANY tool's
 * registry as an argument, so the types cannot live in one tool's generated
 * file. This header is only the data. */

EOF

  ## The file tables first: the asset rows point at them.
  if wants InfiniumAnnotation; then
  rows_of "$cat_dir/InfiniumAnnotation.tsv" | while IFS=$'\t' read -r p _t _b _r _o; do
    emit_file_table "$(slug_of InfiniumAnnotation "$p")" InfiniumAnnotation \
                    "$ia_tag" "$p" "InfiniumAnnotation/$p"
    emit_file_table "$(slug_of InfiniumAnnotation "$p/KYCG")" InfiniumAnnotation \
                    "$ia_tag" "$p/KYCG" "InfiniumAnnotation/$p/KYCG"
  done
  fi
  if wants KYCGKB; then
  rows_of "$cat_dir/KYCGKB.tsv" | while IFS=$'\t' read -r g _tools repo _rest; do
    emit_file_table "$(slug_of KYCGKB "$g")" KYCGKB "$kb_tag" "$g" "KYCGKB/$g" \
                    "" "" "cpg_nocontig.cr"
  done
  fi
  if wants genomes; then
  rows_of "$cat_dir/genomes.tsv" | while IFS=$'\t' read -r g _tools; do
    emit_file_table "$(slug_of genomes "$g")" genomes "$g_tag" "$g" "genomes/$g"
  done
  fi
  ## The coordinate stream ALONE, for a tool that needs to read its own output
  ## against CpG positions but has no use for a knowledgebase. Requested by
  ## methscope, 2026-09-17.
  ##
  ## It comes from zhou-lab/genomes, which is where the file is published as of
  ## genomes v4. Before that it could only be had from the KYCGKB repo, so a
  ## coordinates fetch dragged a knowledgebase's manifest along and the file
  ## had to be LIFTED out of <genome>/KYCG/ into <genome>/. Now the directory
  ## that publishes it is the directory it belongs in, and the same bytes serve
  ## the genomes unit and this one.
  if wants coordinates; then
  rows_of "$cat_dir/genomes.tsv" | while IFS=$'\t' read -r g _tools; do
    emit_file_table "coordinates_$g" genomes "$g_tag" "$g" "genomes/$g" \
                    "cpg_nocontig.cr"
  done
  fi
  if wants methscope; then
  rows_of "$cat_dir/methscope.tsv" | while IFS=$'\t' read -r g sub; do
    emit_file_table "$(slug_of methscope "$g")" methscope "$ms_tag" "$sub" "methscope/$g"
  done
  fi
  if wants methscope_models; then
  rows_of "$cat_dir/methscope_models.tsv" | while IFS=$'\t' read -r t pfx store; do
    emit_file_table "$(slug_of methscope_models "$t")" methscope_models \
                    "$msm_tag" "" "methscope/models" "$pfx"
  done
  fi

  ## Then the ancestries, for the directories that have one. A source only
  ## grows one the first time it is bumped, so most of these are absent and
  ## their rows carry NULL -- which is exactly the old behaviour.
  if wants InfiniumAnnotation; then
  rows_of "$cat_dir/InfiniumAnnotation.tsv" | while IFS=$'\t' read -r p _t _b _r _o; do
    emit_prior_table "$(slug_of InfiniumAnnotation "$p")" InfiniumAnnotation \
                     "$ia_tag" "$p"
    emit_prior_table "$(slug_of InfiniumAnnotation "$p/KYCG")" InfiniumAnnotation \
                     "$ia_tag" "$p/KYCG"
  done
  fi
  if wants KYCGKB; then
  rows_of "$cat_dir/KYCGKB.tsv" | while IFS=$'\t' read -r g _tools repo _rest; do
    emit_prior_table "$(slug_of KYCGKB "$g")" KYCGKB "$kb_tag" "$g"
  done
  fi
  if wants genomes; then
  rows_of "$cat_dir/genomes.tsv" | while IFS=$'\t' read -r g _tools; do
    emit_prior_table "$(slug_of genomes "$g")" genomes "$g_tag" "$g"
  done
  fi
  if wants coordinates; then
  rows_of "$cat_dir/genomes.tsv" | while IFS=$'\t' read -r g _tools; do
    emit_prior_table "coordinates_$g" genomes "$g_tag" "$g"
  done
  fi
  if wants methscope; then
  rows_of "$cat_dir/methscope.tsv" | while IFS=$'\t' read -r g sub; do
    emit_prior_table "$(slug_of methscope "$g")" methscope "$ms_tag" "$sub"
  done
  fi
  if wants methscope_models; then
  rows_of "$cat_dir/methscope_models.tsv" | while IFS=$'\t' read -r t pfx store; do
    emit_prior_table "$(slug_of methscope_models "$t")" methscope_models \
                     "$msm_tag" ""
  done
  fi

  cat <<EOF
static const yame_asset_reg_t YAME_ASSETS[] = {
EOF

  ## Array platforms: the platform directory, then its KYCG/ subdirectory.
  ## Both are fetchable on their own -- sesame wants the first, kycg wants both.
  if wants InfiniumAnnotation; then
  rows_of "$cat_dir/InfiniumAnnotation.tsv" | while IFS=$'\t' read -r p _t _b _r _o; do
    s1=$(slug_of InfiniumAnnotation "$p"); s2=$(slug_of InfiniumAnnotation "$p/KYCG")
    printf '    { "InfiniumAnnotation", "%s", "%s", "%s", "%s", "%s", "%s", YAME_FILES_%s, YAME_NFILES(YAME_FILES_%s), %s },\n' \
      "$p" "$ia_base" "$ia_tag" "$p" "$p" "$(anchor_of InfiniumAnnotation "$ia_tag" "$p")" "$s1" "$s1" \
      "$(prior_ref "$s1" InfiniumAnnotation "$ia_tag" "$p")"
    printf '    { "InfiniumAnnotation", "%s/KYCG", "%s", "%s", "%s/KYCG", "%s/KYCG", "%s", YAME_FILES_%s, YAME_NFILES(YAME_FILES_%s), %s },\n' \
      "$p" "$ia_base" "$ia_tag" "$p" "$p" "$(anchor_of InfiniumAnnotation "$ia_tag" "$p/KYCG")" "$s2" "$s2" \
      "$(prior_ref "$s2" InfiniumAnnotation "$ia_tag" "$p/KYCG")"
  done
  fi

  ## The coordinate stream on its own. Target and store_sub are <genome>, the
  ## same spelling and the same place yame uses, so `methscope fetch
  ## hg38/cpg_nocontig.cr` and `yame fetch hg38/cpg_nocontig.cr` are the same
  ## command over the same file -- and the same one the genomes unit below
  ## writes, so fetching both costs one download.
  if wants coordinates; then
  rows_of "$cat_dir/genomes.tsv" | while IFS=$'\t' read -r g _tools; do
    printf '    { "genomes", "%s", "%s", "%s", "%s", "%s", "%s", YAME_FILES_coordinates_%s, YAME_NFILES(YAME_FILES_coordinates_%s), %s },\n' \
      "$g" "$g_base" "$g_tag" "$g" "$g" "$(anchor_of genomes "$g_tag" "$g")" "$g" "$g" \
      "$(prior_ref "coordinates_$g" genomes "$g_tag" "$g")"
  done
  fi

  ## Whole-genome knowledgebases: one repo each, manifest at the repo root.
  if wants KYCGKB; then
  rows_of "$cat_dir/KYCGKB.tsv" | while IFS=$'\t' read -r g _tools repo _rest; do
    s=$(slug_of KYCGKB "$g")
    printf '    { "KYCGKB", "%s", "%s/%s/raw", "%s", "", "%s/KYCG", "%s", YAME_FILES_%s, YAME_NFILES(YAME_FILES_%s), %s },\n' \
      "$g" "$kb_base" "$repo" "$kb_tag" "$g" "$(anchor_of KYCGKB "$kb_tag" "$g")" "$s" "$s" \
      "$(prior_ref "$s" KYCGKB "$kb_tag" "$g")"
  done
  fi

  ## Genome-level annotation (seqinfo / gaps / cytoband).
  if wants genomes; then
  rows_of "$cat_dir/genomes.tsv" | while IFS=$'\t' read -r g _tools; do
    s=$(slug_of genomes "$g")
    printf '    { "genomes", "%s", "%s", "%s", "%s", "%s", "%s", YAME_FILES_%s, YAME_NFILES(YAME_FILES_%s), %s },\n' \
      "$g" "$g_base" "$g_tag" "$g" "$g" "$(anchor_of genomes "$g_tag" "$g")" "$s" "$s" \
      "$(prior_ref "$s" genomes "$g_tag" "$g")"
  done
  fi

  ## Example query methylomes -- the only source here that is data rather than
  ## annotation. The target carries the genome so the browser files them under
  ## it; the manifest sits at a differently-named path upstream.
  if wants methscope; then
  rows_of "$cat_dir/methscope.tsv" | while IFS=$'\t' read -r g sub; do
    s=$(slug_of methscope "$g")
    printf '    { "methscope", "%s", "%s", "%s", "%s", "%s", "%s", YAME_FILES_%s, YAME_NFILES(YAME_FILES_%s), %s },\n' \
      "$g" "$ms_base" "$ms_tag" "$sub" "$g" "$(anchor_of methscope "$ms_tag" "$sub")" "$s" "$s" \
      "$(prior_ref "$s" methscope "$ms_tag" "$sub")"
  done
  fi

  ## Model bundles: manifest at the repo root, so remote_sub is empty. Mixed
  ## genomes in one directory, so the target is not a genome name.
  if wants methscope_models; then
  rows_of "$cat_dir/methscope_models.tsv" | while IFS=$'\t' read -r t pfx store; do
    s=$(slug_of methscope_models "$t")
    printf '    { "methscope", "%s", "%s", "%s", "", "%s", "%s", YAME_FILES_%s, YAME_NFILES(YAME_FILES_%s), %s },\n' \
      "$t" "$msm_base" "$msm_tag" "$t" "$(anchor_of methscope_models "$msm_tag" "")" "$s" "$s" \
      "$(prior_ref "$s" methscope_models "$msm_tag" "")"
  done
  fi

  cat <<'EOF'
    { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0 }
};

#define YAME_ASSETS_N (sizeof(YAME_ASSETS)/sizeof(YAME_ASSETS[0]) - 1)
EOF

  ## The row-space table serves refstore.c, which identifies a reference from
  ## a file's row count and needs the genome and knowledgebase sources; a tool
  ## without them has no use for it, and its rows would name directories the
  ## tool's registry does not carry.
  if wants genomes && wants KYCGKB; then
  cat <<'EOF'

/**
 * Row-space fingerprints: how many rows a file written against each reference
 * has.
 *
 * A CX record carries its row count and nothing about what those rows ARE, so
 * the count is the only handle on which reference it belongs to -- and it is
 * enough, because no two row spaces in the catalogue share one. That is what
 * lets a command infer -R rather than make the caller repeat what the file
 * already implies.
 *
 * `rows` is hand-pinned in the catalog: it is a property of a file's contents,
 * which no manifest carries. `fetch` is the argument that would put the
 * reference in the store, for the error message when it is not there.
 */
typedef struct yame_ref_rows_s {
    const char *name;        /* hg38, EPIC, ... */
    const char *kind;        /* "genome" (a .cr) or "array" (an ordering) */
    uint64_t    rows;
    const char *store_path;  /* the reference itself, under the store root */
    const char *const *dirs; /* every directory it owns, nearest first */
    const char *fetch;       /* the browser path to give `yame fetch` */
} yame_ref_rows_t;

EOF

  ## Every directory a row space owns, nearest first: a name is looked up in
  ## the unit's own directory before its knowledgebase, because that is the
  ## order someone naming a file means them in. Searching only the
  ## knowledgebase left half the catalogue unreachable by name.
  if wants KYCGKB; then
  rows_of "$cat_dir/KYCGKB.tsv" | while IFS=$'\t' read -r g _rest; do
    printf 'static const char *const YAME_REF_DIRS_%s[] = { "%s", "%s/KYCG", NULL };\n' \
      "$g" "$g" "$g"
  done
  fi
  rows_of "$cat_dir/InfiniumAnnotation.tsv" | cut -f1 | while read -r p; do
    printf 'static const char *const YAME_REF_DIRS_%s[] = { "%s", "%s/KYCG", NULL };\n' \
      "$p" "$p" "$p"
  done
  echo

  cat <<'EOF'
static const yame_ref_rows_t YAME_REF_ROWS[] = {
EOF

  if wants KYCGKB; then
  rows_of "$cat_dir/KYCGKB.tsv" | while IFS=$'\t' read -r g _tools _repo nrows _rest; do
    [ -n "${nrows:-}" ] || continue
    printf '    { "%s", "genome", %s, "%s/cpg_nocontig.cr", YAME_REF_DIRS_%s, "%s" },\n' \
      "$g" "$nrows" "$g" "$g" "$g"
  done
  fi
  rows_of "$cat_dir/InfiniumAnnotation.tsv" |
    while IFS=$'\t' read -r p _tools _beads nrows _rest; do
      [ -n "${nrows:-}" ] || continue
      printf '    { "%s", "array", %s, "%s/%s.ordering.tsv.gz", YAME_REF_DIRS_%s, "%s" },\n' \
        "$p" "$nrows" "$p" "$p" "$p" "$p"
    done

  cat <<'EOF'
    { NULL, NULL, 0, NULL, NULL, NULL }
};

#define YAME_REF_ROWS_N (sizeof(YAME_REF_ROWS)/sizeof(YAME_REF_ROWS[0]) - 1)
EOF
  fi

  cat <<'EOF'

#endif /* _YAME_REGISTRY_H */
EOF
}

emit_sesame() {
  local tag base g_tag g_base
  tag=$(tag_of InfiniumAnnotation); base=$(base_of InfiniumAnnotation)
  g_tag=$(tag_of genomes);          g_base=$(base_of genomes)

  cat <<EOF
/* registry.h -- GENERATED by tools/make_registry.sh. Do not edit.
 *
 * Per platform: the IDAT bead count (for auto-detection), the ordering file
 * name, and the sha256 of <platform>/SHA256SUMS at the pinned tag -- the hard
 * trust anchor. NULL anchor = not published at this tag yet.
 */
#ifndef SESAME_REGISTRY_H
#define SESAME_REGISTRY_H

#define SESAME_BASE_URL     "$base"
#define SESAME_DEFAULT_TAG  "$tag"
#define SESAME_SUMS_FILE    "SHA256SUMS"

typedef struct {
    const char *platform;
    int32_t     beads;        /* IDAT nSNPsRead; 0 = no auto-detect */
    const char *ordering;     /* ordering table filename */
    const char *sums_sha256;  /* sha256 of <platform>/SHA256SUMS; NULL = unpublished */
} sesame_reg_t;

static const sesame_reg_t SESAME_REGISTRY[] = {
EOF

  rows_of "$cat_dir/InfiniumAnnotation.tsv" |
    while IFS=$'\t' read -r p tools beads _rows ord; do
      case ",$tools," in *,sesame,*) ;; *) continue ;; esac
      printf '    { "%s", %s, "%s", "%s" },\n' \
        "$p" "$beads" "$ord" "$(anchor_of InfiniumAnnotation "$tag" "$p")"
    done

  cat <<EOF
    { NULL, 0, NULL, NULL }
};

/* --- Genome-level annotation (seqinfo/gaps/cytoband), hosted separately in
 * zhou-lab/genomes so plotting tools can reuse it. Raw git paths, one repo-wide
 * tag. Layout: <base>/<tag>/<genome>/{SHA256SUMS, seqinfo.tsv.gz, gaps.tsv.gz,
 * cytoband.tsv.gz}. Trust anchor = sha256(<genome>/SHA256SUMS) at the tag. */
#define SESAME_GENOME_BASE_URL  "$g_base"
#define SESAME_GENOME_TAG       "$g_tag"

typedef struct {
    const char *genome;
    const char *sums_sha256;  /* sha256 of <genome>/SHA256SUMS; NULL = unpublished */
} sesame_genome_reg_t;

static const sesame_genome_reg_t SESAME_GENOME_REGISTRY[] = {
EOF

  rows_of "$cat_dir/genomes.tsv" | while IFS=$'\t' read -r g _tools; do
    printf '    { "%s", "%s" },\n' "$g" "$(anchor_of genomes "$g_tag" "$g")"
  done

  cat <<'EOF'
    { NULL, NULL }
};

#endif /* SESAME_REGISTRY_H */
EOF
}

emit_kycg() {
  local tag base kb_base kb_tag
  tag=$(tag_of InfiniumAnnotation); base=$(base_of InfiniumAnnotation)
  kb_base=$(base_of KYCGKB);        kb_tag=$(tag_of KYCGKB)

  cat <<EOF
/* registry.h -- GENERATED by tools/make_registry.sh. Do not edit.
 *
 * Regenerate with:  external/YAME/tools/make_registry.sh --tool=kycg -o src/registry.h
 *
 * Arrays: sha256(<PLATFORM>/KYCG/SHA256SUMS) at the pinned tag is the trust
 * anchor; every file digest chains from it. NULL anchor = not published.
 *
 * Set counts are pinned too. A tag is immutable, so the number of sets it
 * publishes is a fixed fact about it -- there is nothing to discover at run
 * time, and the overview can show have/total without touching the network.
 *
 * Whole genome: the same, anchored on sha256(SHA256SUMS) in KYCGKB_<genome>
 * at its own tag. File sizes are pinned alongside for display only; nothing
 * depends on them, so a set added upstream needs no rebuild.
 *
 * Zenodo keeps the DOI and remains the citable archive. It is recorded here as
 * provenance and is no longer the fetch path.
 *
 * SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2026-present Wanding Zhou
 */
#ifndef _KYCG_REGISTRY_H
#define _KYCG_REGISTRY_H

#include <stdint.h>
#include <stddef.h>

#define KYCG_IA_BASE_URL   "$base"
#define KYCG_IA_TAG        "$tag"
#define KYCG_IA_SUMS_FILE  "SHA256SUMS"
#define KYCG_KB_BASE_URL   "$kb_base"
#define KYCG_ZENODO_BASE   "https://zenodo.org/records"

/* One array platform's KYCG knowledgebase directory. */
typedef struct {
    const char *platform;
    uint64_t    rows;          /* probes in the platform's ordering */
    uint64_t    n_sets;        /* .cm sets published at this tag */
    const char *sums_sha256;   /* sha256 of <platform>/KYCG/SHA256SUMS */
    /* The probe ordering sits one level up, under the platform's own
     * manifest, and is fetched alongside any set: it is what gives a row
     * index a probe name, so a set without it is a column of anonymous bits. */
    const char *plat_sums_sha256; /* sha256 of <platform>/SHA256SUMS */
} kycg_array_reg_t;

static const kycg_array_reg_t KYCG_ARRAY_REGISTRY[] = {
EOF

  rows_of "$cat_dir/InfiniumAnnotation.tsv" |
    while IFS=$'\t' read -r p tools _beads rows _ord; do
      case ",$tools," in *,kycg,*) ;; *) continue ;; esac
      printf '    { "%s", %s, %s, "%s", "%s" },\n' \
        "$p" "$rows" "$(nsets_of InfiniumAnnotation "$tag" "$p/KYCG")" \
        "$(anchor_of InfiniumAnnotation "$tag" "$p/KYCG")" \
        "$(anchor_of InfiniumAnnotation "$tag" "$p")"
    done

  cat <<'EOF'
    { NULL, 0, 0, NULL, NULL }
};

/* A published file size, for display only -- never for correctness. */
typedef struct {
    const char *name;
    uint64_t    size;
} kycg_fsize_t;

EOF

  ## One size table per genome, ASCII-ordered by name (the catalog holds them
  ## in that order already, which is what keeps this emission stable).
  rows_of "$cat_dir/KYCGKB.tsv" | cut -f1 | while read -r g; do
    printf 'static const kycg_fsize_t KYCG_SIZES_%s[] = {\n' "$g"
    rows_of "$cat_dir/file_sizes.tsv" | awk -F'\t' -v g="KYCGKB/$g" '$1==g {
      printf "    { \"%s\", %s },\n", $2, $3 }'
    printf '    { NULL, 0 }\n};\n\n'
  done

  cat <<'EOF'
/* One genome's whole-genome knowledgebase collection.
 *
 * `sums_sha256` is the trust anchor, exactly as for arrays. `record` and `doi`
 * point at the Zenodo deposit, which remains the citable archive but is no
 * longer where kycg fetches from. */
typedef struct {
    const char        *genome;
    uint64_t           rows;      /* CpGs in cpg_nocontig.cr */
    uint64_t           n_sets;    /* .cm sets published at this tag */
    const char        *repo;      /* GitHub repository name */
    const char        *tag;       /* pinned tag */
    const char        *sums_sha256;
    const char        *record;    /* Zenodo record id, provenance only */
    const char        *doi;
    const kycg_fsize_t *sizes;
} kycg_seq_reg_t;

static const kycg_seq_reg_t KYCG_SEQ_REGISTRY[] = {
EOF

  rows_of "$cat_dir/KYCGKB.tsv" |
    while IFS=$'\t' read -r g tools repo rows record doi; do
      case ",$tools," in *,kycg,*) ;; *) continue ;; esac
      printf '    { "%s", %s, %s, "%s", "%s", "%s", "%s", "%s", KYCG_SIZES_%s },\n' \
        "$g" "$rows" "$(nsets_of KYCGKB "$kb_tag" "$g")" "$repo" "$kb_tag" \
        "$(anchor_of KYCGKB "$kb_tag" "$g")" "$record" "$doi" "$g"
    done

  cat <<'EOF'
    { NULL, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL }
};

#endif /* _KYCG_REGISTRY_H */
EOF
}

# ------------------------------------------------------------------ dispatch

render() {
  case $tool in
    yame)      emit_yame ;;
    methscope) emit_yame methscope "methscope methscope_models coordinates" ;;
    kycg)   emit_kycg ;;
    sesame) emit_sesame ;;
    *) echo "make_registry.sh: unknown --tool=$tool (yame|kycg|sesame|methscope)" >&2
       exit 1 ;;
  esac
}

## Atomic write, so a failed run never leaves a half-written registry that
## still compiles.
if [ -n "$out" ]; then
  tmp=$(mktemp)
  render > "$tmp"
  mv "$tmp" "$out"
  echo "wrote $out" >&2
else
  render
fi
