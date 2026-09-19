#!/usr/bin/env bash
## Generate yame's compiled asset registry from the catalog beside this script.
##
##   tools/make_registry.sh              write src/registry.h
##   tools/make_registry.sh --check      verify it matches the catalog (CI)
##   tools/make_registry.sh --refresh [--tag=vN]
##
## THIS GENERATES YAME'S REGISTRY AND NOTHING ELSE. It used to take --tool=,
## emitting kycg's and sesame-cli's headers too -- and with them their public
## struct types, so a field on kycg_seq_reg_t was a commit in this repo. That
## is backwards: kycg depends on YAME, not the reverse. Each tool now owns its
## emitter and sources tools/registry/lib.sh for the catalog lookups.
##
## What the shared generator was actually built to fix was TAG divergence --
## kycg pinned InfiniumAnnotation v8 while sesame-cli pinned v8.1, so the two
## could not share a download. tools/registry/TAGS alone fixes that, and it
## stays: one file chooses a tag, every tool moves at its next submodule bump.
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

check=0
refresh=0
new_tag=""

for arg in "$@"; do
  case $arg in
    --check)   check=1 ;;
    --tag=*)   new_tag=${arg#--tag=} ;;
    --refresh) refresh=1 ;;
    *) echo "make_registry.sh: unknown argument $arg" >&2
       sed -n '2,6p' "$0" | sed 's/^## \{0,1\}//' >&2; exit 1 ;;
  esac
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

## One manifest from upstream into the cache, or the cached copy left alone.
## Downloads beside the target and moves into place only on 200, so a
## transient failure never deletes a manifest that was fine a moment ago --
## the old `curl || rm -f` did exactly that on a flaky link. 429 and 503 are
## "come back later" and get two more tries after a wait (HuggingFace answered
## 429 during the v10 bump); anything else will not improve by asking again.
fetch_manifest() {   ## url dest
  local url=$1 dest=$2 code try
  mkdir -p "$(dirname "$dest")"
  for try in 1 2 3; do
    code=$(curl -sL -o "$dest.part" -w '%{http_code}' "$url" || echo 000)
    if [ "$code" = 200 ]; then mv "$dest.part" "$dest"; return 0; fi
    rm -f "$dest.part"
    case $code in
      429|503) echo "  HTTP $code from $url; waiting 45 s (try $try of 3)" >&2
               sleep 45 ;;
      *)       break ;;
    esac
  done
  echo "  MISS $url (HTTP $code)" >&2
  return 1
}

if [ "$refresh" = 1 ]; then
  ia_tag=${new_tag:-$(tag_of InfiniumAnnotation)}
  ia_base=$(base_of InfiniumAnnotation)
  echo "refreshing InfiniumAnnotation @ $ia_tag" >&2
  for p in $(rows_of "$cat_dir/InfiniumAnnotation.tsv" | cut -f1); do
    for sub in "$p" "$p/KYCG"; do
      fetch_manifest "$ia_base/$ia_tag/$sub/SHA256SUMS" \
                     "$sums_dir/InfiniumAnnotation/$ia_tag/$sub/SHA256SUMS" || true
    done
  done

  g_tag=$(tag_of genomes); g_base=$(base_of genomes)
  echo "refreshing genomes @ $g_tag" >&2
  for g in $(rows_of "$cat_dir/genomes.tsv" | cut -f1); do
    fetch_manifest "$g_base/$g_tag/$g/SHA256SUMS" \
                   "$sums_dir/genomes/$g_tag/$g/SHA256SUMS" || true
  done

  kb_tag=$(tag_of KYCGKB); kb_base=$(base_of KYCGKB)
  echo "refreshing KYCGKB @ $kb_tag" >&2
  rows_of "$cat_dir/KYCGKB.tsv" | while IFS=$'\t' read -r g repo _rest; do
    fetch_manifest "$kb_base/$repo/$kb_tag/SHA256SUMS" \
                   "$sums_dir/KYCGKB/$kb_tag/$g/SHA256SUMS" || true
  done

  ms_tag=$(tag_of methscope); ms_base=$(base_of methscope)
  echo "refreshing methscope @ $ms_tag" >&2
  rows_of "$cat_dir/methscope.tsv" | while IFS=$'\t' read -r _g sub; do
    fetch_manifest "$ms_base/$ms_tag/$sub/SHA256SUMS" \
                   "$sums_dir/methscope/$ms_tag/$sub/SHA256SUMS" || true
  done

  msm_tag=$(tag_of methscope_models); msm_base=$(base_of methscope_models)
  echo "refreshing methscope_models @ $msm_tag" >&2
  ## One manifest at the repo root serves every row, so fetch it ONCE, and
  ## fail the run: a missing models manifest surfaces later as "No such file"
  ## from the emit step, which is the wrong place to learn about a 429.
  msm_url="$msm_base/$msm_tag/SHA256SUMS"
  fetch_manifest "$msm_url" "$sums_dir/methscope_models/$msm_tag/SHA256SUMS" || {
    echo "make_registry.sh: cannot fetch $msm_url." >&2
    echo "  The pin in TAGS is $msm_tag. Check the tag exists, then re-run;" >&2
    echo "  a 429 is HuggingFace rate-limiting and clears on its own." >&2
    exit 1
  }

  ## File sizes: no manifest carries them, so this is the one part of the
  ## catalog that comes from an API -- GitHub's contents API for the git
  ## sources, HuggingFace's tree API for the models. Every source, because a
  ## browser that shows a size for some files and not others reads as if the
  ## others were free.
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
    echo "# so a file added upstream needs no rebuild. Refreshed by --refresh from"
    echo "# the GitHub contents API and the HuggingFace tree API (the one fact in"
    echo "# this catalog that no manifest carries); ASCII-sorted by name within a"
    echo "# scope, the order the emitted tables must keep."
    echo "#"
    echo "# scope is <source>/<subpath>: the directory the file lives in."
    echo "#"
    printf "# scope\tname\tsize\n"

    rows_of "$cat_dir/InfiniumAnnotation.tsv" | cut -f1 | while read -r p; do
      sizes_of "InfiniumAnnotation/$p"      InfiniumAnnotation "$p"      "$ia_tag"
      sizes_of "InfiniumAnnotation/$p/KYCG" InfiniumAnnotation "$p/KYCG" "$ia_tag"
    done
    rows_of "$cat_dir/KYCGKB.tsv" | while IFS=$'\t' read -r g repo _rest; do
      sizes_of "KYCGKB/$g" "$repo" "" "$kb_tag"
    done
    rows_of "$cat_dir/genomes.tsv" | cut -f1 | while read -r g; do
      sizes_of "genomes/$g" genomes "$g" "$g_tag"
    done
    rows_of "$cat_dir/methscope.tsv" | while IFS=$'\t' read -r g sub; do
      sizes_of "methscope/$g" methscope_data "$sub" "$ms_tag"
    done
    ## The models live on HuggingFace. Its tree API lists a revision with
    ## sizes; the repo comes from TAGS, so it is still named in one place.
    hf_repo=${msm_base#https://huggingface.co/}; hf_repo=${hf_repo%/resolve}
    curl -sfL "https://huggingface.co/api/models/$hf_repo/tree/$msm_tag" |
      python3 -c "
import json, sys
try:
    d = json.load(sys.stdin)
except Exception:
    sys.exit(0)
if not isinstance(d, list): sys.exit(0)
for f in sorted(d, key=lambda x: x['path']):
    if f.get('type') != 'file': continue
    if f['path'] in ('.gitattributes', 'README.md', 'SHA256SUMS'): continue
    print('methscope/models\t%s\t%d' % (f['path'], f['size']))
"
  } > "$tmp"
  mv "$tmp" "$cat_dir/file_sizes.tsv"
  echo "done; review the diff under tools/registry/ before committing" >&2
  exit 0
fi



## Per-directory file list, straight out of the cached manifest: the names and
## digests are already there, so compiling them in costs nothing and lets a
## browser offer individual files instead of whole directories. Sizes come from
## catalog/file_sizes.tsv (the contents API); 0 means "not published", and
## nothing depends on it.
emit_file_table() {   ## slug source tag subpath scope [prefix] [skip]
  local slug=$1 src=$2 tg=$3 sub=$4 scope=${5:-} pfx=${6:-} skip=${7:-}
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
    printf '    { "%s", "%s", %s },\n' "$name" "$sha" "$size"
  done < "$p"
  printf '    { NULL, NULL, 0 }\n};\n\n'
}

## A C identifier for "<source>/<target>".
slug_of() { echo "$1/$2" | tr -c 'A-Za-z0-9\n' '_'; }

# ------------------------------------------------------------------ emitters

## Emit a registry in the yame shape (yame_asset_reg_t rows) for TOOL, over
## the SOURCES named -- every source by default, which is yame's own registry;
## a downstream tool that ships fetch over the YAME code it bundles asks for
## only the sources it consumes. One function, so the two cannot drift.
emit_yame() {
  local ia_tag ia_base g_tag g_base kb_tag kb_base ms_tag ms_base msm_tag msm_base
  ia_tag=$(tag_of InfiniumAnnotation); ia_base=$(base_of InfiniumAnnotation)
  g_tag=$(tag_of genomes);             g_base=$(base_of genomes)
  kb_tag=$(tag_of KYCGKB);             kb_base=$(base_of KYCGKB)
  ms_tag=$(tag_of methscope);          ms_base=$(base_of methscope)
  msm_tag=$(tag_of methscope_models);  msm_base=$(base_of methscope_models)

  cat <<EOF
/* registry.h -- GENERATED by tools/make_registry.sh. Do not edit.
 *
 * What \`yame fetch\` knows how to download, and the digest each directory's
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
  rows_of "$cat_dir/InfiniumAnnotation.tsv" | while IFS=$'\t' read -r p _b _r _o; do
    emit_file_table "$(slug_of InfiniumAnnotation "$p")" InfiniumAnnotation \
                    "$ia_tag" "$p" "InfiniumAnnotation/$p"
    emit_file_table "$(slug_of InfiniumAnnotation "$p/KYCG")" InfiniumAnnotation \
                    "$ia_tag" "$p/KYCG" "InfiniumAnnotation/$p/KYCG"
  done
  rows_of "$cat_dir/KYCGKB.tsv" | while IFS=$'\t' read -r g repo _rest; do
    emit_file_table "$(slug_of KYCGKB "$g")" KYCGKB "$kb_tag" "$g" "KYCGKB/$g" \
                    "" "cpg_nocontig.cr"
  done
  rows_of "$cat_dir/genomes.tsv" | while IFS=$'\t' read -r g; do
    emit_file_table "$(slug_of genomes "$g")" genomes "$g_tag" "$g" "genomes/$g"
  done
  rows_of "$cat_dir/methscope.tsv" | while IFS=$'\t' read -r g sub; do
    emit_file_table "$(slug_of methscope "$g")" methscope "$ms_tag" "$sub" "methscope/$g"
  done
  rows_of "$cat_dir/methscope_models.tsv" | while IFS=$'\t' read -r t pfx store; do
    emit_file_table "$(slug_of methscope_models "$t")" methscope_models \
                    "$msm_tag" "" "methscope/models" "$pfx"
  done

  ## Then the ancestries, for the directories that have one. A source only
  ## grows one the first time it is bumped, so most of these are absent and
  ## their rows carry NULL -- which is exactly the old behaviour.
  rows_of "$cat_dir/InfiniumAnnotation.tsv" | while IFS=$'\t' read -r p _b _r _o; do
    emit_prior_table "$(slug_of InfiniumAnnotation "$p")" InfiniumAnnotation \
                     "$ia_tag" "$p"
    emit_prior_table "$(slug_of InfiniumAnnotation "$p/KYCG")" InfiniumAnnotation \
                     "$ia_tag" "$p/KYCG"
  done
  rows_of "$cat_dir/KYCGKB.tsv" | while IFS=$'\t' read -r g repo _rest; do
    emit_prior_table "$(slug_of KYCGKB "$g")" KYCGKB "$kb_tag" "$g"
  done
  rows_of "$cat_dir/genomes.tsv" | while IFS=$'\t' read -r g; do
    emit_prior_table "$(slug_of genomes "$g")" genomes "$g_tag" "$g"
  done
  rows_of "$cat_dir/methscope.tsv" | while IFS=$'\t' read -r g sub; do
    emit_prior_table "$(slug_of methscope "$g")" methscope "$ms_tag" "$sub"
  done
  rows_of "$cat_dir/methscope_models.tsv" | while IFS=$'\t' read -r t pfx store; do
    emit_prior_table "$(slug_of methscope_models "$t")" methscope_models \
                     "$msm_tag" ""
  done

  cat <<EOF
static const yame_asset_reg_t YAME_ASSETS[] = {
EOF

  ## Array platforms: the platform directory, then its KYCG/ subdirectory.
  ## Both are fetchable on their own -- sesame wants the first, kycg wants both.
  rows_of "$cat_dir/InfiniumAnnotation.tsv" | while IFS=$'\t' read -r p _b _r _o; do
    s1=$(slug_of InfiniumAnnotation "$p"); s2=$(slug_of InfiniumAnnotation "$p/KYCG")
    printf '    { "InfiniumAnnotation", "%s", "%s", "%s", "%s", "%s", "%s", YAME_FILES_%s, YAME_NFILES(YAME_FILES_%s), %s },\n' \
      "$p" "$ia_base" "$ia_tag" "$p" "$p" "$(anchor_of InfiniumAnnotation "$ia_tag" "$p")" "$s1" "$s1" \
      "$(prior_ref "$s1" InfiniumAnnotation "$ia_tag" "$p")"
    printf '    { "InfiniumAnnotation", "%s/KYCG", "%s", "%s", "%s/KYCG", "%s/KYCG", "%s", YAME_FILES_%s, YAME_NFILES(YAME_FILES_%s), %s },\n' \
      "$p" "$ia_base" "$ia_tag" "$p" "$p" "$(anchor_of InfiniumAnnotation "$ia_tag" "$p/KYCG")" "$s2" "$s2" \
      "$(prior_ref "$s2" InfiniumAnnotation "$ia_tag" "$p/KYCG")"
  done

  ## The coordinate stream on its own. Target and store_sub are <genome>, the
  ## same spelling and the same place yame uses, so `methscope fetch
  ## hg38/cpg_nocontig.cr` and `yame fetch hg38/cpg_nocontig.cr` are the same
  ## command over the same file -- and the same one the genomes unit below
  ## writes, so fetching both costs one download.

  ## Whole-genome knowledgebases: one repo each, manifest at the repo root.
  rows_of "$cat_dir/KYCGKB.tsv" | while IFS=$'\t' read -r g repo _rest; do
    s=$(slug_of KYCGKB "$g")
    printf '    { "KYCGKB", "%s", "%s/%s", "%s", "", "%s/KYCG", "%s", YAME_FILES_%s, YAME_NFILES(YAME_FILES_%s), %s },\n' \
      "$g" "$kb_base" "$repo" "$kb_tag" "$g" "$(anchor_of KYCGKB "$kb_tag" "$g")" "$s" "$s" \
      "$(prior_ref "$s" KYCGKB "$kb_tag" "$g")"
  done

  ## Genome-level annotation (seqinfo / gaps / cytoband).
  rows_of "$cat_dir/genomes.tsv" | while IFS=$'\t' read -r g; do
    s=$(slug_of genomes "$g")
    printf '    { "genomes", "%s", "%s", "%s", "%s", "%s", "%s", YAME_FILES_%s, YAME_NFILES(YAME_FILES_%s), %s },\n' \
      "$g" "$g_base" "$g_tag" "$g" "$g" "$(anchor_of genomes "$g_tag" "$g")" "$s" "$s" \
      "$(prior_ref "$s" genomes "$g_tag" "$g")"
  done

  ## Example query methylomes -- the only source here that is data rather than
  ## annotation. The target carries the genome so the browser files them under
  ## it; the manifest sits at a differently-named path upstream.
  rows_of "$cat_dir/methscope.tsv" | while IFS=$'\t' read -r g sub; do
    s=$(slug_of methscope "$g")
    printf '    { "methscope", "%s", "%s", "%s", "%s", "%s", "%s", YAME_FILES_%s, YAME_NFILES(YAME_FILES_%s), %s },\n' \
      "$g" "$ms_base" "$ms_tag" "$sub" "$g" "$(anchor_of methscope "$ms_tag" "$sub")" "$s" "$s" \
      "$(prior_ref "$s" methscope "$ms_tag" "$sub")"
  done

  ## Model bundles: manifest at the repo root, so remote_sub is empty. Mixed
  ## genomes in one directory, so the target is not a genome name.
  rows_of "$cat_dir/methscope_models.tsv" | while IFS=$'\t' read -r t pfx store; do
    s=$(slug_of methscope_models "$t")
    printf '    { "methscope", "%s", "%s", "%s", "", "%s", "%s", YAME_FILES_%s, YAME_NFILES(YAME_FILES_%s), %s },\n' \
      "$t" "$msm_base" "$msm_tag" "$t" "$(anchor_of methscope_models "$msm_tag" "")" "$s" "$s" \
      "$(prior_ref "$s" methscope_models "$msm_tag" "")"
  done

  cat <<'EOF'
    { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0 }
};

#define YAME_ASSETS_N (sizeof(YAME_ASSETS)/sizeof(YAME_ASSETS[0]) - 1)
EOF

  ## The row-space table serves refstore.c, which identifies a reference from
  ## a file's row count and needs the genome and knowledgebase sources; a tool
  ## without them has no use for it, and its rows would name directories the
  ## tool's registry does not carry.
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
  rows_of "$cat_dir/KYCGKB.tsv" | while IFS=$'\t' read -r g _rest; do
    printf 'static const char *const YAME_REF_DIRS_%s[] = { "%s", "%s/KYCG", NULL };\n' \
      "$g" "$g" "$g"
  done
  rows_of "$cat_dir/InfiniumAnnotation.tsv" | cut -f1 | while read -r p; do
    printf 'static const char *const YAME_REF_DIRS_%s[] = { "%s", "%s/KYCG", NULL };\n' \
      "$p" "$p" "$p"
  done
  echo

  cat <<'EOF'
static const yame_ref_rows_t YAME_REF_ROWS[] = {
EOF

  rows_of "$cat_dir/KYCGKB.tsv" | while IFS=$'\t' read -r g _repo nrows _rest; do
    [ -n "${nrows:-}" ] || continue
    printf '    { "%s", "genome", %s, "%s/cpg_nocontig.cr", YAME_REF_DIRS_%s, "%s" },\n' \
      "$g" "$nrows" "$g" "$g" "$g"
  done
  rows_of "$cat_dir/InfiniumAnnotation.tsv" |
    while IFS=$'\t' read -r p _beads nrows _rest; do
      [ -n "${nrows:-}" ] || continue
      printf '    { "%s", "array", %s, "%s/%s.ordering.tsv.gz", YAME_REF_DIRS_%s, "%s" },\n' \
        "$p" "$nrows" "$p" "$p" "$p" "$p"
    done

  cat <<'EOF'
    { NULL, NULL, 0, NULL, NULL, NULL }
};

#define YAME_REF_ROWS_N (sizeof(YAME_REF_ROWS)/sizeof(YAME_REF_ROWS[0]) - 1)
EOF

  cat <<'EOF'

#endif /* _YAME_REGISTRY_H */
EOF
}

# ----------------------------------------------------------------- emission

## There is ONE destination, src/registry.h beside this script, and no way to
## name another. A registry is its own repo's file: a generator that could be
## pointed at a path was how this repo came to write, and briefly to own,
## other tools' headers. --check is the read-only form, for CI.
dest=$here/../src/registry.h

if [ "$check" = 1 ]; then
  tmp=$(mktemp); trap 'rm -f "$tmp"' EXIT
  emit_yame > "$tmp"
  if cmp -s "$tmp" "$dest"; then
    echo "src/registry.h is up to date with the catalog" >&2
  else
    echo "src/registry.h does not match the catalog; regenerate it:" >&2
    echo "  tools/make_registry.sh" >&2
    diff "$dest" "$tmp" | head -20 >&2
    exit 1
  fi
else
  ## Atomic write, so a failed run never leaves a half-written registry that
  ## still compiles.
  tmp=$(mktemp)
  emit_yame > "$tmp"
  mv "$tmp" "$dest"
  echo "wrote $dest" >&2
fi
