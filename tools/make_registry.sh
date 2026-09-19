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

## ---- PHASE 1 ONLY. Deleted at the start of phase 2. -------------------------
## The compiled record is still the per-directory yame_asset_reg_t, whose
## anchor is sha256 of the UPSTREAM manifest, so this script still keeps
## those manifests cached under sums/ and still reads the tag and base from
## TAGS to refresh them. Nothing downstream may build against these: they
## are not in lib.sh, and the file-centric record that replaces them needs
## none of it. (Decision 3, 2026-09-19.)
cat_dir=$reg/catalog
sums_dir=$reg/sums
sha256_of() {
  if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | cut -d' ' -f1
  else shasum -a 256 "$1" | cut -d' ' -f1; fi
}
tag_of()  { rows_of "$reg/TAGS" | awk -v s="$1" '$1==s {print $2}'; }
base_of() { rows_of "$reg/TAGS" | awk -v s="$1" '$1==s {print $3}'; }
sums_path() { echo "$sums_dir/$1/$2/$3/SHA256SUMS"; }   ## source tag subpath
anchor_of() {
  local p; p=$(sums_path "$1" "$2" "$3")
  [ -s "$p" ] || { echo "make_registry.sh: no cached manifest at $p" >&2
                   echo "  run --refresh first" >&2; exit 1; }
  sha256_of "$p"
}
## ---------------------------------------------------------------------------

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



## Per-directory file list, straight out of files.tsv: every row whose
## store_path sits in this directory, in table order.
emit_file_table() {   ## slug store_dir
  printf 'static const yame_asset_file_t YAME_FILES_%s[] = {\n' "$1"
  rows_of "$files" | /usr/bin/awk -F'\t' -v d="$2" '
    { p=$2; sub(/\/[^\/]*$/, "", p) }
    p==d { n=$2; sub(/^.*\//, "", n); printf "    { \"%s\", \"%s\", %s },\n", n, $3, ($4==""?0:$4) }'
  printf '    { NULL, NULL, 0 }\n};\n\n'
}

## A C identifier for "<source>/<target>".
slug_of() { echo "$1/$2" | tr -c 'A-Za-z0-9\n' '_'; }

# ------------------------------------------------------------------ emitters

## Emit a registry in the yame shape (yame_asset_reg_t rows) for TOOL, over
## the SOURCES named -- every source by default, which is yame's own registry;
## a downstream tool that ships fetch over the YAME code it bundles asks for
## only the sources it consumes. One function, so the two cannot drift.
## One directory of files.tsv, as the fields today's yame_asset_reg_t row
## needs. Everything comes from the key (source@tag:remote_path) and the
## store_path; the four TRANSITION rules below exist only so that, until
## Phase 2 changes the compiled record, the output stays byte-identical to
## what the per-source catalog produced -- they name the catalog source for
## the cached manifest under sums/, and the label/target/slug spellings the
## old rows carried. They leave with the old record.
unit_fields() {   ## store_dir first_key -> sets: src tag repo catsrc label base rsub target subpath slug
  local key=$2
  local srctag=${key%:*} remote=${key##*:}   ## split on the LAST colon: hf:org/repo@tag:path
  src=${srctag%@*}; tag=${srctag##*@}
  repo=${src#hf:}
  rsub=${remote%/*}; [ "$rsub" = "$remote" ] && rsub=""
  case $src in
    hf:*)                   base="https://huggingface.co/$repo/resolve" ;;
    *)                      base="https://raw.githubusercontent.com/$repo" ;;
  esac
  case $repo in                                  ## TRANSITION: catalog source names
    zhou-lab/InfiniumAnnotation) catsrc=InfiniumAnnotation; label=InfiniumAnnotation ;;
    zhou-lab/KYCGKB_*)           catsrc=KYCGKB;             label=KYCGKB ;;
    zhou-lab/genomes)            catsrc=genomes;            label=genomes ;;
    zhou-lab/methscope_data)     catsrc=methscope;          label=methscope ;;
    zhou-lab/methscope)          catsrc=methscope_models;   label=methscope ;;
    *) echo "make_registry.sh: no transition rule for $repo" >&2; exit 1 ;;
  esac
  target=$1
  [ "$catsrc" = KYCGKB ] && target=${1%/KYCG}   ## TRANSITION: the old row named the genome
  subpath=$rsub
  [ "$catsrc" = KYCGKB ] && subpath=$target      ## TRANSITION: sums/KYCGKB/<tag>/<genome>
  slug=$(slug_of "$catsrc" "$target")
}

## Store directories in table order, each with its first key: one line per
## directory, "<store_dir>\t<key>". The one-source@tag-per-directory rule is
## checked here, since this is the one place every row of a directory passes.
units_of() {
  rows_of "$files" | /usr/bin/awk -F'\t' '
    { d=$2; sub(/\/[^\/]*$/, "", d); st=$1; sub(/:[^:]*$/, "", st) }
    !(d in first) { first[d]=st; order[++n]=d; key[d]=$1; next }
    first[d]!=st  { printf "make_registry.sh: %s mixes %s and %s\n", d, first[d], st > "/dev/stderr"; bad=1 }
    END { if (bad) exit 1; for (i=1;i<=n;i++) print order[i] "\t" key[order[i]] }'
}

emit_yame() {
  local units; units=$(units_of) || exit 1
  local src tag repo catsrc label base rsub target subpath slug

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
  printf '%s\n' "$units" | while IFS=$'\t' read -r d k; do
    unit_fields "$d" "$k"
    emit_file_table "$slug" "$d"
  done

  ## Then the ancestries, for the directories that have one. A source only
  ## grows one the first time it is bumped, so most of these are absent and
  ## their rows carry NULL -- which is exactly the old behaviour.
  printf '%s\n' "$units" | while IFS=$'\t' read -r d k; do
    unit_fields "$d" "$k"
    emit_prior_table "$slug" "$catsrc" "$tag" "$subpath"
  done

  cat <<EOF
static const yame_asset_reg_t YAME_ASSETS[] = {
EOF

  printf '%s\n' "$units" | while IFS=$'\t' read -r d k; do
    unit_fields "$d" "$k"
    printf '    { "%s", "%s", "%s", "%s", "%s", "%s", "%s", YAME_FILES_%s, YAME_NFILES(YAME_FILES_%s), %s },\n' \
      "$label" "$target" "$base" "$tag" "$rsub" "$d" \
      "$(anchor_of "$catsrc" "$tag" "$subpath")" "$slug" "$slug" \
      "$(prior_ref "$slug" "$catsrc" "$tag" "$subpath")"
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
  ## the unit's own directory before its knowledge base, because that is the
  ## order someone naming a file means them in. The row-space files are the
  ## rows of files.tsv with a kind; their name is the store directory that
  ## holds them. Genomes first, then arrays, each in table order.
  refs() { rows_of "$files" | /usr/bin/awk -F'\t' -v k="$1" '$6==k { d=$2; sub(/\/[^\/]*$/, "", d); print d "\t" $5 "\t" $2 }'; }
  for kind in genome array; do
    refs $kind | while IFS=$'\t' read -r name nrows path; do
      printf 'static const char *const YAME_REF_DIRS_%s[] = { "%s", "%s/KYCG", NULL };\n' \
        "$name" "$name" "$name"
    done
  done
  echo

  cat <<'EOF'
static const yame_ref_rows_t YAME_REF_ROWS[] = {
EOF

  for kind in genome array; do
    refs $kind | while IFS=$'\t' read -r name nrows path; do
      printf '    { "%s", "%s", %s, "%s", YAME_REF_DIRS_%s, "%s" },\n' \
        "$name" "$kind" "$nrows" "$path" "$name" "$name"
    done
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
