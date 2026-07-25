#!/usr/bin/env bash
## Generate a tool's compiled asset registry from the shared catalog.
##
##   tools/make_registry.sh --tool=yame   [-o src/registry.h]
##   tools/make_registry.sh --tool=kycg   [-o src/registry.h]
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

here=$(cd "$(dirname "$0")" && pwd)
reg=$here/registry
cat_dir=$reg/catalog
sums_dir=$reg/sums

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

sha256_of() {                      ## hash a file, portable across the lab's boxes
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | cut -d' ' -f1
  else
    shasum -a 256 "$1" | cut -d' ' -f1
  fi
}

## Strip comments and blank lines from a catalog TSV.
rows_of() { grep -v '^#' "$1" | grep -v '^[[:space:]]*$'; }

## A field from TAGS: tag_of <source>, base_of <source>.
tag_of()  { rows_of "$reg/TAGS" | awk -v s="$1" '$1==s {print $2}'; }
base_of() { rows_of "$reg/TAGS" | awk -v s="$1" '$1==s {print $3}'; }

## The cached manifest for one directory, and the facts derived from it.
sums_path() { echo "$sums_dir/$1/$2/$3/SHA256SUMS"; }   ## source tag subpath
anchor_of() {
  local p; p=$(sums_path "$1" "$2" "$3")
  [ -s "$p" ] || { echo "make_registry.sh: no cached manifest at $p" >&2
                   echo "  run --refresh first" >&2; exit 1; }
  sha256_of "$p"
}
nsets_of() {
  local p; p=$(sums_path "$1" "$2" "$3")
  grep -c '\.cm$' "$p" || true
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

  ## File sizes: the GitHub contents API is the only source, so this is the one
  ## part of the catalog that cannot be rebuilt from a manifest.
  echo "refreshing KYCGKB file sizes" >&2
  tmp=$(mktemp)
  {
    echo "# Published file sizes for the KYCGKB whole-genome collections, for"
    echo "# display only -- nothing depends on them, so a set added upstream needs"
    echo "# no rebuild. Refreshed from the GitHub contents API (the one fact in"
    echo "# this catalog that no manifest carries); ASCII-sorted by name, which is"
    echo "# the order the emitted table must keep."
    echo "#"
    printf "# genome\tname\tsize\n"
    rows_of "$cat_dir/KYCGKB.tsv" | while IFS=$'\t' read -r g _tools repo _rest; do
      curl -sfL "https://api.github.com/repos/zhou-lab/$repo/contents/?ref=$kb_tag" \
        | python3 -c "
import json,sys
d = json.load(sys.stdin)
for f in sorted(d, key=lambda x: x['name']):
    if f.get('type') != 'file': continue
    if f['name'] in ('README.md', 'SHA256SUMS'): continue
    print('$g\t%s\t%d' % (f['name'], f['size']))
"
    done
  } > "$tmp"
  mv "$tmp" "$cat_dir/seq_sizes.tsv"
  echo "done; review the diff under tools/registry/ before committing" >&2
  exit 0
fi

[ -n "$tool" ] || { sed -n '2,9p' "$0" | sed 's/^## \{0,1\}//' >&2; exit 1; }

# ------------------------------------------------------------------ emitters

emit_yame() {
  local ia_tag ia_base g_tag g_base kb_tag kb_base
  ia_tag=$(tag_of InfiniumAnnotation); ia_base=$(base_of InfiniumAnnotation)
  g_tag=$(tag_of genomes);             g_base=$(base_of genomes)
  kb_tag=$(tag_of KYCGKB);             kb_base=$(base_of KYCGKB)

  cat <<EOF
/* registry.h -- GENERATED by tools/make_registry.sh --tool=yame. Do not edit.
 *
 * What \`yame fetch\` knows how to download, and the digest each directory's
 * manifest must have. One row per fetchable directory: the store path is the
 * shared layout every tool in the suite reads, so a directory filled here is
 * found by kycg / sesame-cli / methscope-cli without any of them re-fetching.
 *
 * The anchor is sha256(SHA256SUMS) at the pinned tag, hashed from the copy
 * cached under tools/registry/sums/. Verifying a fetched manifest against it
 * is what makes every per-file digest trustworthy, and comparing a STORED
 * manifest against it is how a directory says which tag filled it.
 *
 * SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2021-present Wanding Zhou
 */
#ifndef _YAME_REGISTRY_H
#define _YAME_REGISTRY_H

#include <stddef.h>

typedef struct {
    const char *source;      /* upstream repo family: InfiniumAnnotation, ... */
    const char *target;      /* platform, genome, or "<platform>/KYCG" */
    const char *base_url;    /* <base>/<tag>/<remote_sub>/SHA256SUMS */
    const char *tag;
    const char *remote_sub;  /* "" when the manifest is at the repo root */
    const char *store_sub;   /* path under the store root */
    const char *anchor;      /* sha256 of that directory's SHA256SUMS */
} yame_asset_reg_t;

static const yame_asset_reg_t YAME_ASSETS[] = {
EOF

  ## Array platforms: the platform directory, then its KYCG/ subdirectory.
  ## Both are fetchable on their own -- sesame wants the first, kycg wants both.
  rows_of "$cat_dir/InfiniumAnnotation.tsv" | while IFS=$'\t' read -r p _t _b _r _o; do
    printf '    { "InfiniumAnnotation", "%s", "%s", "%s", "%s", "InfiniumAnnotation/%s", "%s" },\n' \
      "$p" "$ia_base" "$ia_tag" "$p" "$p" "$(anchor_of InfiniumAnnotation "$ia_tag" "$p")"
    printf '    { "InfiniumAnnotation", "%s/KYCG", "%s", "%s", "%s/KYCG", "InfiniumAnnotation/%s/KYCG", "%s" },\n' \
      "$p" "$ia_base" "$ia_tag" "$p" "$p" "$(anchor_of InfiniumAnnotation "$ia_tag" "$p/KYCG")"
  done

  ## Whole-genome knowledgebases: one repo each, manifest at the repo root.
  rows_of "$cat_dir/KYCGKB.tsv" | while IFS=$'\t' read -r g _tools repo _rest; do
    printf '    { "KYCGKB", "%s", "%s/%s/raw", "%s", "", "KYCGKB/%s", "%s" },\n' \
      "$g" "$kb_base" "$repo" "$kb_tag" "$g" "$(anchor_of KYCGKB "$kb_tag" "$g")"
  done

  ## Genome-level annotation (seqinfo / gaps / cytoband).
  rows_of "$cat_dir/genomes.tsv" | while IFS=$'\t' read -r g _tools; do
    printf '    { "genomes", "%s", "%s", "%s", "%s", "genomes/%s", "%s" },\n' \
      "$g" "$g_base" "$g_tag" "$g" "$g" "$(anchor_of genomes "$g_tag" "$g")"
  done

  cat <<'EOF'
    { NULL, NULL, NULL, NULL, NULL, NULL, NULL }
};

#define YAME_ASSETS_N (sizeof(YAME_ASSETS)/sizeof(YAME_ASSETS[0]) - 1)

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
    rows_of "$cat_dir/seq_sizes.tsv" | awk -F'\t' -v g="$g" '$1==g {
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
    yame)   emit_yame ;;
    kycg)   emit_kycg ;;
    sesame) emit_sesame ;;
    *) echo "make_registry.sh: unknown --tool=$tool (yame|kycg|sesame)" >&2
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
