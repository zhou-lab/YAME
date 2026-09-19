#!/usr/bin/env bash
## Generate yame's compiled registry from tools/registry/files.tsv.
##
##   tools/make_registry.sh              write src/registry.h
##   tools/make_registry.sh --check      verify it matches files.tsv (CI)
##   tools/make_registry.sh --refresh    fill blank digests and sizes from
##                                       upstream, verify the filled ones
##
## THIS GENERATES YAME'S REGISTRY AND NOTHING ELSE. Each downstream tool owns
## its own emitter and sources tools/registry/lib.sh for the table; nothing
## here can be pointed at another repo.
##
## files.tsv is the whole input: one row per file, keyed source@tag:remote_path.
## Emission is offline. --refresh is the only mode that touches the network,
## and it never overwrites a filled cell.
set -euo pipefail

## The catalog lookups live beside the catalog, in a file any tool can source.
## Sourcing sets reg, cat_dir and sums_dir as well.
here=$(cd "$(dirname "$0")" && pwd)
. "$here/registry/lib.sh"

check=0
refresh=0
for arg in "$@"; do
  case $arg in
    --check)   check=1 ;;
    --refresh) refresh=1 ;;
    *) echo "make_registry.sh: unknown argument $arg" >&2
       sed -n '2,7p' "$0" | sed 's/^## \{0,1\}//' >&2; exit 1 ;;
  esac
done

## ---- --refresh: fill blank digests and sizes, verify the filled ones ---------
## One upstream manifest per (source@tag, remote directory); every row at that
## address is filled or checked against it. A digest that disagrees is an
## error -- a git tag is immutable -- and nothing is written. Upstream files
## with no row are ignored: the table is the authority for what is offered.
if [ "$refresh" = 1 ]; then
  tmp=$(mktemp); trap 'rm -f "$tmp" "$tmp.new" "$tmp.sums" "$tmp.sizes"' EXIT
  cp "$files" "$tmp"
  bad=0
  rows_of "$files" | /usr/bin/awk -F'\t' '
    { k=$1; sub(/:[^:]*$/, "", k); p=$1; sub(/^.*:/, "", p)
      d=p; if (d ~ /\//) sub(/\/[^\/]*$/, "", d); else d=""; print k "\t" d }' | sort -u > "$tmp.addr"
  while IFS=$'\t' read -r st rdir; do
    src=${st%@*}; tag=${st##*@}
    case $src in
      hf:*) murl="https://huggingface.co/${src#hf:}/resolve/$tag/SHA256SUMS"
            surl="https://huggingface.co/api/models/${src#hf:}/tree/$tag" ;;
      *)    murl="https://raw.githubusercontent.com/$src/$tag/${rdir:+$rdir/}SHA256SUMS"
            surl="repos/$src/contents/$rdir?ref=$tag" ;;
    esac
    echo "refreshing $st${rdir:+ $rdir}" >&2
    if ! curl -sfL --retry 3 --retry-delay 30 -o "$tmp.sums" "$murl"; then
      echo "  cannot fetch $murl; rows at this address left as they were" >&2; continue
    fi
    : > "$tmp.sizes"     ## sizes are best effort: the one fact no manifest carries
    case $src in
      hf:*) curl -sfL "$surl" 2>/dev/null | python3 -c '
import json,sys
try:
    for f in json.load(sys.stdin):
        if f.get("type")=="file": print("%s\t%d" % (f["path"], f["size"]))
except Exception: pass' > "$tmp.sizes" || true ;;
      *) { command -v gh >/dev/null 2>&1 && gh api "$surl" || curl -sfL "https://api.github.com/$surl"; } 2>/dev/null | python3 -c '
import json,sys
try: d=json.load(sys.stdin)
except Exception: sys.exit(0)
for f in (d if isinstance(d,list) else []):
    if f.get("type")=="file": print("%s\t%d" % (f["name"], f["size"]))' > "$tmp.sizes" || true ;;
    esac
    /usr/bin/awk -F'\t' -v OFS='\t' -v st="$st" -v rdir="$rdir" -v sums="$tmp.sums" -v sizes="$tmp.sizes" '
      BEGIN {
        while ((getline l < sums) > 0) { i=index(l, "  "); if (i) sha[substr(l, i+2)] = substr(l, 1, i-1) }
        while ((getline l < sizes) > 0) { i=index(l, "\t"); if (i) size[substr(l, 1, i-1)] = substr(l, i+1) }
      }
      /^#/ || NF < 2 { print; next }
      { k=$1; sub(/:[^:]*$/, "", k); p=$1; sub(/^.*:/, "", p); d=p; if (d ~ /\//) sub(/\/[^\/]*$/, "", d); else d="" }
      k != st || d != rdir { print; next }
      { n=p; sub(/^.*\//, "", n)
        if (!(n in sha)) { printf "  %s is not in the manifest at %s\n", $1, st > "/dev/stderr"; bad=1 }
        else if ($3 == "") $3 = sha[n]
        else if ($3 != sha[n]) { printf "  DIGEST MISMATCH for %s: table %s, upstream %s\n", $1, $3, sha[n] > "/dev/stderr"; bad=1 }
        if ($4 == "" && (n in size)) $4 = size[n]
        print }
      END { exit bad }' "$tmp" > "$tmp.new" || bad=1
    mv "$tmp.new" "$tmp"
  done < "$tmp.addr"
  rm -f "$tmp.addr"
  if [ "$bad" != 0 ]; then echo "make_registry.sh: --refresh found problems; files.tsv NOT written" >&2; exit 1; fi
  if cmp -s "$tmp" "$files"; then echo "files.tsv: nothing to fill" >&2
  else cp "$tmp" "$files"; echo "files.tsv: filled; review the diff before committing" >&2; fi
  exit 0
fi

## The one-source@tag-per-directory rule, checked here because this is the one
## place every row of a directory passes. Exits 1 on a mix.
check_one_tag_per_dir() {
  rows_of "$files" | /usr/bin/awk -F'\t' '
    { d=$2; sub(/\/[^\/]*$/, "", d); st=$1; sub(/:[^:]*$/, "", st) }
    !(d in first) { first[d]=st; next }
    first[d]!=st  { printf "make_registry.sh: %s mixes %s and %s\n", d, first[d], st > "/dev/stderr"; bad=1 }
    END { exit bad }'
}

## A key splits on its LAST colon, so a remote path may not contain one: it
## would parse as part of the source and surface as a wrong URL, not an error.
## (kycg checked all 343 rows on 2026-09-19; this keeps it true.)
check_no_colon_in_path() {
  rows_of "$files" | /usr/bin/awk -F'\t' '
    { k=$1; sub(/^hf:/, "", k); n=gsub(/:/, ":", k) }
    n != 1 { printf "make_registry.sh: key has %d colons after the hf: prefix, want exactly 1: %s\n", n, $1 > "/dev/stderr"; bad=1 }
    END { exit bad }'
}

## src/registry.h: the table, one C row per files.tsv row, in table order;
## then the row-space fingerprints for the files that carry rows/kind.
emit_yame() {
  check_one_tag_per_dir || exit 1
  check_no_colon_in_path || exit 1
  cat <<'EOF'
/* registry.h -- GENERATED by tools/make_registry.sh. Do not edit.
 *
 * What `yame fetch` knows how to download: one row per file, keyed
 * source@tag:remote_path, with where it lands, the digest it must have, its
 * size, and what to call it. Generated from tools/registry/files.tsv and
 * nothing else. There is no directory here -- a directory is the files whose
 * store_path sits in it -- and no anchor: each file is verified against its
 * own digest, and the SHA256SUMS a fetch writes beside it records what was
 * verified.
 *
 * SPDX-License-Identifier: AGPL-3.0-or-later
 * Copyright (C) 2021-present Wanding Zhou
 */
#ifndef _YAME_REGISTRY_H
#define _YAME_REGISTRY_H

#include <stddef.h>
#include <stdint.h>

#include "assets.h"   /* yame_asset_file_t */

static const yame_asset_file_t YAME_FILES[] = {
EOF
  rows_of "$files" | /usr/bin/awk -F'\t' '
    function q(s,   i, c, r) {
      r = ""
      for (i = 1; i <= length(s); i++) { c = substr(s, i, 1); if (c == "\\" || c == "\"") r = r "\\" c; else r = r c }
      return "\"" r "\""
    }
    function url(key,   src, tag, path, st) {
      path=key; sub(/^.*:/, "", path); st=key; sub(/:[^:]*$/, "", st)
      src=st; sub(/@[^@]*$/, "", src); tag=st; sub(/^.*@/, "", tag)
      if (src ~ /^hf:/) { sub(/^hf:/, "", src); return "https://huggingface.co/" src "/resolve/" tag "/" path }
      return "https://raw.githubusercontent.com/" src "/" tag "/" path
    }
    { printf "    { %s, %s, %s, %s, %s, %s, %s, %s, %s },\n",
             q($1), q($2), q(url($1)), q($3), ($4==""?0:$4), q($7), q($8), q($9), q($10) }'
  cat <<'EOF'
    { NULL, NULL, NULL, NULL, 0, NULL, NULL, NULL, NULL }
};

#define YAME_FILES_N (sizeof(YAME_FILES)/sizeof(YAME_FILES[0]) - 1)

/**
 * Row-space fingerprints: how many rows a file written against each reference
 * has.
 *
 * A CX record carries its row count and nothing about what those rows ARE, so
 * the count is the only handle on which reference it belongs to -- and it is
 * enough, because no two row spaces in the table share one. That is what
 * lets a command infer -R rather than make the caller repeat what the file
 * already implies.
 *
 * `rows` is hand-pinned in files.tsv: it is a property of a file's contents,
 * which no manifest carries. `fetch` is the store directory, which is also
 * the argument that would put the reference in the store.
 */
typedef struct yame_ref_rows_s {
    const char *name;        /* hg38, EPIC, ...: the store directory */
    const char *kind;        /* "genome" (a .cr) or "array" (an ordering) */
    uint64_t    rows;
    const char *store_path;  /* the reference itself, under the store root */
    const char *const *dirs; /* every directory it owns, nearest first */
    const char *fetch;       /* the store directory, to give `yame fetch` */
} yame_ref_rows_t;

EOF
  refs() { rows_of "$files" | /usr/bin/awk -F'\t' -v k="$1" '$6==k { d=$2; sub(/\/[^\/]*$/, "", d); print d "\t" $5 "\t" $2 }'; }
  for kind in genome array; do
    refs $kind | while IFS=$'\t' read -r name nrows path; do
      printf 'static const char *const YAME_REF_DIRS_%s[] = { "%s", "%s/KYCG", NULL };\n' "$name" "$name" "$name"
    done
  done
  echo
  cat <<'EOF'
static const yame_ref_rows_t YAME_REF_ROWS[] = {
EOF
  for kind in genome array; do
    refs $kind | while IFS=$'\t' read -r name nrows path; do
      printf '    { "%s", "%s", %s, "%s", YAME_REF_DIRS_%s, "%s" },\n' "$name" "$kind" "$nrows" "$path" "$name" "$name"
    done
  done
  cat <<'EOF'
    { NULL, NULL, 0, NULL, NULL, NULL }
};

#define YAME_REF_ROWS_N (sizeof(YAME_REF_ROWS)/sizeof(YAME_REF_ROWS[0]) - 1)

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
