## Catalog lookups, for any tool that generates a registry from this catalog.
##
## SOURCE THIS, do not run it: it defines functions and paths and does nothing
## else. Where a consumer sources it from, and what it emits afterwards, are
## the consumer's business -- this file states the catalog's facts and stops.
##
## YAME keeps the catalog: which tag each upstream source is pinned at, and the
## cached manifests that make an anchor verifiable with no network. It does not
## keep anyone else's emitter. A tool's registry struct is that tool's public
## type, so adding a field to it must not require a commit here. Each tool owns
## its emitter; this file is the whole of what they share.
##
## WHY THE PATHS RESOLVE OFF BASH_SOURCE. In a sourced file $0 is the caller's
## script, not this one. Only BASH_SOURCE names this file, and the catalog is
## its neighbour, so this is what makes the lookups work from any caller in any
## directory.

reg=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cat_dir=$reg/catalog
sums_dir=$reg/sums

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
  [ -s "$p" ] || { echo "registry/lib.sh: no cached manifest at $p" >&2
                   echo "  run make_registry.sh --refresh first" >&2; exit 1; }
  sha256_of "$p"
}
nsets_of() {
  local p; p=$(sums_path "$1" "$2" "$3")
  grep -c '\.cm$' "$p" || true
}

## ---- files.tsv: one row per file the suite can fetch --------------------------
files=$reg/files.tsv

## The rows whose store_path matches any of the globs given, in table order.
## Plain shell globs, and `*` crosses `/`: '*/KYCG/*' is every knowledge-base
## set, '*/cpg_nocontig.cr' every coordinate stream, 'hg38/data/*' one
## directory. No arguments means every row.
files_of() {
  local row store pat
  rows_of "$files" | while IFS= read -r row; do
    store=${row#*	}; store=${store%%	*}
    if [ $# -eq 0 ]; then printf '%s\n' "$row"; continue; fi
    for pat in "$@"; do
      case $store in $pat) printf '%s\n' "$row"; break ;; esac
    done
  done
}

## One column of one row, by name, so an emitter never re-splits a line:
##   field "$row" sha256
field() {   ## row column
  local i
  case $2 in
    key) i=1 ;; store_path) i=2 ;; sha256) i=3 ;; size) i=4 ;; rows) i=5 ;;
    kind) i=6 ;; recommend) i=7 ;; title) i=8 ;; description) i=9 ;;
    source) i=10 ;; citation) i=11 ;;
    *) echo "registry/lib.sh: no column named $2" >&2; return 1 ;;
  esac
  printf '%s\n' "$1" | cut -f"$i"
}

## The download URL for a key, by the host rule -- the one place it lives:
##   github       https://raw.githubusercontent.com/<org>/<repo>/<tag>/<remote_path>
##   huggingface  https://huggingface.co/<org>/<repo>/resolve/<tag>/<remote_path>
## The key is source@tag:remote_path; the source is <org>/<repo> on GitHub or
## hf:<org>/<repo> on HuggingFace. Split on the LAST colon: the source may
## carry one, the path never does.
url_of() {   ## key
  local srctag=${1%:*} remote=${1##*:}
  local src=${srctag%@*} tag=${srctag##*@}
  case $src in
    hf:*) printf 'https://huggingface.co/%s/resolve/%s/%s\n' "${src#hf:}" "$tag" "$remote" ;;
    *)    printf 'https://raw.githubusercontent.com/%s/%s/%s\n' "$src" "$tag" "$remote" ;;
  esac
}
