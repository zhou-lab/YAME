## The registry, for any tool that generates its own header from it.
##
## SOURCE THIS, do not run it: it defines functions and one path and does
## nothing else. Where a consumer sources it from, and what it emits
## afterwards, are the consumer's business -- this file states the table's
## facts and stops. It is the whole of what the suite's tools share.
##
## The table is tools/registry/files.tsv: one row per file the suite can
## fetch, keyed source@tag:remote_path. The paths resolve off BASH_SOURCE
## rather than $0, because in a sourced file $0 is the caller's script.

reg=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
files=$reg/files.tsv

## Strip comments and blank lines from a TSV.
rows_of() { grep -v '^#' "$1" | grep -v '^[[:space:]]*$'; }

## The rows whose store_path matches any of the globs given, in table order;
## no globs means every row. Ordinary shell globbing: `*` and `?` stop at a
## slash, `**` crosses it. So 'EPICv2/*' is that directory's own files --
## the same set `yame fetch EPICv2` means -- and 'EPICv2/**' includes
## EPICv2/KYCG/ beneath it.
files_of() {
  local pats; pats=$(printf '%s\t' "$@")
  rows_of "$files" | /usr/bin/awk -F'\t' -v pats="$pats" '
    function g2re(g,   r, i, c, n) {
      r = "^"; n = length(g)
      for (i = 1; i <= n; i++) {
        c = substr(g, i, 1)
        if (c == "*") { if (substr(g, i+1, 1) == "*") { r = r ".*"; i++ } else r = r "[^/]*" }
        else if (c == "?") r = r "[^/]"
        else if (index(".[]()+^$|\\{}", c)) r = r "\\" c
        else r = r c
      }
      return r "$"
    }
    BEGIN { np = split(pats, p, "\t"); for (i = 1; i <= np; i++) if (p[i] != "") re[++n] = g2re(p[i]) }
    n == 0 { print; next }
    { for (i = 1; i <= n; i++) if ($2 ~ re[i]) { print; next } }'
}

## One column of one row, by name, so an emitter never re-splits a line:
##   field "$row" sha256
field() {   ## row column
  local i
  case $2 in
    key) i=1 ;; store_path) i=2 ;; sha256) i=3 ;; size) i=4 ;; rows) i=5 ;;
    kind) i=6 ;; title) i=7 ;; description) i=8 ;;
    source) i=9 ;; citation) i=10 ;;
    *) echo "registry/lib.sh: no column named $2" >&2; return 1 ;;
  esac
  printf '%s\n' "$1" | cut -f"$i"
}

## The three parts of a key, source@tag:remote_path, split on the LAST @ and
## the LAST colon: a source may carry a colon (hf:org/repo), a path never
## does. One parser, so three repos cannot each get it subtly wrong.
key_source() { local s=${1%:*}; printf '%s\n' "${s%@*}"; }   ## <org>/<repo>, or hf:<org>/<repo>
key_tag()    { local s=${1%:*}; printf '%s\n' "${s##*@}"; }  ## the immutable upstream tag
key_path()   { printf '%s\n' "${1##*:}"; }                   ## the path within that tag

## The download URL for a key, by the host rule -- the one place it lives:
##   github       https://raw.githubusercontent.com/<org>/<repo>/<tag>/<remote_path>
##   huggingface  https://huggingface.co/<org>/<repo>/resolve/<tag>/<remote_path>
url_of() {   ## key
  local src; src=$(key_source "$1")
  case $src in
    hf:*) printf 'https://huggingface.co/%s/resolve/%s/%s\n' "${src#hf:}" "$(key_tag "$1")" "$(key_path "$1")" ;;
    *)    printf 'https://raw.githubusercontent.com/%s/%s/%s\n' "$src" "$(key_tag "$1")" "$(key_path "$1")" ;;
  esac
}
