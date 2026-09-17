#!/bin/bash
## The docs against the binary, and against each other. Offline and fast, so
## it runs on every push rather than only at a release.
##
## docs/llms.txt and docs/index.html are hand-maintained and nothing
## regenerates them. So a renamed flag, a dropped subcommand or a default that
## moved stays wrong in the docs until a reader hits it. t_docs.sh catches an
## example that FAILS; this catches the ones that never run, and prose.
##
## Three checks, all mechanical:
##   1. every subcommand the binary lists is described in llms.txt
##   2. every flag llms.txt shows for a subcommand exists in that
##      subcommand's own -h; a flag that does not is simply wrong
##   3. every subcommand llms.txt describes still exists in the binary
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/.." && pwd)
llms=$root/docs/llms.txt
[ -f "$llms" ] || { echo "no docs/llms.txt"; exit 1; }

## The banner exits 1 by convention, so capture it before reading.
banner=$("$YAME" </dev/null 2>&1) || true
subs=$(printf '%s\n' "$banner" | /usr/bin/awk '/^  [a-z]/ && $1 != "yame" {print $1}' | sort -u)
[ -n "$subs" ] || { echo "could not read the subcommand list from the banner"; exit 1; }

## the subcommand names as an alternation, for the block boundary above
subs_re=$(printf '%s\n' "$subs" | paste -sd'|' -)

bad=0

## ---- 1 & 2. each subcommand: described, and its documented flags exist -----
for sub in $subs; do
  ## The block starts at "  <sub>" and ends at the NEXT command header or the
  ## next "## " heading. Bounding it by "any indented lowercase word" ran the
  ## fetch block to the end of the file and swept up every flag in the recipes.
  block=$(/usr/bin/awk -v s="$sub" -v all="$subs_re" '
    function is_header(l,   i) { return (l ~ ("^  (" all ")([ \t]|$)")) }
    $0 ~ ("^  " s "([ \t]|$)") { on = 1; print; next }
    on && (/^## / || is_header($0)) { on = 0 }
    on { print }' "$llms")

  if [ -z "$block" ]; then
    echo "llms.txt does not describe the subcommand: $sub"
    bad=1
    continue
  fi

  ## Flags the docs show for it, and flags its own -h shows.
  ## Only the USAGE line and the option lines, never prose. Prose mentions
  ## flags it is denying ("The mask is positional, not -m") and flags that
  ## belong to other commands ("-m, and -R elsewhere"), and both read as
  ## errors otherwise. An option line is one whose first character is a dash.
  ##
  ## `|| true`: a command with no documented flags makes grep exit 1, which
  ## under `set -e` would end the run rather than report anything.
  doc_flags=$(printf '%s\n' "$block" |
              /usr/bin/awk 'NR == 1 || $1 ~ /^-[A-Za-z]/' |
              grep -oE '(^|[^-A-Za-z0-9])-[A-Za-z]([^A-Za-z]|$)' |
              grep -oE '\-[A-Za-z]' | sort -u || true)
  help=$("$YAME" "$sub" -h </dev/null 2>&1) || true
  real_flags=$(printf '%s\n' "$help" | grep -oE '(^|[^-A-Za-z0-9])-[A-Za-z]([^A-Za-z]|$)' |
               grep -oE '\-[A-Za-z]' | sort -u || true)

  ## Membership by `case`, NOT by `printf | grep -q`. `grep -q` exits on the
  ## first match, printf then dies of a broken pipe, and under `pipefail` the
  ## pipeline reports failure -- so a flag that IS present reads as missing.
  ## It is a race, so it failed about one run in eight.
  real_set=" $(printf '%s ' $real_flags)"
  for f in $doc_flags; do
    case $real_set in
      *" $f "*) ;;
      *) echo "llms.txt shows $sub $f, but $sub -h does not have it"; bad=1 ;;
    esac
  done
done

## ---- 3. nothing documented has been removed from the binary ---------------
subs_set=" $(printf '%s ' $subs)"
for sub in $(/usr/bin/awk '/^  [a-z][a-z]*[ ]+[[<]/ {print $1}' "$llms" | sort -u || true); do
  case $subs_set in
    *" $sub "*) ;;
    *) echo "llms.txt describes $sub, which the binary no longer lists"; bad=1 ;;
  esac
done

[ "$bad" -eq 0 ] || { echo "the docs and the binary disagree"; exit 1; }
echo "ok: t_docs_consistency ($(printf '%s\n' "$subs" | grep -c .) subcommands)"
