#!/bin/bash
## Every command the documentation shows, actually run.
##
## Nothing regenerates docs/llms.txt or docs/index.html, so a flag that is
## renamed, a default that moves, or a command that breaks only on the path the
## docs take stays wrong until a reader hits it. `binarize -t 0.5 -c 3 -o
## calls.cg <indexed input>` wrote a SHORT file for several releases exactly
## that way: every test in this suite used the stdout form instead.
##
## The commands are read straight out of the two doc files, in the order they
## appear, and run in one sandbox. A command inherits what earlier ones wrote,
## which is how the docs read: calls.cg exists because binarize made it.
##
## A command is SKIPPED, not failed, when an input it names is not there.
## Plenty of the examples are illustrations over placeholder names (in.cg,
## a.cg, query.cg) and were never runnable. The count of both is printed, and
## a fall in the ran count is the signal that something stopped being covered.
##
## Opt-in, like layer 5: it needs the network or a populated store, and about
## 90 MB of fixtures. The release tests set YAME_TEST_DOCS=1.
set -uo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/.." && pwd)

if [ -z "${YAME_TEST_DOCS:-}" ]; then
  echo "skip: docs gate runs with YAME_TEST_DOCS=1 (release checks)"; exit 0
fi
command -v python3 >/dev/null || { echo "skip: no python3 to read the docs" >&2; exit 0; }

## Remember the caller's store before pointing YAME_DATA_HOME at the sandbox,
## so the seeding above can copy out of it.
export YAME_DATA_HOME_ORIG=${YAME_DATA_HOME:-}
d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
cd "$d"
mkdir -p store
export YAME_DATA_HOME="$d/store"

## ---- the commands, in document order ---------------------------------------
python3 - "$root/docs/llms.txt" "$root/docs/index.html" > cmds.txt <<'PY'
import sys, re, html
out = []
for p in sys.argv[1:]:
    s = open(p).read()
    if p.endswith('.html'):
        blocks = re.findall(r'<pre><code>(.*?)</code></pre>', s, re.S)
        text = "\n".join(html.unescape(re.sub(r'<[^>]+>', '', b)) for b in blocks)
    else:
        text = s
    for ln in text.split('\n'):
        t = ln.strip()
        ## a continuation or a comment line is not a command
        if not t.startswith('yame '): continue
        if t.endswith('\\'): continue          # multi-line: not handled, skipped
        out.append(t)
seen = set()
for c in out:
    if c in seen: continue
    seen.add(c); print(c)
PY
total=$(grep -c . cmds.txt)
[ "$total" -gt 20 ] || { echo "only $total commands found in the docs; the reader is broken"; exit 1; }

## ---- seed the sandbox ------------------------------------------------------
## The coordinates go in the STORE, because the docs name them as `-R hg38`.
## The data files go in the CWD, because the docs name them bare.
##
## A file already in the caller's store is COPIED rather than downloaded. The
## fixtures are about 90 MB and this runs on every release, so re-fetching them
## each time would be a download nobody needs and a burst the server may
## rate-limit. `have` is the caller's store, read but never written.
have=${YAME_DATA_HOME_ORIG:-}
## An index rides with its data file, the way a fetch delivers it. Copying the
## .cg alone left `subset` with no index and three documented commands failed.
seed_store() {                      # seed_store <registry path>
  local rel=$1
  if [ -n "$have" ] && [ -f "$have/$rel" ]; then
    mkdir -p "$(dirname "$YAME_DATA_HOME/$rel")"
    cp "$have/$rel" "$YAME_DATA_HOME/$rel" || return 1
    [ -f "$have/$rel.idx" ] && cp "$have/$rel.idx" "$YAME_DATA_HOME/$rel.idx"
    return 0
  fi
  "$YAME" fetch -y "$rel" >/dev/null 2>&1
}
seed_cwd() {                        # seed_cwd <registry path>
  local rel=$1 base=${1##*/}
  if [ -n "$have" ] && [ -f "$have/$rel" ]; then
    cp "$have/$rel" "$base" || return 1
    [ -f "$have/$rel.idx" ] && cp "$have/$rel.idx" "$base.idx"
    return 0
  fi
  "$YAME" fetch -y -c "$rel" >/dev/null 2>&1
}

seed_store hg38/cpg_nocontig.cr ||
  { echo "skip: cannot get the coordinates (offline?)" >&2; exit 0; }
for f in human_hg38_test.cg human_hg38_immune_mixture.cg human_hg38_celltypes.cg \
         human_hg38_40_celltypes_chr20.cg; do
  seed_cwd "hg38/data/$f" ||
    { echo "skip: cannot get $f (offline?)" >&2; exit 0; }
done
## The masks the docs fetch in their own examples. Those fetch lines are
## skipped below, so the commands that USE a mask would skip too without this.
seed_store hg38/KYCG/CGI.20220904.cm
seed_cwd  hg38/KYCG/Blacklist.20220304.cm

## ---- run them --------------------------------------------------------------
ran=0; skipped=0; failed=0
: > skipped.txt
while IFS= read -r cmd; do
  [ -n "$cmd" ] || continue

  ## fetch is covered by t_fetch.sh and t_fetch_http.sh, and running the
  ## examples that really download would pull whole directories. The offline
  ## forms still run, since those are the ones a reader copies to look around.
  case $cmd in
    "yame fetch"*)
      case $cmd in
        *" -l"*|*" -n"*) ;;
        *) skipped=$((skipped+1)); echo "fetch: $cmd" >> skipped.txt; continue ;;
      esac ;;
  esac

  ## Every token that looks like an input file must exist. Placeholder names
  ## (in.cg, a.cg) never will, and those examples are illustrations.
  missing=
  for tok in $(printf '%s\n' "$cmd" | grep -oE '[A-Za-z0-9_./-]+\.(cg|cm|cr|cx|txt|bed|tsv)\b'); do
    case $cmd in
      *"-o $tok"*|*"> $tok"*|*">$tok"*) continue ;;   # an output, not an input
    esac
    ## pack, rowop and dsample name their output as the LAST word, so a
    ## trailing file there is written, not read. Only those three: pairwise
    ## and mask also end in a file and that one IS an input.
    case $cmd in
      "yame pack "*|"yame rowop "*|"yame dsample "*)
        case $cmd in
          *" $tok") case $tok in *.cg|*.cx) continue ;; esac ;;
        esac ;;
    esac
    [ -e "$tok" ] || missing="$tok"
  done
  if [ -n "$missing" ]; then
    skipped=$((skipped+1)); echo "no $missing: $cmd" >> skipped.txt; continue
  fi

  ## `yame` is on no PATH here, so point the name at the binary under test.
  ## Run it in a plain subshell with NO pipefail: the docs end pipelines with
  ## `head`, which closes the pipe and leaves the producer on SIGPIPE (141).
  ## That is correct behaviour for the reader and must not read as a failure.
  real=${cmd/#yame /"$YAME" }
  if ! out=$(bash -c "$real" 2>&1 >/dev/null); then
    failed=$((failed+1))
    echo "FAIL: $cmd"
    printf '%s\n' "$out" | head -5 | sed 's/^/      /'
  else
    ran=$((ran+1))
  fi
done < cmds.txt

echo "docs: $ran ran, $skipped skipped, $failed failed (of $total)"
[ "$failed" -eq 0 ] || { echo "--- skipped, for reference:"; sed 's/^/  /' skipped.txt; exit 1; }
## A floor, so the gate cannot quietly stop covering anything.
[ "$ran" -ge 12 ] || { echo "only $ran commands ran; the gate has stopped covering the docs";
                       sed 's/^/  /' skipped.txt; exit 1; }
echo "ok: t_docs"
