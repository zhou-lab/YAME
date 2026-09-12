#!/bin/bash
## fetch, against a loopback mirror. The base URLs and manifest anchors are
## compiled in, so the test serves a tree that mirrors the real path layout:
## the REAL cached manifest for one directory (tools/registry/sums, so the
## compiled anchor matches) and one real 1 KB file from it -- the only data
## fixture this suite commits. Everything fetch checks, it checks for real:
## the manifest against its anchor, the file against the manifest's digest.
##
## Needs python3 for the server and a build with libcurl; skips otherwise.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/.." && pwd)

command -v python3 >/dev/null || { echo "skip: no python3 for the mirror" >&2; exit 0; }
## -h exits 1 by convention, so capture first: piping it straight into grep
## would fail the pipeline under `set -o pipefail` whatever grep found.
help=$("$YAME" fetch -h </dev/null 2>&1) || true
printf '%s\n' "$help" | grep -q 'fetch available' ||
  { echo "skip: built without libcurl" >&2; exit 0; }

d=$(mktemp -d); trap 'kill $srv 2>/dev/null; rm -rf "$d"' EXIT
cd "$d"

## ---- the mirror tree: /zhou-lab/InfiniumAnnotation/raw/v8.1/EPIC/KYCG/ ----
scope=EPIC/KYCG; tag=v8.1
sums="$root/tools/registry/sums/InfiniumAnnotation/$tag/$scope/SHA256SUMS"
[ -f "$sums" ] || { echo "no cached manifest at $sums"; exit 1; }
tree="mirror/zhou-lab/InfiniumAnnotation/raw/$tag/$scope"
mkdir -p "$tree"
cp "$sums" "$tree/SHA256SUMS"
cp "$here/fixtures/Blacklist.20220304.cm" "$tree/"
## the manifest must actually list our file with the digest we serve
grep -q "$(sha256sum "$tree/Blacklist.20220304.cm" | cut -c1-64)  Blacklist.20220304.cm" "$tree/SHA256SUMS" ||
  { echo "fixture no longer matches the cached manifest; refresh test/fixtures/"; exit 1; }

## ---- serve it on a free loopback port --------------------------------------
port=$(python3 -c 'import socket; s=socket.socket(); s.bind(("127.0.0.1",0)); print(s.getsockname()[1])')
( cd mirror && exec python3 -m http.server --bind 127.0.0.1 "$port" ) >server.log 2>&1 &
srv=$!
for i in $(seq 1 50); do
  python3 -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:$port/', timeout=1)" 2>/dev/null && break
  sleep 0.1
done
export YAME_ASSETS_MIRROR="http://127.0.0.1:$port"
export YAME_DATA_HOME="$d/store"
mkdir -p "$YAME_DATA_HOME"
asset="$scope/Blacklist.20220304.cm"

## ---- 1. the registry lists, offline -----------------------------------------
"$YAME" fetch -l </dev/null > list.txt 2>&1 || { echo "fetch -l failed"; cat list.txt; exit 1; }
grep -q 'InfiniumAnnotation' list.txt || { echo "fetch -l does not list InfiniumAnnotation"; exit 1; }

## ---- 2. a named single file fetches, verifies, and lands in the store ------
"$YAME" fetch -y "$asset" </dev/null > f1.log 2>&1 || { echo "fetch of $asset failed"; cat f1.log; exit 1; }
got="$YAME_DATA_HOME/$asset"
[ -f "$got" ] || { echo "fetched file is not at $got"; find "$YAME_DATA_HOME" -type f; exit 1; }
cmp -s "$got" "$here/fixtures/Blacklist.20220304.cm" || { echo "fetched bytes differ from the fixture"; exit 1; }
[ -f "$YAME_DATA_HOME/$scope/SHA256SUMS" ] || { echo "no manifest written beside the file"; exit 1; }
## the mirror was actually used
grep -q "Blacklist.20220304.cm" server.log || { echo "the fetch did not go through the mirror"; cat server.log; exit 1; }

## ---- 3. fetching again is a no-op: no second download ----------------------
n1=$(grep -c 'Blacklist.20220304.cm' server.log)
"$YAME" fetch -y "$asset" </dev/null > f2.log 2>&1 || { echo "re-fetch failed"; cat f2.log; exit 1; }
n2=$(grep -c 'Blacklist.20220304.cm' server.log)
[ "$n2" -eq "$n1" ] || { echo "re-fetch downloaded the file again ($n1 -> $n2 requests)"; exit 1; }

## ---- 4. -c fetches into the current directory, no manifest -----------------
mkdir cwd && ( cd cwd && "$YAME" fetch -y -c "$asset" </dev/null > ../f3.log 2>&1 ) ||
  { echo "fetch -c failed"; cat f3.log; exit 1; }
cmp -s cwd/Blacklist.20220304.cm "$here/fixtures/Blacklist.20220304.cm" || { echo "fetch -c bytes differ"; exit 1; }
[ ! -e cwd/SHA256SUMS ] || { echo "fetch -c wrote a manifest it says it does not"; exit 1; }

## ---- 5. a file the mirror corrupts is refused, and not kept -----------------
rm -rf "$YAME_DATA_HOME"/*
printf 'not the blacklist\n' > "$tree/Blacklist.20220304.cm"
if "$YAME" fetch -y "$asset" </dev/null > f4.log 2>&1; then
  echo "a corrupted download was accepted"; cat f4.log; exit 1
fi
grep -qi 'sha256\|digest\|mismatch\|verif' f4.log || { echo "corruption refused but not explained"; cat f4.log; exit 1; }
[ ! -s "$YAME_DATA_HOME/$asset" ] || { echo "a corrupted file was left in the store"; exit 1; }
cp "$here/fixtures/Blacklist.20220304.cm" "$tree/"

## ---- 6. a name the registry does not know ----------------------------------
if "$YAME" fetch -y "$scope/NoSuchFile.cm" </dev/null > f5.log 2>&1; then
  echo "an unknown asset name was accepted"; exit 1
fi
[ -s f5.log ] || { echo "unknown asset refused silently"; exit 1; }

## ---- 7. a whole directory asks before fetching; one named file does not ----
## `fetch -c <dir>` prints the plan and exits 1 off a terminal without -y;
## a name that resolves to ONE file confirms itself (the v1.39 rule).
if "$YAME" fetch -c "$scope" </dev/null > f6.log 2>&1; then
  echo "a directory fetch ran without -y"; exit 1
fi
grep -qi 'refus\|-y' f6.log || { echo "directory refusal did not mention -y"; cat f6.log; exit 1; }
rm -rf cwd2 && mkdir cwd2 && ( cd cwd2 && "$YAME" fetch -c "$asset" </dev/null > ../f7.log 2>&1 ) ||
  { echo "a single named file needed -y"; cat f7.log; exit 1; }

## ---- 8. the non-fetching option surface ------------------------------------
## -n says what it would do and stops, successfully, downloading nothing.
before=$(grep -c 'Blacklist.20220304.cm' server.log || true)
rm -rf "$YAME_DATA_HOME"/*
"$YAME" fetch -n "$scope" </dev/null > dry.log 2>&1 || { echo "fetch -n exited non-zero"; cat dry.log; exit 1; }
[ -s dry.log ] || { echo "fetch -n said nothing"; exit 1; }
[ -z "$(find "$YAME_DATA_HOME" -type f 2>/dev/null)" ] || { echo "fetch -n wrote files"; exit 1; }
after=$(grep -c 'Blacklist.20220304.cm' server.log || true)
[ "$after" -eq "$before" ] || { echo "fetch -n downloaded something"; exit 1; }

## -l is a TSV dump: a header and one row per file, no network at all.
"$YAME" fetch -l </dev/null > tsv.txt 2>&1
head -1 tsv.txt | grep -q 'target' || { echo "fetch -l header is not the TSV header"; head -1 tsv.txt; exit 1; }
[ "$(wc -l < tsv.txt)" -gt 100 ] || { echo "fetch -l listed only $(wc -l < tsv.txt) rows"; exit 1; }

## -g filters that listing by term, and every term must match.
"$YAME" fetch -l -g KYCG </dev/null > g1.txt 2>&1
[ "$(wc -l < g1.txt)" -lt "$(wc -l < tsv.txt)" ] || { echo "-g KYCG did not narrow the listing"; exit 1; }
tail -n +2 g1.txt | grep -qv 'KYCG' && { echo "-g KYCG returned a row without KYCG"; exit 1; }
"$YAME" fetch -l -g KYCG,EPIC </dev/null > g2.txt 2>&1
[ "$(wc -l < g2.txt)" -le "$(wc -l < g1.txt)" ] || { echo "a second -g term widened the listing"; exit 1; }
## a term nothing matches
"$YAME" fetch -l -g zzzznope </dev/null > g3.txt 2>&1 || true
[ "$(tail -n +2 g3.txt | wc -l)" -eq 0 ] || { echo "-g zzzznope matched something"; exit 1; }

## -d overrides the store root
alt="$d/altstore"; mkdir -p "$alt"
"$YAME" fetch -y -d "$alt" "$asset" </dev/null > alt.log 2>&1 ||
  { echo "fetch -d failed"; cat alt.log; exit 1; }
[ -f "$alt/$asset" ] || { echo "fetch -d did not use the given root"; find "$alt" -type f; exit 1; }

## -q silences progress but still fetches
rm -rf "$YAME_DATA_HOME"/*
"$YAME" fetch -q -y "$asset" </dev/null > quiet.log 2>&1 ||
  { echo "fetch -q failed"; cat quiet.log; exit 1; }
[ -f "$YAME_DATA_HOME/$asset" ] || { echo "fetch -q did not fetch"; exit 1; }

## -f re-downloads what is already present
n_before=$(grep -c 'Blacklist.20220304.cm' server.log || true)
"$YAME" fetch -f -y "$asset" </dev/null > force.log 2>&1 ||
  { echo "fetch -f failed"; cat force.log; exit 1; }
n_after=$(grep -c 'Blacklist.20220304.cm' server.log || true)
[ "$n_after" -gt "$n_before" ] || { echo "fetch -f did not re-download"; exit 1; }

## ---- 9. the single-file form: -u with -s and -o -----------------------------
## A URL, the digest it must have, and where it goes. This path shares the
## download and verify code but none of the registry.
url="$YAME_ASSETS_MIRROR/zhou-lab/InfiniumAnnotation/raw/$tag/$scope/Blacklist.20220304.cm"
want=$(sha256sum "$here/fixtures/Blacklist.20220304.cm" | cut -c1-64)
"$YAME" fetch -u "$url" -s "$want" -o direct.cm </dev/null > u1.log 2>&1 ||
  { echo "fetch -u failed"; cat u1.log; exit 1; }
cmp -s direct.cm "$here/fixtures/Blacklist.20220304.cm" || { echo "fetch -u bytes differ"; exit 1; }
## the wrong digest is refused and the file not kept
rm -f bad.cm
if "$YAME" fetch -u "$url" -s "${want%??}00" -o bad.cm </dev/null > u2.log 2>&1; then
  echo "fetch -u accepted a wrong digest"; exit 1
fi
[ ! -s bad.cm ] || { echo "fetch -u left a file that failed its digest"; exit 1; }
## a URL that is not there
if "$YAME" fetch -u "$YAME_ASSETS_MIRROR/nope/missing.cm" -s "$want" -o gone.cm </dev/null > u3.log 2>&1; then
  echo "fetch -u accepted a 404"; exit 1
fi

## ---- 10. -t and -k: an overriding tag has no compiled anchor ---------------
## Without -k an unpinned tag must be refused rather than fetched blind.
rm -rf "$YAME_DATA_HOME"/*
if "$YAME" fetch -y -t v0.0 "$asset" </dev/null > t1.log 2>&1; then
  echo "an unpinned -t tag was fetched without -k"; exit 1
fi
grep -qiE 'tag|anchor|pin|-k' t1.log || { echo "the -t refusal does not explain itself"; cat t1.log; exit 1; }

## ---- 11. the browser's fetch confirmation ----------------------------------
## The tree lists DIRECTORIES, not files, so f always proposes a whole unit.
## Declining is the half that needs no bytes: f must show the plan and ask,
## and n must return to the browser having downloaded nothing.
rm -rf "$YAME_DATA_HOME"/*
python3 - "$YAME" <<'PY' || { echo "the fetch confirmation did not behave"; exit 1; }
import os, pty, sys, time, select, signal
YAME = sys.argv[1]
pid, fd = pty.fork()
if pid == 0:
    os.environ["TERM"] = "xterm"; os.environ["COLUMNS"] = "110"; os.environ["LINES"] = "34"
    os.execv(YAME, ["yame", "fetch"])
out = b""
def pump(t, quiet=0.03, until=None):
    """Read until quiet for `quiet` s, or `until` (bytes) appears, or `t` s."""
    global out
    end = time.time() + t; last = time.time()
    while time.time() < end:
        if until and until in out: return
        r, _, _ = select.select([fd], [], [], 0.01)
        if r:
            try: out += os.read(fd, 65536); last = time.time()
            except OSError: return
        elif time.time() - last >= quiet and out and not until: return
pump(2.0)
for k in [b"/"] + [bytes([c]) for c in b"Blacklist"] + [b"\r"]:
    os.write(fd, k); pump(0.5)
os.write(fd, b" "); pump(0.5)              # select the highlighted directory
mark = len(out)
os.write(fd, b"f"); pump(3.0, until=b"Proceed")              # propose the fetch
prompt = out[mark:].decode("utf8", "replace")
if "Proceed" not in prompt:
    sys.exit("f did not ask before fetching:\n" + prompt[-300:])
if "Fetch" not in prompt:
    sys.exit("the prompt does not say what it would fetch")
os.write(fd, b"n"); pump(1.0)              # decline
os.write(fd, b"q"); pump(1.0)
st = None; t0 = time.time()
while time.time() - t0 < 15:
    p, s_ = os.waitpid(pid, os.WNOHANG)
    if p: st = s_; break
    pump(0.05, quiet=0.01)
if st is None:
    os.kill(pid, signal.SIGKILL); os.waitpid(pid, 0)
    sys.exit("the browser did not exit after declining")
if os.WIFSIGNALED(st):
    sys.exit(f"the browser died on signal {os.WTERMSIG(st)}")
PY
[ -z "$(find "$YAME_DATA_HOME" -type f 2>/dev/null)" ] ||
  { echo "declining the prompt still downloaded something"; exit 1; }

## ---- 11b. and accepting it, when a real unit can be mirrored ---------------
## The progress and settle rendering only runs during an actual download. The
## smallest unit in the catalogue is Mammal40/KYCG at 0.3 MB, which can be
## mirrored from the lab's shared store when it is present; where it is not
## (CI), this half skips and the rest of the file still runs.
SHARED=${YAME_SHARED_STORE:-/mnt/isilon/zhou_lab/projects/20191221_references/YAME}
if [ -d "$SHARED/Mammal40/KYCG" ] && [ -f "$SHARED/Mammal40/KYCG/SHA256SUMS" ]; then
  ## Selecting the unit takes its knowledgebase AND the platform directory
  ## above it -- "the index at the top of the list is fetched with anything
  ## else taken from this unit" -- so mirror both, or the run reports failures.
  m2="mirror/zhou-lab/InfiniumAnnotation/raw/$tag/Mammal40"
  mkdir -p "$m2/KYCG"
  ## every regular file, not a glob: SHA256SUMS has no dot in its name, and
  ## without the manifest the whole directory fails verification
  find "$SHARED/Mammal40" -maxdepth 1 -type f -exec cp {} "$m2"/ \;
  find "$SHARED/Mammal40/KYCG" -maxdepth 1 -type f -exec cp {} "$m2/KYCG"/ \;
  rm -rf "$YAME_DATA_HOME"/*
  python3 - "$YAME" <<'PY' || { echo "the browser fetch did not complete"; exit 1; }
import os, pty, sys, time, select
YAME = sys.argv[1]
pid, fd = pty.fork()
if pid == 0:
    os.environ["TERM"] = "xterm"; os.environ["COLUMNS"] = "110"; os.environ["LINES"] = "34"
    os.execv(YAME, ["yame", "fetch"])
out = b""
def pump(t, quiet=0.03, until=None):
    """Read until quiet for `quiet` s, or `until` (bytes) appears, or `t` s."""
    global out
    end = time.time() + t; last = time.time()
    while time.time() < end:
        if until and until in out: return
        r, _, _ = select.select([fd], [], [], 0.01)
        if r:
            try: out += os.read(fd, 65536); last = time.time()
            except OSError: return
        elif time.time() - last >= quiet and out and not until: return
pump(2.0)
for k in [b"/"] + [bytes([c]) for c in b"Mammal40"] + [b"\r"]:
    os.write(fd, k); pump(0.5)
os.write(fd, b" "); pump(0.5)
os.write(fd, b"f"); pump(3.0, until=b"Proceed")
mark = len(out)
os.write(fd, b"y"); pump(30.0, until=b"press any key")   # accept; wait for the summary, not a timer
report = out[mark:].decode("utf8", "replace")
if "fetched" not in report:
    sys.exit("no completion line after the fetch:\n" + report[-400:])
if "failed" in report:
    sys.exit("the browser fetch reported failures:\n" + report[-400:])
## It ends on "press any key to return to the catalogue", so the first key
## goes to that screen and only the SECOND one reaches the browser.
os.write(fd, b" "); pump(1.0)
os.write(fd, b"q"); pump(1.0)
st = None; t0 = time.time()
while time.time() - t0 < 30:
    p, s_ = os.waitpid(pid, os.WNOHANG)
    if p: st = s_; break
    pump(0.05, quiet=0.01)
if st is None:
    os.kill(pid, 9); os.waitpid(pid, 0); sys.exit("the browser hung during the fetch")
if os.WIFSIGNALED(st): sys.exit(f"the browser died on signal {os.WTERMSIG(st)}")
PY
  [ -n "$(find "$YAME_DATA_HOME" -name '*.cm' 2>/dev/null | head -1)" ] ||
    { echo "the browser fetch landed no files"; find "$YAME_DATA_HOME" -type f | head; exit 1; }
  ## whatever it fetched must verify against the manifest it wrote
  for f in $(find "$YAME_DATA_HOME" -name SHA256SUMS); do
    ( cd "$(dirname "$f")" && sha256sum -c --quiet SHA256SUMS 2>/dev/null ) ||
      { echo "a browser-fetched directory does not verify against its own manifest"; exit 1; }
  done
fi

## ---- 12. a store that is behind says so, once, on stderr --------------------
## Bare `fetch` and `fetch -l` print one line per directory whose manifest is
## from an earlier tag this binary knows, naming the command that repairs it.
## A named fetch prints nothing extra, since it IS the repair, and -q is quiet.
##
## Staged with the real mechanism: methscope_models carries prior anchors and
## its earlier manifests are cached, so a directory holding a cached prior
## manifest is exactly a store filled at that earlier tag.
rm -rf "$YAME_DATA_HOME"/*
"$YAME" fetch -l </dev/null 2> quiet.err >/dev/null
grep -q 'earlier tag' quiet.err && { echo "an empty store was reported as behind"; cat quiet.err; exit 1; }

cur=$(grep -o 'methscope_models *v[0-9]*' "$root/tools/registry/TAGS" | /usr/bin/awk '{print $2}')
prev="v$(( ${cur#v} - 1 ))"
prevsums="$root/tools/registry/sums/methscope_models/$prev/SHA256SUMS"
if [ -f "$prevsums" ]; then
  mkdir -p "$YAME_DATA_HOME/hg38/models"
  cp "$prevsums" "$YAME_DATA_HOME/hg38/models/"
  "$YAME" fetch -l </dev/null 2> behind.err >/dev/null
  grep -q "hg38/models was fetched at an earlier tag than this yame; run: yame fetch hg38/models" behind.err ||
    { echo "-l did not report the $prev-filled hg38/models as behind"; cat behind.err; exit 1; }
  [ "$(grep -c 'hg38/models' behind.err)" -eq 1 ] ||
    { echo "hg38/models was reported more than once"; cat behind.err; exit 1; }
  "$YAME" fetch </dev/null 2> bare.err >/dev/null
  grep -q 'earlier tag' bare.err || { echo "bare fetch did not report the store as behind"; exit 1; }
  ## naming any target suppresses the report: the fetch is the repair
  "$YAME" fetch -n "$asset" </dev/null 2> named.err >/dev/null
  grep -q 'earlier tag' named.err && { echo "a named fetch printed the report"; cat named.err; exit 1; }
  "$YAME" fetch -q -l </dev/null 2> q.err >/dev/null
  grep -q 'earlier tag' q.err && { echo "-q did not silence the report"; exit 1; }
else
  echo "skip: no cached $prev manifest to stage an old-tag store" >&2
fi
rm -rf "$YAME_DATA_HOME"/*

## ---- 13. a stale local manifest is a pin conflict, not silently overwritten --
rm -rf "$YAME_DATA_HOME"/*
mkdir -p "$YAME_DATA_HOME/$scope"
printf '%s  Blacklist.20220304.cm\n' "$(printf 'x%.0s' $(seq 1 64))" > "$YAME_DATA_HOME/$scope/SHA256SUMS"
if "$YAME" fetch -y "$asset" </dev/null > f8.log 2>&1; then
  ## acceptable only if it says what it did about the stale manifest
  grep -qi 'manifest\|pin\|stale\|conflict\|upgrad' f8.log ||
    { echo "a stale manifest was overwritten without a word"; cat f8.log; exit 1; }
else
  grep -qi 'manifest\|pin\|conflict\|-f' f8.log || { echo "pin conflict refused but not explained"; cat f8.log; exit 1; }
fi
