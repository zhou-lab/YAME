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

## ---- 8. a stale local manifest is a pin conflict, not silently overwritten --
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
