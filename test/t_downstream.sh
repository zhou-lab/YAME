#!/bin/bash
## fetch over a downstream tool's own registry. probe_fetch.c is a complete
## downstream tool -- a registry of one directory and one file, a name, a
## store variable -- linked against libyame.a. Everything fetch prints must
## speak in that tool's voice and list that tool's catalogue, and the store
## must resolve through that tool's variable ahead of YAME_DATA_HOME.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/.." && pwd)
[ -f "$root/libyame.a" ] && [ -x "$root/yame-config" ] || { echo "skip: no libyame.a (run 'make lib')" >&2; exit 0; }

d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
cd "$d"
${CC:-cc} -O1 -g -std=gnu99 $("$root/yame-config" --cflags) -o probe_fetch "$here/probe_fetch.c" \
  $("$root/yame-config" --libs) 2>cc.err || { echo "probe_fetch did not build"; cat cc.err; exit 1; }

export METHPROBE_DATA_HOME="$d/mpstore"; mkdir -p "$METHPROBE_DATA_HOME"
export YAME_DATA_HOME="$d/ystore";        mkdir -p "$YAME_DATA_HOME"

## ---- 1. the usage speaks as the tool, and lists its store variable first --
./probe_fetch -h </dev/null > h.txt 2>&1 || true
grep -q '^  methprobe fetch ' h.txt || { echo "usage does not say 'methprobe fetch'"; head -3 h.txt; exit 1; }
grep -q 'yame fetch' h.txt && { echo "usage still says 'yame fetch' somewhere"; grep 'yame fetch' h.txt; exit 1; }
grep -q 'Resolved in order: -d, \$METHPROBE_DATA_HOME, \$YAME_DATA_HOME' h.txt ||
  { echo "the store order does not list the tool's variable first"; grep Resolved h.txt; exit 1; }
grep -q "METHPROBE_DATA_HOME: $METHPROBE_DATA_HOME" h.txt ||
  { echo "the resolved root is not the tool's store"; grep DATA_HOME h.txt; exit 1; }

## ---- 2. -l lists the tool's catalogue and nothing of yame's ----------------
./probe_fetch -l </dev/null > l.txt 2>/dev/null
[ "$(tail -n +2 l.txt | wc -l)" -eq 1 ] || { echo "-l listed $(tail -n +2 l.txt | wc -l) rows, want the tool's 1"; cat l.txt; exit 1; }
grep -q 'hg38/probe.*one.cm' l.txt || { echo "-l does not list the tool's file"; cat l.txt; exit 1; }
grep -q 'InfiniumAnnotation' l.txt && { echo "-l leaked yame's catalogue"; exit 1; }

## ---- 3. errors are in the tool's voice --------------------------------------
./probe_fetch nosuch </dev/null > e.txt 2>&1 || true
grep -q '^methprobe fetch: nothing in the catalogue is called nosuch' e.txt ||
  { echo "an unknown name was not refused as methprobe"; head -2 e.txt; exit 1; }
grep -q 'methprobe fetch -l' e.txt || { echo "the hint does not name the tool's -l"; exit 1; }

## ---- 4. the stale-store report is in the tool's voice too -------------------
## Stage the tool's directory with its file on disk, recorded at a digest its
## registry does not pin: stale, and the line says so with the tool's verb,
## naming the -y -f command that repairs it.
mkdir -p "$METHPROBE_DATA_HOME/hg38/probe"
: > "$METHPROBE_DATA_HOME/hg38/probe/one.cm"
printf 'ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff  one.cm\n' > "$METHPROBE_DATA_HOME/hg38/probe/SHA256SUMS"
./probe_fetch -l </dev/null 2> r.txt >/dev/null
grep -q '^\[methprobe fetch\] hg38/probe: 1 of 1 files come from an earlier release of zhou-lab/probe than this methprobe pins (one.cm); replace them with: methprobe fetch -y -f hg38/probe' r.txt ||
  { echo "the stale-store line is not in the tool's voice"; cat r.txt; exit 1; }
## and the tool's store variable wins over YAME_DATA_HOME: nothing was read from ystore
[ -z "$(find "$YAME_DATA_HOME" -type f)" ] || { echo "the tool read yame's store"; exit 1; }

## ---- 5. the library is self-contained ---------------------------------------
## A static archive links whole objects, so an object in libyame.a that
## references a symbol only the yame EXECUTABLE defines breaks every downstream
## that pulls that object -- v1.43's first cut did it through summary.o. Linking
## the WHOLE archive into a trivial program is the definitive test: every
## undefined reference in every object must resolve from the archive itself or
## its declared dependencies, or this fails at link time.
printf 'int main(void) { return 0; }\n' > whole.c
case "$(uname -s)" in
  Darwin) whole="-Wl,-all_load $root/libyame.a" ;;
  *)      whole="-Wl,--whole-archive $root/libyame.a -Wl,--no-whole-archive" ;;
esac
libs=$("$root/yame-config" --libs | sed "s|$root/libyame.a||")
${CC:-cc} -std=gnu99 $("$root/yame-config" --cflags) -o whole whole.c $whole $libs 2>whole.err ||
  { echo "libyame.a references a symbol it does not provide (an object depends on the executable):"
    grep -E 'undefined reference|Undefined symbols|referenced from' whole.err | head -5; exit 1; }
