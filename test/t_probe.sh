#!/bin/bash
## Layer 4: build a small consumer against libyame.a and check the library
## contract it depends on. Skips (passing) if the archive is absent, since a
## plain `make` builds only the CLI.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/.." && pwd)
cfg="$root/yame-config"
lib="$root/libyame.a"

if [ ! -f "$lib" ] || [ ! -x "$cfg" ]; then
  echo "skip: no libyame.a (run 'make lib')" >&2
  exit 0
fi

d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
cd "$d"
printf '3\t1\n0\t0\n5\t5\n2\t8\n' > one.txt
"$YAME" pack -f m one.txt > one.cg
for i in 1 2 3; do awk -v v=$i 'BEGIN { for (r = 0; r < 4; r++) print v "\t" (9 - v) }'; done |
  split -l 4 - part.
for f in part.*; do "$YAME" pack -f m "$f" >> three.cg; done

## A CX stream that is a prefix of a larger file: one record, then bytes that
## are not a BGZF block. This is the methscope bundle shape.
cp one.cg bundle.bin
limit=$(wc -c < bundle.bin)
printf 'MSBNDL1\0\3\0\0\0mrmp' >> bundle.bin

cc -O1 -g -std=gnu99 $("$cfg" --cflags) -o probe "$here/probe.c" $("$cfg" --libs) 2>cc.err ||
  { echo "probe did not build"; cat cc.err; exit 1; }
./probe one.cg three.cg bundle.bin "$limit"
