#!/bin/bash
## Memory errors, locally. The CI matrix runs the suite under
## AddressSanitizer, which is stricter and faster, but ASan needs a runtime
## that is not installed everywhere; valgrind usually is, so this is the check
## a developer actually gets before pushing.
##
## Skips (passing) when valgrind is absent, and when the binary was built with
## a sanitizer -- the two do not coexist.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}

command -v valgrind >/dev/null || { echo "skip: no valgrind" >&2; exit 0; }
if strings "$YAME" 2>/dev/null | grep -q '__asan_'; then
  echo "skip: binary is ASan-instrumented; CI covers that leg" >&2; exit 0
fi

d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
cd "$d"

awk 'BEGIN { for (i = 0; i < 64; i++) print (i % 7) "\t" (7 - (i % 7)) }' > mu.txt
"$YAME" pack -f m mu.txt > a.cg
awk 'BEGIN { for (i = 0; i < 64; i++) print ((i + 3) % 7) "\t" (7 - ((i + 3) % 7)) }' |
  "$YAME" pack -f m - > b.cg
cat a.cg b.cg > two.cg
printf 'a\nb\n' > names.txt
"$YAME" index -s names.txt two.cg
awk 'BEGIN { for (i = 0; i < 64; i++) print (i % 4 == 0) ? 1 : 0 }' | "$YAME" pack -f b - > m.cg

## 64 rows is deliberate: dsample sized its bitset with a ceiling and memset
## it with floor-plus-one, so it wrote one byte past the heap on exactly the
## inputs whose row count divides by 8.
## Memory ERRORS only -- invalid reads and writes, uninitialised branches,
## double frees. Leaks at exit are deliberately not failed on: several
## subcommands end with a few bytes outstanding because the process is about
## to die anyway, and failing on those would bury the class that actually
## corrupts data. The CI sanitizer leg runs with detect_leaks=0 for the same
## reason.
run() {
  valgrind -q --error-exitcode=9 --leak-check=no \
           "$YAME" "$@" >/dev/null 2>vg.txt || {
    echo "valgrind flagged: yame $*"
    head -20 vg.txt
    exit 1
  }
}
run info two.cg
run unpack -a -f 1 two.cg
run unpack -c -s 8 -f 1 a.cg
run subset two.cg b
run rowsub -m m.cg a.cg
run rowsub -B 8_24 a.cg
run binarize a.cg
run mask a.cg m.cg
run dsample -s 1 -N 16 a.cg
run rowop -o binasum two.cg
run rowop -o stat two.cg
run summary -m m.cg a.cg
run pairwise a.cg b.cg
run pairwise -S a.cg b.cg
run pairwise -S -1 a two.cg two.cg
run hprint -c -g a.cg
