#!/bin/bash
## Adversarial input. The contract is one sentence:
##
##   every input either yields the right answer, or exits non-zero with a
##   message -- exit 0 with wrong data is the only forbidden outcome.
##
## Asserted as an invariant rather than as a golden map, so it holds for any
## fixture: a truncated store either fails, or reads as a PREFIX of the
## records the whole store holds. v1.40 hardened this path after junk, a
## truncated .cg and an empty file all read as an empty store with rc 0.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
cd "$d"

## ---- a 4-record store, each record from its own pack so records sit on
## ---- block boundaries the way a concatenated store's do
for i in 1 2 3 4; do
  awk -v v=$i 'BEGIN { for (r = 0; r < 25; r++) print v "\t" (10 - v) }' > s$i.txt
  "$YAME" pack -f m s$i.txt > s$i.cg
done
cat s1.cg s2.cg s3.cg s4.cg > whole.cg
printf 's1\ns2\ns3\ns4\n' > names.txt
"$YAME" index -s names.txt whole.cg
## Nrow, Format, UnitBytes -- the structural columns. Sample and NSample are
## deliberately left out: a truncated copy has no .idx beside it, so those two
## read NA there and would differ for a reason that is not corruption.
"$YAME" info whole.cg | tail -n +2 | cut -f4-6 > info.whole
[ "$(wc -l < info.whole)" -eq 4 ] || { echo "fixture: expected 4 records"; exit 1; }
size=$(wc -c < whole.cg)

## ---- 1. truncation at every interior byte ----------------------------------
## Either the command fails, or what it reports is a prefix of the whole
## store's records. Nothing may report a record the file cannot hold.
bad=0; clean=0; failed=0
for ((n = 1; n < size; n++)); do
  head -c "$n" whole.cg > cut.cg
  if out=$("$YAME" info cut.cg 2>/dev/null); then
    clean=$((clean + 1))
    printf '%s\n' "$out" | tail -n +2 | cut -f4-6 > info.cut
    k=$(wc -l < info.cut)
    if [ "$k" -gt 0 ] && ! head -n "$k" info.whole | diff -q - info.cut >/dev/null; then
      echo "truncation at $n: exit 0 but records are not a prefix of the whole store"
      diff <(head -n "$k" info.whole) info.cut | head -4
      bad=$((bad + 1))
    fi
  else
    failed=$((failed + 1))
  fi
done
[ "$bad" -eq 0 ] || { echo "$bad truncation(s) returned wrong data with exit 0"; exit 1; }
## Both outcomes must actually occur, or the loop is not exercising the path.
[ "$failed" -gt 0 ] || { echo "no truncation was detected at all"; exit 1; }
[ "$clean" -gt 0 ] || { echo "no truncation read as a clean prefix"; exit 1; }

## ---- 2. inputs that are not a CX stream at all -----------------------------
## Each must exit non-zero and say something; none may print a record row.
notcx() {                       # <label> <file>
  local out rc
  out=$("$YAME" info "$2" 2>err.txt) && rc=0 || rc=$?
  [ "$rc" -ne 0 ] || { echo "$1: expected non-zero exit"; echo "$out"; exit 1; }
  [ -s err.txt ] || { echo "$1: failed silently"; exit 1; }
}
head -c 200 /dev/urandom > noise.bin
notcx "random bytes" noise.bin
printf 'chr1\t1\t2\n' > text.txt
notcx "a text file" text.txt
## gzip that is not BGZF: a real gzip member with no BC extra field
printf 'hello\n' | gzip > plain.gz
notcx "plain gzip" plain.gz

## An EMPTY file is not corrupt: zero records, exit 0. Documented, and the
## v1.40 hardening deliberately kept it.
: > empty.cg
"$YAME" info empty.cg >/dev/null 2>&1 || { echo "empty file should read as zero records"; exit 1; }

## ---- 3. a store followed by bytes that are not a BGZF block ----------------
## The bundle shape (methscope MSBNDL1). Unbounded this must fail rather than
## silently stop -- the bound is opt-in, through the library (see t_probe.sh).
cat whole.cg > tail.cg; printf 'MSBNDL1\0\3\0\0\0mrmp' >> tail.cg
notcx "store + foreign tail" tail.cg

## ---- 4. two end markers in a row ------------------------------------------
## A concatenated store normally carries one empty BGZF member between
## records; two in a row used to exhaust the reader's single-member tolerance
## and make every later record silently invisible.
"$YAME" pack -f m s1.txt > e1.cg      # ends with its own empty member
cat e1.cg e1.cg s2.cg > dbl.cg        # e1's trailer meets e1's trailer
n_dbl=$("$YAME" info dbl.cg 2>/dev/null | tail -n +2 | wc -l) || n_dbl=fail
[ "$n_dbl" = 3 ] || { echo "two end markers in a row: saw $n_dbl records, want 3"; exit 1; }

## ---- 5. a header that promises more than the file holds --------------------
## Cut inside the final record's data: the record header is intact and says
## how many bytes follow, so this is caught by the length check, not by BGZF.
head -c $((size - 40)) whole.cg > promise.cg
notcx "header promises more than follows" promise.cg

## ---- 6. an index that does not match the file ------------------------------
## A wrong .idx must not be believed into a bad seek. Offsets past end of
## file, and every offset identical (what conda's index -1 produced before
## v1.42), both have to be refused rather than returning the wrong record.
## Binary goes to a file, never through a shell variable: command
## substitution strips NUL bytes, which a .cg is full of.
idxbad() {                      # <label>
  local rc
  "$YAME" subset whole.cg s3 > got.cg 2>err.txt && rc=0 || rc=$?
  if [ "$rc" -eq 0 ]; then
    ## Exit 0 is allowed only if the bytes really are s3.
    "$YAME" unpack got.cg > got.txt 2>/dev/null ||
      { echo "$1: exit 0 but the output is not a readable store"; exit 1; }
    "$YAME" unpack s3.cg > want.txt
    diff -q want.txt got.txt >/dev/null ||
      { echo "$1: exit 0 but returned the wrong record"; exit 1; }
  else
    [ -s err.txt ] || { echo "$1: failed silently"; exit 1; }
  fi
}
cp whole.cg.idx idx.good
awk -F'\t' 'BEGIN {OFS = "\t"} {print $1, 999999999}' idx.good > whole.cg.idx
idxbad "index past end of file"
awk -F'\t' 'BEGIN {OFS = "\t"} {print $1, 0}' idx.good > whole.cg.idx
idxbad "index with every offset identical"
cp idx.good whole.cg.idx

## ---- 7. a failed allocation is loud ---------------------------------------
## An out-of-memory malloc returns NULL and the caller writes through it: on a
## busy machine that is a segfault with an empty stderr, and sometimes a
## TRUNCATED stream with exit 0, which downstream cannot tell from a genuinely
## short one. A `paste` pipeline scored two models on 1.93M and 189K rows that
## way and the asymmetry read as a real difference between them.
##
## Under a hard address-space cap the command must name what it could not
## allocate and exit non-zero. Skipped where `ulimit -v` does not bite.
if ( ulimit -v 20000 2>/dev/null ); then
  awk 'BEGIN { for (i = 0; i < 200000; i++) print (i % 10) "\t" (10 - (i % 10)) }' |
    "$YAME" pack -f m - > big.cg 2>/dev/null
  ## `&& rc=0 || rc=$?`: the subshell is EXPECTED to fail, and a bare `; rc=$?`
  ## would let `set -e` kill the script before the assignment ran.
  rc=0
  ( ulimit -v 20000; "$YAME" unpack -f 1 big.cg > oom.out 2> oom.err ) || rc=$?
  if [ "$rc" -eq 0 ]; then
    echo "a starved unpack exited 0 (wrote $(wc -l < oom.out) lines)"; exit 1
  fi
  [ "$rc" -lt 128 ] || { echo "a starved unpack died on signal $((rc - 128))"; exit 1; }
  ## which allocation fails first depends on the build (-O0 lays the heap out
  ## differently), so accept the shared phrase from any of them
  grep -qiE 'out of memory|cannot allocate' oom.err ||
    { echo "a starved unpack did not say it ran out of memory"; head -3 oom.err; exit 1; }
fi
