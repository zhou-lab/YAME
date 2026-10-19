#!/bin/bash
## rowsub on formats 0 and 6 without decompressing the record. Both are fixed
## width -- 1 and 2 bits a row, stored as is -- so row r sits at a known byte
## and rowsub reads only the BGZF blocks holding the rows asked for, through
## the file's block map; no index column is needed. The bytes out must be
## the ones a full decode gives, which the pipe form always does.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
cd "$d"

## row counts that are not a multiple of 8 or 4, so the last byte is partial
awk 'BEGIN { srand(1); for (i = 0; i < 1200003; i++) print (rand() < 0.4) }' |
  "$YAME" pack -f b - > f0.cg
awk 'BEGIN { srand(2); for (i = 0; i < 1200001; i++) {
  u = (rand() < 0.7); print (u && rand() < 0.5) "\t" u } }' |
  "$YAME" pack -f d - > f6.cg
awk 'BEGIN { srand(3); for (i = 0; i < 1200001; i++) print int(rand() * 9) "\t" int(rand() * 9) }' |
  "$YAME" pack -f m - > f3.cg
cat f0.cg f6.cg f3.cg > s.cg                  # f0 is the longest record
printf 'f0\nf6\nf3\n' > names.txt
"$YAME" index -s names.txt s.cg

## the full read: the same bytes through a pipe, which cannot seek
full() { cat s.cg | "$YAME" rowsub "$@" -; }
path() {   # path <want> <rowsub args...>
  local want=$1; shift
  local got
  got=$(YAME_ROWSUB_PATH=1 "$YAME" rowsub "$@" 2>&1 >/dev/null | grep '^\[rowsub\] row-index')
  [ "$got" = "[rowsub] $want" ] ||
    { echo "rowsub $*: expected '$want', got '$got'"; exit 1; }
}
same() {   # same <label> <rowsub args without the input>
  local l=$1; shift
  full "$@" > want.cg
  "$YAME" rowsub "$@" s.cg > got.cg
  cmp want.cg got.cg || { echo "$l: differs from a full decode"; exit 1; }
}

## rows on and around byte boundaries (8 rows a byte for f0, 4 for f6), the
## first and last rows, repeats, out of order, and a random scatter
awk 'BEGIN { srand(5); for (i = 0; i < 400; i++) print int(rand() * 1200001) + 1
  print 1; print 1200001; print 8; print 9; print 4; print 5; print 65536 * 4
  print 65536 * 4 + 1; print 7; print 7; print 1200000 }' > rows.txt
same "-l" -l rows.txt
path "row-index 1, full 0, fixed-width 2" -l rows.txt s.cg
for sel in "-B 0_1" "-B 1199990_1200001" "-B 100003_950007" "-I 11_100000" "-B 1199999"; do
  same "$sel" $sel
  path "row-index 1, full 0, fixed-width 2" $sel s.cg
done
## -S: chosen records in the order asked
printf 'f6\nf3\nf0\nf6\n' > pick.txt
"$YAME" rowsub -S pick.txt -l rows.txt s.cg > got.cg
## the pipe has no -S; the reference is subset | rowsub
"$YAME" subset -l pick.txt s.cg | "$YAME" rowsub -l rows.txt - > want.cg
cmp want.cg got.cg || { echo "-S: differs from subset | rowsub"; exit 1; }
path "row-index 1, full 0, fixed-width 3" -S pick.txt -l rows.txt s.cg

## a row past the end of the shorter records: the same refusal as a full read
printf '1200003\n' > high.txt
if "$YAME" rowsub -l high.txt s.cg > /dev/null 2> e1.txt; then
  echo "a row past the end was accepted"; exit 1; fi
full -l high.txt > /dev/null 2> e2.txt || true
diff <(sed 's/s.cg/X/' e1.txt) <(sed 's/s.cg/X/' e2.txt) ||
  { echo "out-of-range: a different error from a full read"; exit 1; }
## a mask reads every row, so it takes the full path
awk 'BEGIN { for (i = 0; i < 1200003; i++) print (i % 7 == 0) }' | "$YAME" pack -f b - > m.cg
cat f0.cg > one.cg; "$YAME" index one.cg
path "row-index 0, full 1, fixed-width 0" -m m.cg one.cg

## ---- it reads only what it needs --------------------------------------------
## Damage a BGZF block in the middle of f6. A full decode must fail on it;
## asked for rows in f6's first bytes, the fixed-width path never reads it.
cp s.cg dmg.cg; cp s.cg.idx dmg.cg.idx
python3 - dmg.cg f0.cg <<'PY'
import os, struct, sys
f = sys.argv[1]; b = bytearray(open(f, "rb").read())
start = os.path.getsize(sys.argv[2])          # f6 begins where f0 ends
offs, o = [], start
while o < len(b):
    offs.append(o); o += struct.unpack("<H", b[o+16:o+18])[0] + 1
mid = offs[2]                                 # f6's third block
for i in range(mid + 40, mid + 80): b[i] ^= 0xFF
open(f, "wb").write(b)
PY
printf 'f6\n' > f6.txt; printf '1\n2\n17\n' > early.txt
"$YAME" rowsub -S f6.txt -l early.txt dmg.cg > got.cg ||
  { echo "fixed-width path read a block it did not need"; exit 1; }
"$YAME" rowsub -S f6.txt -l early.txt s.cg > want.cg
cmp want.cg got.cg || { echo "damaged file: rows differ"; exit 1; }
if cat dmg.cg | "$YAME" rowsub -l early.txt - > /dev/null 2>&1; then
  echo "fixture: the damage did not break a full read"; exit 1; fi

## ---- records starting and ending mid-block ------------------------------------
python3 - s.cg mid.cg <<'PY'
import gzip, struct, sys, zlib
raw = gzip.open(sys.argv[1], "rb").read()
def blk(b):
    c = zlib.compressobj(6, zlib.DEFLATED, -15); co = c.compress(b) + c.flush()
    return struct.pack("<BBBBIBBHBBHH", 31, 139, 8, 4, 0, 0, 255, 6, 66, 67, 2,
                       len(co) + 25) + co + struct.pack("<II", zlib.crc32(b), len(b))
with open(sys.argv[2], "wb") as f:
    for i in range(0, len(raw), 7777): f.write(blk(raw[i:i + 7777]))
    f.write(bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000"))
PY
"$YAME" index -s names.txt mid.cg
for sel in "-l rows.txt" "-B 100003_950007"; do
  cat mid.cg | "$YAME" rowsub $sel - > want.cg
  "$YAME" rowsub $sel mid.cg > got.cg
  cmp want.cg got.cg || { echo "mid-block $sel: differs from a full read"; exit 1; }
  path "row-index 1, full 0, fixed-width 2" $sel mid.cg
done

## ---- the block map in the .idx --------------------------------------------
## `yame index` writes formats 0 and 6 a third column too, under "f1.": the
## first byte of each block, so rowsub never reads a block header to find a
## row. (Cold on network storage a header costs about what reading and
## inflating its block does.) Without it, rowsub walks only the record's own
## block headers.
awk -F'\t' '$1 == "f0" || $1 == "f6"' s.cg.idx | awk -F'\t' '$3 !~ /^f1\.[0-9a-z._]+\.k[0-9a-z]+$/ {bad = 1} END {exit bad}' ||
  { echo "index: no f1. column on f0/f6"; cut -c1-60 s.cg.idx; exit 1; }
## Break the header of a block in the middle of f6. With the column nothing
## reads that header; with two columns the walk over f6's blocks stops at
## it and the record is read in full, which fails there.
cp s.cg hdr.cg; cp s.cg.idx hdr.cg.idx
python3 - hdr.cg f0.cg <<'PY'
import os, struct, sys
f = sys.argv[1]; b = bytearray(open(f, "rb").read())
o = os.path.getsize(sys.argv[2]); k = 0
while k < 3: o += struct.unpack("<H", b[o+16:o+18])[0] + 1; k += 1
b[o] ^= 0xFF                                  # f6's fourth block: bad magic
open(f, "wb").write(b)
PY
"$YAME" rowsub -S f6.txt -l early.txt hdr.cg > got.cg ||
  { echo "with the column, rowsub read a block header it did not need"; exit 1; }
"$YAME" rowsub -S f6.txt -l early.txt s.cg > want.cg
cmp want.cg got.cg || { echo "header-damaged file: rows differ"; exit 1; }
path "row-index 0, full 0, fixed-width 1" -S f6.txt -l early.txt hdr.cg
cut -f1,2 s.cg.idx > hdr.cg.idx
if "$YAME" rowsub -S f6.txt -l early.txt hdr.cg > /dev/null 2>&1; then
  echo "fixture: without the column the damaged header went unread"; exit 1; fi

## a stale column (the file rewritten behind it) and a column re-signed with
## a wrong block address are refused and the record located from its headers
cat f0.cg f6.cg f3.cg > re.cg; cp s.cg.idx re.cg.idx
awk 'BEGIN { srand(9); for (i = 0; i < 1200001; i++) { u = (rand() < 0.5); print (u && rand() < 0.5) "\t" u } }' |
  "$YAME" pack -f d - > f6b.cg
cmp -s <(stat -c %s f6.cg 2>/dev/null || stat -f %z f6.cg) <(stat -c %s f6b.cg 2>/dev/null || stat -f %z f6b.cg) && echo "  note: same size"
cat f0.cg f6b.cg f3.cg > re.cg
cat re.cg | "$YAME" rowsub -l rows.txt - > want.cg
"$YAME" rowsub -l rows.txt re.cg > got.cg 2> err.txt
cmp want.cg got.cg || { echo "stale f1. column: rows differ from a full read"; exit 1; }
python3 - s.cg.idx bad.idx <<'PY'
import sys, zlib
D = "0123456789abcdefghijklmnopqrstuvwxyz"
def e36(n):
    o = ""
    while True:
        o = D[n % 36] + o; n //= 36
        if not n: return o
out = []
for line in open(sys.argv[1]):
    f = line.rstrip("\n").split("\t")
    if f[0] == "f6":
        p = f[2].split(".")[:-1]
        ck = [[int(x, 36) for x in c.split("_")] for c in p[4:]]
        ck[1][3] += 7                            # a block's in-block offset
        p[4:] = ["_".join(e36(x) for x in c) for c in ck]
        body = ".".join(p)
        f[2] = body + ".k" + e36(zlib.crc32(body.encode()) & 0xffffffff)
    out.append("\t".join(f) + "\n")
open(sys.argv[2], "w").write("".join(out))
PY
## rows inside that block only, so the read seeks to the edited address
awk -F'\t' '$1 == "f6" {print $3}' s.cg.idx | python3 -c '
import sys; p = sys.stdin.read().strip().split(".")[4:-1]
sp = [0]
for c in p[1:]: sp.append(sp[-1] + int(c.split("_")[1], 36))
print("\n".join(str(4 * (sp[1] + j) + 1) for j in (5, 9, 300)))' > rows_k.txt
"$YAME" rowsub -S f6.txt -l rows_k.txt s.cg > want.cg
cp bad.idx s.cg.idx
"$YAME" rowsub -S f6.txt -l rows_k.txt s.cg > got.cg 2> err.txt
cmp want.cg got.cg || { echo "re-signed f1. column: rows differ"; exit 1; }
grep -q 'no longer match' err.txt || { echo "re-signed f1. column: not called stale"; exit 1; }
"$YAME" index -s names.txt s.cg

## writers carry it: subset -o copies records raw, index -1 keeps listed ones
printf 'f6\nf3\nf0\n' > pick3.txt
"$YAME" subset -o sub.cg -l pick3.txt s.cg > /dev/null
[ "$(awk -F'\t' '$3 ~ /^f1\./' sub.cg.idx | wc -l)" -eq 2 ] ||
  { echo "subset -o dropped an f1. column"; cut -c1-50 sub.cg.idx; exit 1; }
path "row-index 1, full 0, fixed-width 2" -l rows.txt sub.cg
cp f6.cg ap.cg; "$YAME" index -1 one ap.cg; cat f0.cg >> ap.cg; "$YAME" index -1 two ap.cg
[ "$(awk -F'\t' '$3 ~ /^f1\./' ap.cg.idx | wc -l)" -eq 2 ] ||
  { echo "index -1 lost an f1. column"; exit 1; }
cat ap.cg | "$YAME" rowsub -B 0_1200001 - > want.cg
"$YAME" rowsub -B 0_1200001 ap.cg > got.cg
cmp want.cg got.cg || { echo "index -1 file: rows differ"; exit 1; }

echo "ok: rowsub on formats 0 and 6 without a full decode"
