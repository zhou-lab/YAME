#!/bin/bash
## Format 7, coordinates: delta-encoded positions, a varint per row.
##
## A delta up to 0x3fff fits in two bytes; anything larger takes eight. Every
## fixture in the suite spaced its rows 100 bp apart, so the 8-byte encoding
## -- which every real genome needs, at every gap over 16 kb -- was never
## written or read by a test, and neither was a stream cut off inside one.
set -euo pipefail
YAME=${YAME:?export YAME=/path/to/yame}
d=$(mktemp -d); trap 'rm -rf "$d"' EXIT
cd "$d"

## ---- 1. deltas of every width round-trip ----------------------------------
## 100 -> 20000 is 19900 (8 bytes); 20000 -> 5e9 is past 32 bits; on chr2,
## 10 -> 16393 is 16383, the largest 2-byte delta
printf 'chr1\t100\t102\nchr1\t20000\t20002\nchr1\t5000000000\t5000000002\nchr2\t10\t12\nchr2\t16393\t16395\n' > big.bed
"$YAME" pack -v -f r big.bed big.cr 2>v.txt
grep -q 'Format 7' v.txt || { echo "pack -v -f r said nothing of the format"; exit 1; }
diff big.bed <("$YAME" unpack big.cr 2>/dev/null) || { echo "large deltas did not round-trip"; exit 1; }
## and a slice of it decodes the same rows
diff <(sed -n 2,3p big.bed) <("$YAME" rowsub -B 1_3 big.cr 2>/dev/null | "$YAME" unpack - 2>/dev/null) ||
  { echo "rowsub across an 8-byte delta returned the wrong rows"; exit 1; }
if "$YAME" rowsub -B 9_10 big.cr >/dev/null 2>err.txt; then
  echo "rowsub -B past the end of a coordinate stream was accepted"; exit 1
fi
grep -q 'bigger than the data vector size' err.txt || { echo "-B past the end failed without saying why"; exit 1; }

## ---- 2. what pack refuses ----------------------------------------------------
if printf 'chr1\tabc\t3\n' | "$YAME" pack -f r - bad.cr 2>err.txt; then
  echo "pack -f r accepted a start that is not a number"; exit 1
fi
grep -q 'not a nonnegative integer' err.txt || { echo "a bad start failed without saying why"; exit 1; }
if printf 'chr1\n' | "$YAME" pack -f r - bad.cr 2>err.txt; then
  echo "pack -f r accepted a line with one field"; exit 1
fi
grep -q 'fewer than 2 fields' err.txt || { echo "a one-field line failed without saying why"; exit 1; }

## ---- 3. a stream that stops inside an 8-byte delta ---------------------------
## The record header must agree with the short body, or the reader stops at
## the header instead; this is corruption the header cannot see, and it must
## be named rather than read as a shorter coordinate list.
if command -v python3 >/dev/null; then
  printf 'chr1\t100\t102\nchr1\t5000000000\t5000000002\n' > two.bed
  "$YAME" pack -f r two.bed two.cr
  python3 - <<'PY'
import gzip, zlib, struct
p = gzip.open("two.cr").read()
sig, fmt, n = p[:8], p[8:9], struct.unpack("<Q", p[9:17])[0]
body = p[17:17 + n]
cut = body[:5 + 1 + 4]        # "chr1\0", a 1-byte delta, then 4 of the 8 bytes
rec = sig + fmt + struct.pack("<Q", len(cut)) + cut
c = zlib.compressobj(6, zlib.DEFLATED, -15); co = c.compress(rec) + c.flush()
bs = 12 + 6 + len(co) + 8
hdr = struct.pack("<BBBBIBBHBBHH", 31, 139, 8, 4, 0, 0, 255, 6, 66, 67, 2, bs - 1)
eof = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
open("mid.cr", "wb").write(hdr + co + struct.pack("<II", zlib.crc32(rec) & 0xffffffff, len(rec)) + eof)
PY
  if "$YAME" unpack mid.cr >/dev/null 2>err.txt; then
    echo "a stream cut inside an 8-byte delta was read as if whole"; exit 1
  fi
  grep -q 'ends mid-record: 8-byte delta' err.txt ||
    { echo "a truncated 8-byte delta failed without naming it"; cat err.txt; exit 1; }
fi

## ---- 4. a coordinate stream as a summary query: one row per chromosome ------
[ "$("$YAME" summary -T big.cr 2>/dev/null | tail -n +2 | cut -f2,6 | paste -sd' ' -)" = \
  "$(printf '1-chr1\t3 1-chr2\t2')" ] || { echo "-T on a coordinate query"; exit 1; }
echo "ok: t_format7"
