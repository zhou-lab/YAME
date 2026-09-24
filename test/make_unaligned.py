#!/usr/bin/env python3
"""Join CX files into one whose records do NOT start on block boundaries.

Every yame writer flushes before each record, so its records are block-aligned
and `subset` / `split` can copy their compressed bytes verbatim. A store built
by appending to one shared writer is not, and those commands must notice and
fall back to decode/encode. Nothing yame writes today has that shape, so this
makes one: it inflates the inputs and deflates all their records into a
single BGZF member, followed by the empty member every store ends with.

    make_unaligned.py <out.cx> <in1.cx> [in2.cx ...]
"""
import gzip, struct, sys, zlib

out, ins = sys.argv[1], sys.argv[2:]
payload = b"".join(gzip.open(f).read() for f in ins)
c = zlib.compressobj(6, zlib.DEFLATED, -15)
co = c.compress(payload) + c.flush()
## BGZF: a gzip member whose BC extra field states the member's own size
bsize = 12 + 6 + len(co) + 8
hdr = struct.pack("<BBBBIBBHBBHH", 31, 139, 8, 4, 0, 0, 255, 6, 66, 67, 2, bsize - 1)
eof = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
with open(out, "wb") as fh:
    fh.write(hdr + co + struct.pack("<II", zlib.crc32(payload) & 0xffffffff,
                                    len(payload) & 0xffffffff) + eof)
