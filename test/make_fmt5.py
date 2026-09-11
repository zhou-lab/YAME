import struct, sys, zlib
#!/usr/bin/env python3
"""Write a format-5 CX record.

Format 5 is obsolete: yame still DECODES it, and nothing can pack one any
more, so the only way to test the decoder is to write the bytes by hand.
This encodes a ternary vector (0 off, 1 on, 2 NA) exactly as format5.c
documents -- NA runs as a byte with the top bit clear, 0/1 values packed four
to a byte with the top bit set -- wraps it in a CX record header, and puts
that in a BGZF member with the trailing empty member every store ends with.

    make_fmt5.py <out.cg> [values...]      default: a mixed vector
"""

## Build a BGZF member the crude but valid way: a gzip member carrying the
## BC extra field that marks its own size.
def bgzf_block(payload):
    """One BGZF member: a gzip member carrying the BC extra field that states
    its own compressed size, which is what makes it seekable."""
    c = zlib.compressobj(6, zlib.DEFLATED, -15)
    co = c.compress(payload) + c.flush()
    xlen = 6
    bsize = 12 + xlen + len(co) + 8
    hdr = struct.pack("<BBBBIBBHBBHH", 31, 139, 8, 4, 0, 0, 255, xlen,
                      66, 67, 2, bsize - 1)
    return hdr + co + struct.pack("<II", zlib.crc32(payload) & 0xffffffff,
                                  len(payload) & 0xffffffff)
CDSIG = 266563789635
BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")

if len(sys.argv) > 2:
    vals = [int(a) for a in sys.argv[2:]]
else:
    ## leading NA run, packed 0/1 groups that fill a byte and that do not,
    ## an NA run in the middle, and a tail that ends mid-byte
    vals = [2, 2, 2, 0, 1, 0, 2, 1, 1, 0, 2, 2, 0, 1, 1, 1]
## encode: NA runs as a byte with MSB 0, 0/1 runs packed 4 per byte MSB 1
out = bytearray(); i = 0
while i < len(vals):
    if vals[i] == 2:
        run = 0
        while i < len(vals) and vals[i] == 2 and run < 127: run += 1; i += 1
        out.append(run)
    else:
        slots = []
        while i < len(vals) and vals[i] != 2 and len(slots) < 4:
            slots.append(vals[i]); i += 1
        b = 0x80
        for k, v in enumerate(slots):
            off = 6 - 2*k
            b |= (0x2 << off) | (v << off)
        out.append(b)
payload = struct.pack("<Q", CDSIG) + b"5" + struct.pack("<Q", len(out)) + bytes(out)
open(sys.argv[1], "wb").write(bgzf_block(payload) + BGZF_EOF)
print(" ".join("NA" if v == 2 else str(v) for v in vals))
