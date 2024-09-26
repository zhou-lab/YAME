---
title: Methylation Data Storage
nav_order: 2
---

# Methylation data pack and unpack

## yame pack

`yame papck` provides the functionality of packaging different inputs into `.cx` file for easier downstream analysis.  
Note: please make sure that the input file match the dimension and order of the reference bed file. You can use [bedtools intersect](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html) to match with the reference bed file.

`yame pack` has the following format specification `-f`

For binary data:
    0. 1 byte for 8 binary CpGs
    1. Value (1 byte) + Run-Length Encoding (RLE) (2 bytes)

For state data:
    2. State text + Index RLE (Best for chromatin states)

For sequencing data:
    3. MU RLE + Ladder byte (Input: 2-column text, M and U)

For fraction data:
    4. Fraction / NA-RLE (32 bytes)

For differential meth data:
    5. 2-bits + NA-RLE (Input: only 0, 1, 2 values)
    6. 2-bits boolean for S (set) and U (universe)

For referene coordinates:
    7. Compressed BED format for CGs

For more help with `pack`, run `yame pack` in the terminal or check out the
[pack help page]({% link docs/subcommands/YAME_pack.markdown %}).

## yame unpack

`yame unpack` provides the decoding functionality from `yame pack`.

Example usage:

```bash
yame unpack -a -f 1 yourfile.cg
```

This command will unpack the .cg files for all the samples it contains, and output the fraction of methylation/coverage with coverage >= 1.

For more help with `unpack`, run `yame unpack` in the terminal or check out the
[unpack help page]({% link docs/subcommands/YAME_unpack.markdown %}).