---
title: yame pack
parent: YAME Subcommands
nav_order: 1
---

# yame pack
```bash

Usage: yame pack [options] <in.txt> <out.cx>
The input text file must match the dimension and order of the reference CpG bed file.

Options:
    -f [char] Format specification (choose one character or number):
              (b) Binary data. Format default to 0 or 1 depending on size:
                  0 - 1 byte for 8 binary CpGs
                  1 - Value (1 byte) + Run-Length Encoding (RLE) (2 bytes)
              (s) State data. Format default to 2:
                  2 - State text + Index RLE (Best for chromatin states).
                      Use format 0 for sequence context.
              (m) Sequencing MU data. Format default to 3:
                  3 - MU RLE + Ladder byte (Input: 2-column text, M and U).
              (d) Differential meth data. Format default to 6:
                  5 - 2-bits + NA-RLE (Input: only 0, 1, 2 values).
                  6 - 2-bits boolean for S (set) and U (universe).
                      (Input: 2-column text, S and U).
              (n) Fraction data. Format default to 4:
                  4 - Fraction / NA-RLE (32 bytes).
              (r) Reference coordinates. Format default to 7:
                  7 - Compressed BED format for CGs.
    -u [int]  Number of bytes per unit data when inflated (1-8).
              Lower values are more memory efficient but may be lossier.
              0 - Inferred from data.
    -v        Verbose mode
    -h        Display this help message
```
