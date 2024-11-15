---
title: yame rowsub
parent: YAME Subcommands
nav_order: 4
---

# yame rowsub
```bash

Usage: yame rowsub [options] <in.cx>
This function outputs to stdout.
The 0 in [beg0] below means 0-based. Similarly, [beg1], [end1], [index1], etc.
The number in (), e.g., [blockIndex0]_(blockSize), is optional with a default.

Options:
    -v        verbose
    -l [PATH] rows in a plain text of [index1] on each row. index1: 1-based. No sorting requirement.
    -L [PATH] rows in a plain text of [chrm]_[beg1] on each row. Requires -R. No sorting requirement.
    -R [PATH] row coordinates to use. Required by -L.
    -1        The row coordinate (from -R) will be added to output as the first dataset.
    -m [PATH] rows in a mask file (format 1 or 2).
    -B [STR]  a row index range [rowIndexBeg0]_(rowIndexEnd1). By default, rowIndexEnd1=rowIndexBeg0+1.
    -I [STR]  a row index range [blockIndex0]_(blockSize). By default, blockSize=1000000.
    -h        This help

```
