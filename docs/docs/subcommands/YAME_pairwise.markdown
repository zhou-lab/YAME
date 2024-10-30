---
title: yame pairwise
parent: YAME Subcommands
nav_order: 14
---

# yame pairwise
```bash

Usage: yame pairwise [options] <MU1.cx> (<MU2.cx>)
Return a format 6 set that represent differential methylation between MU1 and MU2.
If MU2 is not given, use the top 2 samples in MU1.cx.

Options:
    -o        output cx file name. if missing, output to stdout without index.
    -H        1: higher meth level in sample 1 than 2 (default).
              2: higher meth level in sample 2 than 1.
              others: diff levels, i.e., 1 and 2 combined.
    -c        minimum coverage (default: 1)
    -d        minimum delta meth level/effect size (default: 0)
    -h        This help

```
