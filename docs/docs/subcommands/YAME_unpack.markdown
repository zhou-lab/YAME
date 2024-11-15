---
title: yame unpack
parent: YAME Subcommands
nav_order: 2
---

# yame unpack
```bash

Usage: yame unpack [options] <in.cx> [[sample 1], [sample 2], ...]

Options:
    -a        Process all samples
    -C        Output column names
    -R [PATH] Row coordinate .cr file name.
    -r        0: Row coordinate output in chrm-beg0-end1 (default, for cg).
              1: Row coordinate output in chrm-beg0-end0 (for allc).
              other: Row coordinate output in chrm_beg1.
    -l [PATH] Path to the sample list. Ignored if sample names are provided on the command line.
    -H [N]    Process N samples from the start of the list, where N is less than or equal to the
              total number of samples.
    -T [N]    Process N samples from the end of the list, where N is less than or equal to the
              total number of samples. Requires index.
    -f [N]    Display format for data format 3. Options are:
                   N == 0: Compound MU
                   N <  0: M<tab>U
                   N >  0: Fraction (with number for the min coverage)
    -c        Enable chunk process
    -s        Chunk size (default is 1M)
    -u [int]  number of bytes for each unit data while inflated. Lower number needs less memory
              efficient but could be lossier. Can only be 1-8.
              0 means this will be inferred from data.
    -h        Display this help message.

```
