---
title: yame subset
parent: YAME Subcommands
nav_order: 3
---

# yame subset
```bash

Usage: yame subset [options] <in.cx> [<sample1> <sample2> ...]
If -o <out.cx>, an index will also be generated. Otherwise, output .cx to stdout without index.

Options:
    -v        verbose
    -o        output cx file name. if missing, output to stdout without index.
    -l        Path to the sample list. Ignored if sample names are provided on the command line.
    -s        Filter format 2 <in.cx> instead of samples in files. Output format 0.
    -H [N]    Process N samples from the start of the list, where N is less than or equal to the total number of samples.
    -T [N]    Process N samples from the end of the list, where N is less than or equal to the total number of samples. Requires index.
    -h        This help

```

