---
title: yame rowop
parent: YAME Subcommands
nav_order: 11
---

# yame rowop
```bash

Usage: yame rowop [options] <in.cx> <out>
Options:
    -o        Operations (choose one):
              binasum     Sum binary data to M and U (format 3).
                          Output: new cx file.
              musum       Sum M and U separately (format 3).
                          Output: new cx file.
              mean        Mean beta and counts of data points (format 3).
                          Output: plain text (two columns).
              std         Standard deviation. Requires format 3 cx.
                          Output: plain text (std, counts).
              binstring   Binarize data to row-wise string (format 3).
                          Output: plain text file with binary strings.
              cometh      Co-methylation of neighboring CGs.
                          Output: plain text in uint64_t U0U1-U0M1-M0U1-M0M1,
                          U0U2-U0M2-M0U2-M0M2, etc. '0' is target CG,
                          followed by 1, 2, etc for neighboring CGs.
                          Each pair occupies 16 bits. For visual, use -v.
                          Intermediate methylations (0.3-0.7) are excluded.
    -w        Number of neighboring CGs for cometh (default: 5).
    -c        Minimum sequencing depth for rowops (default 1).
    -v        Verbose mode
    -h        Display this help message


```
