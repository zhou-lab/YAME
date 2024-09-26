---
title: yame summary
parent: YAME Subcommands
nav_order: 6
---

# yame summary
```bash

Usage: yame summary [options] <query.cm>
Query should be of format 0,1,2,3, can be a multi-sample set.

Options:
    -m        Mask feature (.cx) file, can be multi-sample.
              If '-', the whole sample will bed kept in memory, same as -M.
    -M        All masks will be loaded to memory. This save disk IO.
    -u        Optional universe set as a .cx file. If given, the masks and queries are both subset.
    -H        Suppress header printing.
    -q        The backup query file name if the query file name is '-'.
    -F        Use full feature/query file name instead of base name.
    -T        State features always show section names.
    -s        Sample list provided to override the query index file. Only applies to the first query.
    -h        This help.


```
