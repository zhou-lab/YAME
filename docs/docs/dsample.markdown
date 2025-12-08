---
title: Downsampling
nav_order: 8
---

# Downsampling of methylation sites

Sometimes, we want to test computational methods across different sparisity levels, this can be achieved with `yame dsample` function.

Below shows an example where we only want to keep 10,000 CG sites under a fixed seed.

```bash
yame dsample -s 1 -N 10000 yourfile.cg > output.cg
``` 

For more help with `dsample`, run `yame dsample` in the terminal or check out the
[dsample help page]({% link docs/subcommands/YAME_dsample.markdown %}).
