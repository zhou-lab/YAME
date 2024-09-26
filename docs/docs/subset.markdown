---
title: Methylation Data Subset
nav_order: 3
---

# Subset of the packed .cx file

## Subset samples

`yame subset` will subset the samples from the multi sample `.cx` file. We can store the samples we want to keep with a text file or simply providing them in the command line. Make sure the sample names provided is in the sample names of `.cx` file.

```bash
yame subset -l sample.txt yourfile.cg
```

For more help with `subset`, run `yame subset` in the terminal or check out the
[subset help page]({% link docs/subcommands/YAME_subset.markdown %}).

## Subset cx sites

If we want to obtain methylation sites from specific regions, we can use this `yame rowsub` function. 
If we have a list of CpG sites we want to subset, for example 

```
chr16_18300002
chr16_18300046
chr16_18300140
chr16_18300162
chr16_18300172
```
We can use `-L` with `-R` together (we provided two row coordiante .cr files [mm10](https://github.com/zhou-lab/KYCGKB_mm10)/[hg38](https://github.com/zhou-lab/KYCGKB_hg38)). 

```bash
yame rowsub -R cpg_nocontig.cr -L CpG_sites.tsv yourfile.cg > subset.cg
```

For more help with `rowsub`, run `yame rowsub` in the terminal or check out the
[rowsub help page]({% link docs/subcommands/YAME_rowsub.markdown %}).