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

## Example of generating a feature file using yame pack

Here's an example showing how to generate the chrommHMM full stack hg38 annotation feature file. First, we have an input bed file `hg38_genome_100_segments.bed.gz` that might look like this:

```
chr1    10000   10400   2_GapArtf2
chr1    10400   10600   27_Acet1
chr1    10600   10800   38_EnhWk4
chr1    10800   12800   1_GapArtf1
chr1    12800   13000   38_EnhWk4
chr1    13000   13200   37_EnhWk3
chr1    13200   14800   6_Quies3
chr1    14800   15200   72_TxWk2
chr1    15200   16000   6_Quies3
chr1    16000   16200   83_TxEx3
chr1    16200   16400   82_TxEx2
```

We can run the following commands in the terminal, where we use [bedtools intersect](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html) with the reference CpG .cr file [mm10](https://github.com/zhou-lab/KYCGKB_mm10)/[hg38](https://github.com/zhou-lab/KYCGKB_hg38) to obtain the feature each CpG belongs to. 

```bash
yame unpack cpg_nocontig.cr | gzip > cpg_nocontig.bed.gz

zcat hg38_genome_100_segments.bed.gz | sortbed | bedtools intersect -a cpg_nocontig.bed.gz -b - -loj -sorted | bedtools groupby -g 1-3 -c 7 -o first -i - | cut -f4 | yame pack -f s - ChromHMMfullStack.cm
```

The output .cm feature file can be used to run enrichment and obtain aggregated methylation levels over different features, see [`enrichment`]({% link docs/enrichment.markdown %}). 

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