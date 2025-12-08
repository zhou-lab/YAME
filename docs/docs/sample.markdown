---
title: Sample Combination/Indexing
nav_order: 7
---

# Sample operations

`yame index` generates the index file for the `.cx` input. This is helpful when `.cx` contains multiple samples.

## Example of merging a feature file with multiple .cm files

Here's an example showing how to merge multiple histone modification ChIP-Seq peak bed files into one merged `.cm` feature file using `yame index`.
First, obtain a `.tsv` file describing the index and path for each of the bed file like `controlfiles.tsv` below:

```
268     GSM648494       human_hm/268_sort_peaks.narrowPeak.bed
269     GSM648495       human_hm/269_sort_peaks.narrowPeak.bed
272     GSM575295       human_hm/272_b_sort_peaks.broadPeak.bed
273     GSM575280       human_hm/273_sort_peaks.narrowPeak.bed
274     GSM575296       human_hm/274_b_sort_peaks.broadPeak.bed
275     GSM575281       human_hm/275_sort_peaks.narrowPeak.bed
367     GSM575223       human_hm/367_sort_peaks.narrowPeak.bed
368     GSM575222       human_hm/368_sort_peaks.narrowPeak.bed
382     GSM610328       human_hm/382_sort_peaks.narrowPeak.bed
```

Then, make individual `.cm` files for each feature, see [`yame pack`]({% link docs/pack_unpack.markdown %}).

```bash
cat controlfiles.tsv | parallel --colsep '\t' -j 72 'id={1};path={3}; sortbed $path | bedtools intersect -a cpg_nocontig.bed.gz -b - -sorted -c | cut -f4 | yame pack -f b - $id.cm'
```

Then, run the following commands and tailor the threshold for quality control.

```bash
awk '{print ""$1".cm", $2";"$4;}' controlfiles.tsv | while read fn anno; do yame summary $fn; done > qc.txt

awk '$1!~/QFile/ && $6>5000' qc.txt | awk 'NR==FNR{a[$1]=1;}NR!=FNR&&($1".cm" in a){print $0;}' - controlfiles.tsv | awk '{print ""$1".cm", $2";"$4;}' | sort -k2,2 | while read fn anno; do cat $fn >> merged.cm; yame index -1 $anno merged.cm; done
```

`yame split` splits the samples when provided with `-s` sample name list.

For more help with `split`, run `yame split` in the terminal or check out the
[split help page]({% link docs/subcommands/YAME_split.markdown %}).



