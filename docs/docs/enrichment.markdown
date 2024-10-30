---
title: Methylation Data Enrichment
nav_order: 1
---

# Tutorial: Enrichment Testing with YAME

## Prerequisites

1. A bed file containing the output significant coordinates from differential analysis
2. Installed [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html) on your system
3. A reference coordinate bed file (We have provided hg38 and mm10 CpG reference coordinate annotations .cr on KYCG github) [mm10](https://github.com/zhou-lab/KYCGKB_mm10)/[hg38](https://github.com/zhou-lab/KYCGKB_hg38). 

```bash
yame unpack cpg_nocontig.cr | gzip > cpg_nocontig.bed.gz
```

## File preparation

First, we will pack the bedfile into a .cx format. If the input bedfile is already sorted, you can start with the intersect step. Check out the [bedtools instersect help page](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html) if you encounter any problems at this step.

```bash
bedtools sort -i yourfile.bed | bedtools intersect -a cpg_nocontig.bed.gz -b - -sorted -c | cut -f4 | yame pack -fb - > yourfile.cg   
```
## Enrichment testing

Then we simply run [`yame summary`]({% link docs/summarize.markdown %}) with `-m` feature file for enrichment testing. We have provided comprehensive enrichment feature files, and you can download them from th KYCG github page [mm10](https://github.com/zhou-lab/KYCGKB_mm10)/[hg38](https://github.com/zhou-lab/KYCGKB_hg38). You can also create your own feature file with [`yame pack`]({% link docs/pack_unpack.markdown %}).

```bash
yame summary -m feature.cm yourfile.cg > yourfile.txt
```

## Output interpretations

Detailed information of the output columns can be found on the [`yame summary`]({% link docs/summarize.markdown %}) page. Basically, a higher log2oddsratio indicates a stronger association between the feature being tested and the query set. Generally, a large log2 odds ratio is typically considered to be around 2 or greater, with values between 1 and 2 often being viewed as potentially important and worthy of further investigation, while values around 0.5 might be considered a small effect size. For significance testing, [seasame](https://www.bioconductor.org/packages/release/bioc/html/sesame.html) R package provided the testEnrichmentFisherN function, which is also provided in the yame github R page. The four input parameters correspond to the four columns from yame summary output.
```
ND = N_mask
NQ = N_query
NDQ = N_overlap
NU = N_universe
```

## built in differential calling function

```bash
yame pairwise -H 1 -c 10 <(yame subset sample1.cg sample1) <(yame subset sample2.cg sample2) -o output.cg
```
-H controls directionality and -c controls minimum coverage. 

## Enrichment testing with background (new)

Selecting the appropriate background for enrichment testing is crucial because it can significantly impact the interpretation of the results. Usually, we use the background set that is measured in the experiment under different conditions. 

```bash
yame mask -c query.cg universe.cg | yame summary -m feature.cm - > yourfile.txt
```

## Enrichment testing with background (old)

To enable enrichment testing with background, we need to prepare the .cx file to include the two bed files (one is for the query set, one is for the background set). 

```bash
bedtools intersect -a cpg_nocontig.bed.gz -b query.bed -sorted -c | cut -f4 > query_intersect.bed
bedtools intersect -a cpg_nocontig.bed.gz -b background.bed -sorted -c | cut -f4 > background_intersect.bed
```
And then we put them together where the first column is the query and the second column is the background set. 

```bash
paste query_intersect.bed background_intersect.bed | yame pack -f6 - > query_background.cg
```

Conducting enrichment testing with background using `yame summary` is the same procedure as shown above. 