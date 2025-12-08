---
title: Enrichment Testing
nav_order: 2
---

# Enrichment Testing with YAME
{: .no_toc }

Test methylation enrichment across genomic features to identify biological associations.
{: .fs-6 .fw-300 }

## Table of Contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Overview

Enrichment testing identifies whether methylation patterns are statistically associated with specific genomic features (e.g., chromatin states, histone modifications, transcription factor binding sites, regulatory elements). YAME's `summary` command performs rapid enrichment testing by comparing observed overlaps with expected overlaps based on background frequencies.

### What You'll Learn

- Prepare differential methylation results for enrichment testing
- Convert BED files to CX format
- Perform enrichment testing with and without custom backgrounds
- Interpret enrichment results and statistical significance
- Use built-in differential calling functionality

---

## Prerequisites

Before starting, ensure you have:

1. **Differential methylation results** - A BED file containing significant CpG coordinates
2. **bedtools** - Install from [bedtools.readthedocs.io](https://bedtools.readthedocs.io/en/latest/content/installation.html)
3. **Reference CpG coordinates** - Download `.cr` files from KYCG:
   - [hg38 reference](https://github.com/zhou-lab/KYCGKB_hg38)
   - [mm10 reference](https://github.com/zhou-lab/KYCGKB_mm10)
4. **Feature annotations** - Download `.cm` feature files from KYCG (chromatin states, histone marks, etc.)

### Unpack Reference Coordinates

First, convert the reference `.cr` file to a BED format:

```bash
yame unpack cpg_nocontig.cr | gzip -c > cpg_nocontig.bed.gz
```

---

## Basic Workflow

### Step 1: Prepare Your Query Set

Convert your BED file of differential methylation sites to CX format. This involves:
1. Sorting your BED file
2. Intersecting with reference CpG coordinates
3. Packing into binary format

```bash
bedtools sort -i yourfile.bed | \
  bedtools intersect -a cpg_nocontig.bed.gz -b - -sorted -c | \
  cut -f4 | \
  yame pack -fb - > yourfile.cg
```

**Explanation:**
- `bedtools sort` - Ensures BED file is sorted for efficient intersection
- `bedtools intersect -c` - Counts overlaps with each reference CpG
- `cut -f4` - Extracts the CpG identifier column
- `yame pack -fb` - Packs binary data (format b) into CX format

**Note:** If your BED file is already sorted, you can skip the `bedtools sort` step.

### Step 2: Run Enrichment Testing

Test enrichment against genomic features using the `-m` option:

```bash
yame summary -m feature.cm yourfile.cg > enrichment_results.txt
```

**Feature files** (`.cm` format) are available from KYCG and include:
- Chromatin states (ChromHMM, Roadmap Epigenomics)
- Histone modifications (ChIP-seq peaks)
- Transcription factor binding sites
- Gene annotations (promoters, exons, introns, UTRs)
- Regulatory elements (enhancers, insulators)
- CpG islands and shores
- Repetitive elements

You can also create custom feature files using [`yame pack`]({% link docs/pack_unpack.markdown %}).

---

## Understanding Output

The enrichment results contain the following columns:

| Column | Description |
|--------|-------------|
| **QFile** | Query file name |
| **Query** | Sample name in query file |
| **MFile** | Feature/mask file name |
| **Mask** | Feature name being tested |
| **N_univ** | Total CpGs in universe (reference) |
| **N_query** | CpGs covered in query |
| **N_mask** | CpGs in feature |
| **N_overlap** | CpGs in both query and feature |
| **Log2OddsRatio** | Enrichment strength (log₂ scale) |
| **Beta** | Average methylation level in feature |
| **Depth** | Approximate coverage depth in feature |

### Interpreting Log2OddsRatio

The log₂ odds ratio quantifies enrichment strength:

- **≥ 2.0** - Strong enrichment (4-fold or greater)
- **1.0 - 2.0** - Moderate enrichment (2-4 fold)
- **0.5 - 1.0** - Weak enrichment (1.4-2 fold)
- **≈ 0** - No enrichment (random overlap)
- **< 0** - Depletion (under-representation)

Higher positive values indicate stronger association between your query set and the feature being tested.

### Statistical Significance

For formal significance testing, use the Fisher's exact test implementation from the [sesame R package](https://www.bioconductor.org/packages/release/bioc/html/sesame.html):

```r
# In R using sesame package
library(sesame)

# Map YAME output columns to Fisher's test parameters
ND <- N_mask      # Feature size
NQ <- N_query     # Query size  
NDQ <- N_overlap  # Overlap count
NU <- N_univ      # Universe size

# Perform Fisher's exact test
testEnrichmentFisherN(ND, NQ, NDQ, NU)
```

An alternative R implementation is also available in the YAME GitHub repository under the `R/` folder.

---

## Built-in Differential Calling

YAME includes a built-in pairwise differential methylation caller:

```bash
yame pairwise -H 1 -c 10 \
  <(yame subset sample1.cg sample1) \
  <(yame subset sample2.cg sample2) \
  -o differential_output.cg
```

**Parameters:**
- `-H 1` - Controls directionality (1 = hypermethylation in first sample)
- `-c 10` - Minimum coverage threshold (default: 10)
- `-o` - Output file name

**Directionality options for `-H`:**
- `1` - Hypermethylated in first sample
- `-1` - Hypomethylated in first sample  
- `0` - Either direction

The output `.cg` file can be directly used for enrichment testing as described above.

---

## Enrichment with Custom Background

### Why Use a Custom Background?

The choice of background is crucial for accurate enrichment testing:

- **Default behavior**: Compares against all CpGs in the reference
- **Custom background**: Compares against CpGs measured in your specific experiment

Using a custom background (e.g., all CpGs covered by your assay) provides more accurate enrichment estimates by accounting for technical biases like:
- Array probe coverage
- Sequencing accessibility
- Coverage variation across genomic regions

### Method 1: Using Mask (Recommended)

The simplest approach uses the `mask` command to restrict analysis to your universe of interest:

```bash
yame mask -c query.cg universe.cg | \
  yame summary -m feature.cm - > enrichment_results.txt
```

**How it works:**
- `query.cg` - Your differential methylation sites
- `universe.cg` - All CpGs measured in your experiment
- The mask operation restricts enrichment testing to only those CpGs in your universe

### Method 2: Two-Column Format (Legacy)

For more complex scenarios, you can create a two-column CX file:

**Step 1: Intersect both query and background with reference**

```bash
# Query set intersection
bedtools intersect -a cpg_nocontig.bed.gz -b query.bed -sorted -c | \
  cut -f4 > query_intersect.bed

# Background set intersection  
bedtools intersect -a cpg_nocontig.bed.gz -b background.bed -sorted -c | \
  cut -f4 > background_intersect.bed
```

**Step 2: Combine into two-column format**

```bash
paste query_intersect.bed background_intersect.bed | \
  yame pack -f6 - > query_background.cg
```

**Step 3: Run enrichment testing**

```bash
yame summary -m feature.cm query_background.cg > enrichment_results.txt
```

In this format:
- **Column 1**: Query set (sites of interest)
- **Column 2**: Background/universe (all measured sites)

---

## Creating Custom Feature Files

You can create your own feature annotations using `yame pack`. See the detailed example in the [pack_unpack documentation]({% link docs/pack_unpack.markdown %}#example-generating-feature-files).

Quick example for a custom BED file:

```bash
bedtools intersect -a cpg_nocontig.bed.gz -b custom_features.bed -loj -sorted | \
  bedtools groupby -g 1-3 -c 7 -o first | \
  cut -f4 | \
  yame pack -f2 - > custom_features.cm
```

---

## Best Practices

1. **Always use sorted BED files** - Ensures efficient bedtools operations
2. **Match reference coordinates** - Use the same genome build (hg38/mm10) throughout
3. **Choose appropriate background** - Use experiment-specific backgrounds when possible
4. **Consider multiple testing correction** - Apply FDR/Bonferroni when testing many features
5. **Validate enrichments** - Cross-reference with biological literature and databases
6. **Check coverage requirements** - Ensure sufficient coverage for meaningful statistics

---

## Troubleshooting

**Problem**: "No overlaps found"
- Ensure genome builds match (hg38 vs hg19, mm10 vs mm9)
- Verify BED files are properly formatted (0-based coordinates)
- Check that CpG identifiers match the reference format

**Problem**: "All features show similar enrichment"
- May indicate improper background selection
- Consider using experiment-specific background

**Problem**: "Very high log2 odds ratios (>10)"
- May indicate low coverage or small sample sizes
- Review N_overlap and N_mask values for adequacy

---

## Additional Resources

- [KYCG hg38 features](https://github.com/zhou-lab/KYCGKB_hg38)
- [KYCG mm10 features](https://github.com/zhou-lab/KYCGKB_mm10)
- [bedtools documentation](https://bedtools.readthedocs.io/)
- [sesame R package](https://bioconductor.org/packages/sesame/)

---

