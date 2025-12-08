---
# Feel free to add content and custom Front Matter to this file.
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults

layout: default
title: Home
description: "YAME (Yet Another Methylation Encoder) is a fast and lightweight toolkit for storing, manipulating, and analyzing large-scale DNA methylation data at the sequence level. It introduces the CX format family for ultra-compact representation of methylation values, states, and genomic coordinates."
nav_order: 0
---

# YAME - Yet Another Methylation Encoder
{: .fs-9 }

A fast and lightweight toolkit for sequence-level DNA methylation analysis
{: .fs-6 .fw-300 }

[Get Started](#quick-start){: .btn .btn-primary .fs-5 .mb-4 .mb-md-0 .mr-2 }
[View on GitHub](https://github.com/zhou-lab/YAME){: .btn .fs-5 .mb-4 .mb-md-0 }

---

## Overview

YAME (Yet Another Methylation Encoder) is designed for efficient sequence-level DNA methylation data management, capable of handling both bulk and single-cell DNA methylome workflows. It introduces a family of compact binary formats (**CX formats**) that represent methylation values, MU counts, categorical states, fraction data, masks, and genomic coordinates in a uniform compressed structure.

### Key Features

- **🚀 Ultra-fast performance** - Efficient C implementation with BGZF compression
- **💾 High compression ratios** - Compact binary CX formats dramatically reduce storage requirements
- **📈 Massive scalability** - Handle hundreds of thousands of single cells
- **🔧 Comprehensive toolkit** - Pack, unpack, subset, downsample, merge, and analyze methylation data
- **🧬 Versatile data support** - MU counts, binary methylation, chromatin states, fractions, differential calls, and CpG coordinates
- **🔗 Seamless integration** - Works with bedtools, KYCG knowledge base, and standard methylation workflows

### Use Cases

- **Single-cell methylome analysis** - Efficiently store and analyze sparse single-cell data
- **Bulk methylome processing** - Fast operations on large cohorts
- **Enrichment testing** - Test methylation enrichment across genomic features
- **Feature aggregation** - Summarize methylation over bins, chromatin states, or custom regions
- **Pseudobulk generation** - Merge single cells into cluster-level pseudobulks
- **Differential methylation** - Identify and test differentially methylated sites

---

## Quick Start

### Installation

**Option 1: Install via Conda (Recommended)**
```bash
conda install yame -c bioconda
```

**Option 2: Build from Source**
```bash
git clone https://github.com/zhou-lab/YAME.git
cd YAME
make
```

To run YAME from anywhere, move it to your PATH:
```bash
mv ./yame ~/bin/yame
```

### Basic Usage Example

**1. Pack methylation data into CX format**
```bash
# Pack binary methylation data
yame pack -fb yourfile.bed > yourfile.cg

# Pack MU count data (M and U columns)
yame pack -f3 methylation_counts.txt > yourfile.cg
```

**2. Summarize methylation data**
```bash
yame summary yourfile.cg
```

**3. Test enrichment over genomic features**
```bash
yame summary -m ChromHMM_states.cm yourfile.cg > enrichment_results.txt
```

**4. Subset samples from multi-sample data**
```bash
yame subset -l sample_list.txt yourfile.cg > subset.cg
```

**5. Merge single cells into pseudobulks**
```bash
yame subset -l cluster1_cells.txt single_cell.cg | yame rowop - -o binasum > pseudobulk.cg
```

---

## CX Format Family

YAME uses the **CX format family** to efficiently store different types of methylation-related data:

| Format | Extension | Data Type | Use Case |
|--------|-----------|-----------|----------|
| Format 0 | `.cg` | Binary (1 byte/8 CpGs) | Presence/absence of methylation |
| Format 1 | `.cg` | Value + RLE | Compressed methylation values |
| Format 2 | `.cm` | State text + Index RLE | Chromatin states, genomic features |
| Format 3 | `.cg` | MU RLE + Ladder byte | M and U sequencing counts |
| Format 4 | `.cg` | Fraction + NA-RLE | Methylation fractions |
| Format 5 | `.cg` | 2-bits + NA-RLE | Differential methylation (0,1,2) |
| Format 6 | `.cg` | 2-bits boolean | Set and universe indicators |
| Format 7 | `.cr` | Compressed BED | CpG reference coordinates |

All CX formats use BGZF compression and share a consistent internal structure (`cdata_t` blocks), enabling uniform handling across different data types.

---

## Command Reference

The following commands provide the core functionality of YAME. Click on each command for detailed documentation.

### Data Storage & Conversion
- [`pack`]({% link docs/pack_unpack.markdown %}) - Pack data into CX format
- [`unpack`]({% link docs/pack_unpack.markdown %}) - Unpack CX files to text format

### Data Subsetting
- [`subset`]({% link docs/subset.markdown %}) - Extract specific samples from multi-sample files
- [`rowsub`]({% link docs/subset.markdown %}) - Extract specific genomic regions or CpG sites

### Data Summarization & Analysis
- [`info`]({% link docs/summarize.markdown %}) - Display basic parameters of CX files
- [`summary`]({% link docs/summarize.markdown %}) - Calculate statistics and test enrichment over features

### Sample Operations
- [`index`]({% link docs/sample.markdown %}) - Index samples in multi-sample files
- [`split`]({% link docs/sample.markdown %}) - Split multi-sample data into individual files
- [`chunk`]({% link docs/sample.markdown %}) - Divide data into smaller fragments
- [`chunkchar`]({% link docs/sample.markdown %}) - Chunk text data into fragments

### Row Operations
- [`rowop`]({% link docs/rowop.markdown %}) - Perform row-wise operations (mean, sum, std, binarization, co-methylation)

### Data Manipulation
- [`mask`]({% link docs/dsample.markdown %}) - Mask methylation data by setting masked records to M=U=0
- [`dsample`]({% link docs/dsample.markdown %}) - Downsample to a fixed number of sites

### Differential Analysis
- [`pairwise`]({% link docs/enrichment.markdown %}) - Built-in differential methylation calling

---

## Tutorials & Workflows

- **[Enrichment Testing]({% link docs/enrichment.markdown %})** - Test methylation enrichment across genomic features
- **[Data Storage & Format Conversion]({% link docs/pack_unpack.markdown %})** - Working with CX formats
- **[Subsetting Data]({% link docs/subset.markdown %})** - Extract samples and regions
- **[Summarization & Analysis]({% link docs/summarize.markdown %})** - Calculate statistics and aggregations
- **[Row Operations]({% link docs/rowop.markdown %})** - Merge pseudobulks and perform calculations
- **[Sample Management]({% link docs/sample.markdown %})** - Handle multi-sample datasets
- **[Downsampling]({% link docs/dsample.markdown %})** - Test methods at different sparsity levels

---

## Reference

**Goldberg\*, Fu\*, Atkins, Moyer, Lee, Deng, Zhou† (2025)**  
KnowYourCG: Facilitating Base-level Sparse Methylome Interpretation.  
*Science Advances*  
DOI: [10.1126/sciadv.adw3027](https://www.science.org/doi/10.1126/sciadv.adw3027)

---

## Acknowledgements

This work is supported by **NIH/NIGMS 5R35GM146978**.

YAME integrates with the [KYCG knowledge base](https://github.com/zhou-lab/KYCG) for comprehensive methylation feature analysis.

---

## Getting Help

- **Documentation**: Browse the menu for detailed guides
- **Command Help**: Run `yame` or `yame <command>` for usage information
- **Issues**: Report bugs on [GitHub Issues](https://github.com/zhou-lab/YAME/issues)
- **KYCG Resources**: Download reference coordinates and feature files from [KYCG hg38](https://github.com/zhou-lab/KYCGKB_hg38) or [KYCG mm10](https://github.com/zhou-lab/KYCGKB_mm10)

---
