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

YAME (Yet Another Methylation Encoder) is designed for efficient sequence-level DNA methylation data management, capable of handling both bulk and single-cell DNA methylome workflows. It introduces a family of compact binary formats (**CX formats**) that represent methylation values, MU counts, categorical states, fraction data, masks, and genomic coordinates in a uniform compressed structure. Use Cases include:
    
- *Single-cell methylome analysis* - Efficiently store and analyze sparse single-cell data
- *Bulk methylome processing* - Fast operations on large cohorts
- *Enrichment testing* - Test methylation enrichment across genomic features
- *Feature aggregation* - Summarize methylation over bins, chromatin states, or custom regions
- *Pseudobulk generation* - Merge single cells into cluster-level pseudobulks
- *Differential methylation* - Identify and test differentially methylated sites

---

## Tutorials & Workflows

- **[Storage & Format]({% link docs/storage.markdown %})** - Working with CX formats
- **[Summarize & Encode]({% link docs/summarize.markdown %})** - Calculate statistics and aggregations
- **[Test Enrichment]({% link docs/enrichment.markdown %})** - Test methylation enrichment across genomic features
- **[Subset Rows]({% link docs/sample.markdown %})** - Extract samples and regions
- **[Aggregate Row-wise]({% link docs/rowop.markdown %})** - Merge pseudobulks and perform calculations
- **[Combine, Split & Index]({% link docs/sample.markdown %})** - Handle multi-sample datasets
- **[Mask Data]({% link docs/dsample.markdown %})** - Test methods at different sparsity levels

## Installation

```bash
# Option 1: Install via Conda (Recommended)
conda install yame -c bioconda
# Option 2: Build from Source
git clone https://github.com/zhou-lab/YAME.git
cd YAME
make
```

## Quick Start

```bash
# Pack binary methylation data
yame pack -fb yourfile.bed > yourfile.cg
# Pack MU count data (M and U columns)
yame pack -f3 methylation_counts.txt > yourfile.cg
# Summarize methylation data
yame summary yourfile.cg
# Test enrichment over genomic features
yame summary -m ChromHMM_states.cm yourfile.cg > enrichment_results.txt
# Subset samples from multi-sample data
yame subset -l sample_list.txt yourfile.cg > subset.cg
# Merge single cells into pseudobulks
yame subset -l cluster1_cells.txt single_cell.cg | yame rowop - -o binasum > pseudobulk.cg
```

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

- **Command Help**: Run `yame` or `yame <command>` for usage information
- **Issues**: Report bugs on [GitHub Issues](https://github.com/zhou-lab/YAME/issues)
- **KYCG Resources**: Download reference coordinates and feature files from [KYCG hg38](https://github.com/zhou-lab/KYCGKB_hg38) or [KYCG mm10](https://github.com/zhou-lab/KYCGKB_mm10)

