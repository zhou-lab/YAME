---
# Feel free to add content and custom Front Matter to this file.
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults

layout: default
title: Home
description: "YAME (Yet Another MEthylation tool) is a methylation toolset designed for sequence-level DNA methylation data management. It is a command-line C program capable of performing sequence-level enrichment testing, row operations (such as merging pseudobulks), downsampling, and other related tasks with ultra fast speed."
nav_order: 0
---

# YAME - Understand Sequencing Level Methylation Data
{: .fs-9 }
[View it on GitHub](https://github.com/zhou-lab/YAME){: .btn .fs-5 .mb-4 .mb-md-0 }

---

YAME (Yet Another MEthylation tool) is a methylation toolset designed for sequence-level DNA methylation data management. It is a command-line C program capable of performing sequence-level enrichment testing, row operations (such as merging pseudobulks), downsampling, and other related tasks with ultra fast speed.

## Reference

David Goldberg, Hongxiang Fu, *et al.*,
KnowYourCG: supervised learning of sparse DNA methylomes,
*Manuscript in submission*

## Installation Guide
The source can be retrieved by:

```bash
# git
git clone https://github.com/zhou-lab/YAME.git
cd YAME
```
After retrieving the source code, building YAME simply proceeds as follows:

```bash
make
```

If you want to run yame from anywhere without specifying the path, move it to a directory that is already in your PATH, such as ~/bin

```bash
mv ./yame ~/bin/yame
```

## Overview of Functionalities

The following list provides an overview of the various functionalities provided by
`YAME`. You can also find much of this by typing `YAME` in the terminal after installation.

  - [`pack`]({% link docs/pack_unpack.markdown %}) Pack data into a cx file. 
  - [`unpack`]({% link docs/pack_unpack.markdown %}) Unpack data from a cx file.
  - [`subset`]({% link docs/subset.markdown %}) Subset samples from a cx file.
  - [`rowsub`]({% link docs/subset.markdown %}) Subset rows a cx file using an index list file.
  - [`info`]({% link docs/summarize.markdown %})   Display basic parameter of the cx file.
  - [`summary`]({% link docs/summarize.markdown %}) calculate summary, with or without masks.
  - [`index`]({% link docs/sample.markdown %}) Index samples in a cx file.
  - [`split`]({% link docs/sample.markdown %}) Split multi-sample data into single-sample data.
  - [`chunk`]({% link docs/sample.markdown %}) Chunk data into smaller fragments.
  - [`chunkchar`]({% link docs/sample.markdown %}) Chunk text data into smaller fragments.
  - [`rowop`]({% link docs/rowop.markdown %}) Perform operations on rows, e.g., sum binary values.
  - [`mask`]({% link docs/dsample.markdown %}) Mask methylation data by setting masked record to M=U=0.
  - [`dsample`]({% link docs/dsample.markdown %}) Downsample methylation data by setting unsampled records to M=U=0.

## Acknowledgement
  - This work is supported by NIH/NIGMS 5R35GM146978.

---