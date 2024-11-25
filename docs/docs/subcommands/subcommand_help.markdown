---
title: YAME Subcommands
nav_order: 9
has_children: true
permalink: docs/subcommands
---

# YAME Subcommands

Usage for all YAME subcommands (printed with `YAME`).

```bash

Program: yame (Yet Another Methylation tool) - whole genome DNA methylation data management.
Version: 0.3.20230904
Contact: Wanding Zhou<wanding.zhou@pennmedicine.upenn.edu>

Usage:   yame <command> [options]

Available commands:
     pack         - Pack data into a cx file.
     unpack       - Unpack data from a cx file.
     subset       - Subset samples from a cx file.
     rowsub       - Subset rows a cx file using an index list file.
     info         - Display basic parameter of the cx file.
     summary      - calculate summary, with or without masks.
     index        - Index samples in a cx file.
     split        - Split multi-sample data into single-sample data.
     chunk        - Chunk data into smaller fragments.
     chunkchar    - Chunk text data into smaller fragments.
     rowop        - Perform operations on rows, e.g., sum binary values.
     mask         - Mask methylation data by setting masked record to M=U=0.
     dsample      - Downsample methylation data by setting unsampled records to M=U=0.

```
