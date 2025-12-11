---
title: 1. Storage & Format
nav_order: 1
---

# 1. Storage & Format
{: .no_toc }

Convert between text and compressed CX binary formats for efficient storage and analysis.
{: .fs-5 .fw-300 }

---

## Overview

YAME’s `pack` and `unpack` commands provide bidirectional conversion between human-readable text formats (BED/TSV/etc.) and the compressed CX binary formats.

### CX formats:

- Dramatically reduce storage requirements (often 10–100×)
- Are optimized for different methylation / feature data types
- Are interoperable across all YAME commands

### Most workflows follow this pattern:

1. **Align** your data to a CpG reference coordinate file (`.cr`, format 7)  
2. **Pack** into an appropriate CX format (0/2/3/4/5/6)  
3. Use other YAME commands (`summary`, `rowsub`, `rowop`, `subset`, etc.)  
4. Optionally **unpack** to text for external tools

---

## Format Overview

YAME currently supports the following CX format family:

| Format                                                             | Code      | Typical Ext | Best for                                                                     |
|--------------------------------------------------------------------|-----------|-------------|------------------------------------------------------------------------------|
| [**Format 0**]({% link docs/formats/format0_binary.markdown %})    | `0` / `b` | `.cg`       | Binary presence/absence (DMR sites, ChIP-seq peaks, generic 0/1 tracks)      |
| [**Format 1**]({% link docs/formats/format1_integer.markdown %})   | `1`       | `.cg`       | Integer values with RLE (count tracks, QC metrics, per-CpG integer signals)  |
| [**Format 2**]({% link docs/formats/format2_states.markdown %})    | `2` / `s` | `.cm`       | Chromatin states, genomic annotations, gene features, windows/bins           |
| [**Format 3**]({% link docs/formats/format3_mu.markdown %})        | `3` / `m` | `.cg`       | M/U read counts from bisulfite sequencing                                    |
| [**Format 4**]({% link docs/formats/format4_beta.markdown %})      | `4`       | `.cg`       | Continuous methylation values (beta/fraction), array/WGBS imputed values     |
| [**Format 6**]({% link docs/formats/format6_su.markdown %})        | `6`       | `.cx`       | Set + universe representation for enrichment; sparse single-cell methylation |
| [**Format 7**]({% link docs/formats/format7_reference.markdown %}) | `7`       | `.cr`       | CpG genomic reference coordinates (required for all `.cg` / `.cm` files)     |

All CX formats use BGZF compression and share a consistent internal structure (`cdata_t` blocks), enabling uniform handling across different data types.

## Pack / Unpack Basics

The core commands are:

```bash
# Pack text → CX
yame pack -f<format> [options] <input.txt> [output.cx]
# Read from stdin:
cat data.txt | yame pack -fb - data.cg
# Unpack CX → text
yame unpack [options] <input.cx>
````
  
### ⚠️ **Your input MUST match reference CpG coordinates exactly:**

1. *Same number of rows* as reference CpGs
2. *Same order* as reference CpGs
3. *One value per CpG* in reference

**How to ensure alignment:**

```bash
# Step 1: Get reference coordinates
yame unpack cpg_nocontig.cr | gzip > cpg_ref.bed.gz
# Step 2: Intersect your BED file with reference
bedtools intersect -a cpg_ref.bed.gz -b your_data.bed -sorted -c | \
  cut -f4 | \
  yame pack -fb - > aligned_output.cg
# the following two output should match in dimension
yame info cpg_nocontig.cr
yame info aligned_output.cg
```

A more complete workflow (including reference coordinate alignment) is described in each format’s page.


---

## Practical Workflows

### Workflow 1: Process Bisulfite Sequencing Data

Complete pipeline from bisulfite sequencing output to analysis:

```bash
# Assuming you have M and U counts from your pipeline
# Format: chr start end M U

# 1. Extract M and U, align with reference
bedtools intersect -a cpg_ref.bed.gz -b methylation_calls.bed -loj -sorted | \
  awk '{if ($8==".") print "0\t0"; else print $8"\t"$9}' | \
  yame pack -f3 - > sample.cg

# 2. Verify data quality
yame info sample.cg
yame summary sample.cg

# 3. Perform enrichment analysis
yame summary -m ChromHMM.cm sample.cg > chromatin_enrichment.txt
yame summary -m genes.cm sample.cg > gene_enrichment.txt

# 4. If needed, unpack for external tools
yame unpack -f 5 sample.cg > sample_cov5.txt
```

---

### Workflow 2: Create Multi-Feature Database

Build comprehensive feature database:

```bash
#!/bin/bash
# Create comprehensive feature database for hg38

# 1. Prepare reference
yame unpack cpg_nocontig.cr | gzip > cpg_ref.bed.gz

# 2. Download and process multiple features
declare -A features=(
  ["ChromHMM_15"]="ChromHMM_15state.bed.gz"
  ["CpG_Islands"]="cpgIslandExt.bed.gz"
  ["Promoters"]="promoters_2kb.bed.gz"
  ["Enhancers"]="enhancers_merged.bed.gz"
  ["TFBS"]="tfbs_combined.bed.gz"
)

# 3. Create individual feature files
for name in "${!features[@]}"; do
  file="${features[$name]}"
  echo "Processing $name..."
  
  zcat "$file" | bedtools sort | \
    bedtools intersect -a cpg_ref.bed.gz -b - -loj -sorted | \
    bedtools groupby -g 1-3 -c 7 -o first | \
    cut -f4 | \
    yame pack -f2 - > "${name}.cm"
done

echo "Feature files created:"
ls -lh *.cm
```

---

### Workflow 3: Convert Array Data to CX Format

Convert Illumina array data:

```bash
#!/bin/bash
# Convert 450k/EPIC array data to CX format

# Assuming you have: array_data.txt with columns: ProbeID, Beta

# 1. Get array manifest with probe coordinates
# manifest.bed format: chr start end ProbeID

# 2. Intersect with reference CpGs
bedtools intersect -a cpg_ref.bed.gz -b manifest.bed -loj -sorted > probe_cpg_map.bed

# 3. Map betas to CpG positions
join -1 4 -2 1 -t$'\t' \
  <(sort -k4,4 probe_cpg_map.bed) \
  <(sort -k1,1 array_data.txt) | \
  sort -k2,2 -k3,3n | \
  cut -f7 | \
  awk '{if ($1=="") print "NA"; else print $1}' | \
  yame pack -f4 - > array_sample.cg

echo "Array data converted to CX format"
yame info array_sample.cg
```

---

## Best Practices

1. **Always align with reference first**
   ```bash
   bedtools intersect -a cpg_ref.bed.gz -b your_data.bed -sorted -c
   ```

2. **Verify dimensions match**
   ```bash
   # Count reference CpGs
   zcat cpg_ref.bed.gz | wc -l
   
   # Count your data rows
   wc -l your_data.txt
   ```

3. **Handle missing values properly**
   - Format 3: Use M=0, U=0 for no coverage
   - Format 4: Use "NA" for missing values
   - Format 2: Use "." for unassigned CpGs

