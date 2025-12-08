---
title: Data Storage & Format Conversion
nav_order: 2
---

# Methylation Data Packing and Unpacking
{: .no_toc }

Convert between text and compressed CX binary formats for efficient storage and analysis.
{: .fs-6 .fw-300 }

## Table of Contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Overview

YAME's `pack` and `unpack` commands provide bidirectional conversion between human-readable text formats and the compressed CX binary format family. The CX formats dramatically reduce storage requirements (often 10-100x compression) while enabling rapid analysis operations.

### Why Use CX Formats?

- **Massive compression** - Reduce storage from gigabytes to megabytes
- **Fast I/O** - BGZF compression enables rapid reading/writing
- **Uniform structure** - All formats share consistent internal architecture
- **Optimized for methylation** - Each format designed for specific data types
- **Seamless integration** - Works with all YAME commands and external tools

### Key Concepts

- **Packing** - Convert text data (BED, TSV, CSV) to compressed binary CX format
- **Unpacking** - Convert CX format back to human-readable text
- **Format specification** - Different CX formats optimized for different data types
- **Reference coordinates** - All data must align with reference CpG coordinate order
- **BGZF compression** - Block-compressed format allowing random access

---

## Understanding CX Formats

### The CX Format Family

YAME supports 8 different CX formats, each optimized for specific data types:

| Format | Type Code | File Ext | Data Type | Compression Method | Best For |
|--------|-----------|----------|-----------|-------------------|----------|
| **Format 0** | `b` or `0` | `.cg` | Binary | 1 byte per 8 CpGs | Presence/absence, DMR calls |
| **Format 1** | `1` | `.cg` | Value + RLE | Value (1B) + RLE (2B) | Compressed integer counts |
| **Format 2** | `s` or `2` | `.cm` | State + RLE | State labels + Index RLE | Chromatin states, genomic features |
| **Format 3** | `3` | `.cg` | MU RLE | M/U counts + Ladder byte | Bisulfite sequencing (M,U) |
| **Format 4** | `4` | `.cg` | Fraction + NA-RLE | Float fraction + NA encoding | Methylation beta values |
| **Format 5** | `5` | `.cg` | 2-bit + NA-RLE | Ternary values {0,1,2} | Differential calls (hypo/unchanged/hyper) |
| **Format 6** | `6` | `.cg` | 2-bit boolean | Set and universe flags | Query + background for enrichment |
| **Format 7** | `7` | `.cr` | Compressed BED | Delta-encoded coordinates | CpG reference coordinates |

### Choosing the Right Format

**Quick Decision Tree:**

```
What type of data do you have?

├─ Binary (yes/no, present/absent)
│  └─> Format 0 (b)
│
├─ Chromatin states or categorical features
│  └─> Format 2 (s)
│
├─ M and U count pairs from bisulfite sequencing
│  └─> Format 3
│
├─ Beta values or methylation fractions (0.0-1.0)
│  └─> Format 4
│
├─ Differential calls (hypomethylated/unchanged/hypermethylated)
│  └─> Format 5
│
├─ Query + background for enrichment testing
│  └─> Format 6
│
└─ Reference genomic coordinates
   └─> Format 7
```

---

## yame pack

### Basic Syntax

```bash
yame pack -f<format> [options] <input_file> [output_file]
```

### Essential Options

| Option | Description | Example |
|--------|-------------|---------|
| `-f<format>` | Format specification (required) | `-fb`, `-f2`, `-f3` |
| `-` | Read from stdin | `cat data.txt \| yame pack -fb -` |
| `output.cg` | Output filename (optional) | `yame pack -fb input.txt output.cg` |

### Critical Input Requirement

⚠️ **Your input MUST match reference CpG coordinates exactly:**

1. **Same number of rows** as reference CpGs
2. **Same order** as reference CpGs
3. **One value per CpG** in reference

**How to ensure alignment:**

```bash
# Step 1: Get reference coordinates
yame unpack cpg_nocontig.cr | gzip > cpg_ref.bed.gz

# Step 2: Intersect your BED file with reference
bedtools intersect -a cpg_ref.bed.gz -b your_data.bed -sorted -c | \
  cut -f4 | \
  yame pack -fb - > aligned_output.cg
```

---

## Format-Specific Packing Examples

### Format 0/b - Binary Data

**Use case:** Presence/absence, binary methylation calls, DMR indicators

**Input:** One line per CpG, values 0 or 1

```bash
# Example input (binary_data.txt):
# 1
# 0
# 1
# 1
# 0

# Pack from text file
yame pack -fb binary_data.txt > binary_output.cg

# Pack from BED file via bedtools
bedtools sort -i dmr_sites.bed | \
  bedtools intersect -a cpg_ref.bed.gz -b - -sorted -c | \
  cut -f4 | \
  yame pack -fb - > dmr_binary.cg
```

**Compression:** ~32x reduction (8 CpGs per byte)

---

### Format 2/s - Chromatin States

**Use case:** Categorical annotations, chromatin states, genomic features

**Input:** One line per CpG, containing state label (text string)

```bash
# Example input (states.txt):
# TssA
# TssFlnk
# TxFlnk
# Tx
# Tx
# EnhG
# Enh

# Pack state data
yame pack -f2 states.txt > chromatin_states.cm

# Or use 's' shorthand
yame pack -fs states.txt > chromatin_states.cm
```

**Special property:** Excellent for feature files with many categories

---

### Format 3 - M/U Sequencing Counts

**Use case:** Bisulfite sequencing methylated (M) and unmethylated (U) read counts

**Input:** Two tab-separated columns per line (M count, U count)

```bash
# Example input (mu_counts.txt):
# 12    3
# 8     8
# 0     15
# 20    2
# 0     0
# 5     10

# Pack M/U count data
yame pack -f3 mu_counts.txt > sequencing_data.cg
```

**Compression:** Highly efficient for typical sequencing coverage distributions

**Note:** M and U can be any non-negative integers

---

### Format 4 - Methylation Fractions

**Use case:** Beta values, methylation fractions, continuous values 0.0-1.0

**Input:** One line per CpG, float values between 0 and 1 (or NA)

```bash
# Example input (beta_values.txt):
# 0.75
# 0.33
# 0.88
# NA
# 0.50
# 0.92

# Pack fraction data
yame pack -f4 beta_values.txt > methylation_fractions.cg
```

**Special handling:** 
- NA values are properly encoded
- Preserves precision for methylation fractions
- Efficient for array data

---

### Format 5 - Differential Methylation

**Use case:** Three-state differential calls

**Input:** One line per CpG, values {0, 1, 2} only
- **0** = Hypomethylated
- **1** = Unchanged / Not significant
- **2** = Hypermethylated

```bash
# Example input (diff_calls.txt):
# 0
# 1
# 2
# 1
# 0
# 2

# Pack differential calls
yame pack -f5 diff_calls.txt > differential.cg
```

**Compression:** 2 bits per CpG (4 CpGs per byte)

---

### Format 6 - Set and Universe

**Use case:** Enrichment testing with query set and background/universe

**Input:** Two tab-separated columns
- Column 1: Query set (0 or 1)
- Column 2: Universe/background (0 or 1)

```bash
# Example input (query_universe.txt):
# 1    1   (in query, in universe)
# 0    1   (not in query, in universe)
# 1    1   (in query, in universe)
# 0    0   (not in query, not in universe)

# Pack query + universe
yame pack -f6 query_universe.txt > query_bg.cg
```

**Use case:** More precise enrichment testing (see [enrichment guide]({% link docs/enrichment.markdown %}))

---

### Format 7 - Reference Coordinates

**Use case:** CpG reference coordinate files

**Input:** BED format with CpG identifiers

```bash
# Example input (cpg_coords.bed):
# chr1    10468   10469   chr1_10469
# chr1    10470   10471   chr1_10471
# chr1    10483   10484   chr1_10484

# Pack reference coordinates
yame pack -f7 cpg_coords.bed > cpg_reference.cr
```

**Output:** `.cr` file (CpG reference)

**Compression:** Delta encoding of positions for maximum compression

---

## Creating Feature Files - Complete Guide

Feature files (format 2) are essential for enrichment testing and feature-based analysis. Here's the complete workflow.

### Example 1: ChromHMM Chromatin States

This detailed example shows how to create a ChromHMM feature file from full-stack chromatin state annotations.

#### Input Data Format

Your ChromHMM BED file (`hg38_chromhmm_100_segments.bed.gz`) should look like:

```
chr1    10000   10400   1_TssA
chr1    10400   10600   2_TssFlnk
chr1    10600   10800   8_EnhW1
chr1    10800   12800   15_Quies
chr1    12800   13000   8_EnhW1
chr1    13000   13200   7_EnhG2
chr1    13200   14800   15_Quies
```

**Column format:**
1. Chromosome
2. Start position
3. End position  
4. Chromatin state label

#### Step-by-Step Workflow

**Step 1: Prepare reference CpG coordinates**

```bash
# Unpack reference .cr file to BED format
yame unpack cpg_nocontig.cr | gzip > cpg_nocontig.bed.gz

# Verify it worked
zcat cpg_nocontig.bed.gz | head -3
# chr1    10468   10469   chr1_10469
# chr1    10470   10471   chr1_10471
# chr1    10483   10484   chr1_10484
```

**Step 2: Intersect CpGs with chromatin states**

```bash
zcat hg38_chromhmm_100_segments.bed.gz | \
  bedtools sort | \
  bedtools intersect -a cpg_nocontig.bed.gz -b - -loj -sorted | \
  bedtools groupby -g 1-3 -c 7 -o first | \
  cut -f4 | \
  yame pack -f2 - > ChromHMM_FullStack.cm
```

**Command breakdown:**
- `zcat ... | bedtools sort` - Ensure chromatin state BED is sorted
- `bedtools intersect -loj` - Left outer join: keep ALL reference CpGs
  - CpGs overlapping states get the state label
  - CpGs not overlapping any state get "." (NA)
- `bedtools groupby -g 1-3 -c 7 -o first` - For CpGs in multiple states, take first
- `cut -f4` - Extract the state label (4th column is state name)
- `yame pack -f2 -` - Pack as format 2 (categorical state data)

**Step 3: Verify the output**

```bash
# Check file info
yame info ChromHMM_FullStack.cm

# Preview state distribution
yame unpack ChromHMM_FullStack.cm | sort | uniq -c | sort -rn | head
```

**Step 4: Use for enrichment testing**

```bash
yame summary -m ChromHMM_FullStack.cm query.cg > enrichment_results.txt
```

---

### Example 2: Histone Modification Peaks

Create feature file from ChIP-seq peaks:

```bash
# Single histone mark
bedtools intersect -a cpg_ref.bed.gz -b H3K4me3_peaks.bed -sorted -c | \
  cut -f4 | \
  yame pack -fb - > H3K4me3.cm

# Multiple marks - see sample.markdown for merging multiple features
```

---

### Example 3: Gene Annotations

Create features for different gene regions:

```bash
# Promoters (2kb upstream of TSS)
bedtools intersect -a cpg_ref.bed.gz -b promoters.bed -sorted -c | \
  cut -f4 | \
  yame pack -fb - > promoters.cm

# Exons
bedtools intersect -a cpg_ref.bed.gz -b exons.bed -sorted -c | \
  cut -f4 | \
  yame pack -fb - > exons.cm

# Gene bodies
bedtools intersect -a cpg_ref.bed.gz -b gene_bodies.bed -sorted -c | \
  cut -f4 | \
  yame pack -fb - > gene_bodies.cm
```

---

### Example 4: Genomic Bins/Windows

Create fixed-size genomic windows:

```bash
# 1. Create 100kb windows
bedtools makewindows -g hg38.genome -w 100000 | \
  awk '{print $0"\tchr"$1":"$2/100000}' > windows_100kb.bed

# 2. Assign CpGs to windows
yame unpack cpg_ref.cr | \
  bedtools intersect -a - -b windows_100kb.bed -loj -sorted | \
  bedtools groupby -g 1-3 -c 7 -o first | \
  cut -f4 | \
  yame pack -f2 - > bins_100kb.cm

# 3. Use for single-cell aggregation
yame summary -m bins_100kb.cm single_cell.cg > binned_methylation.txt
```

---

### Example 5: Custom Region Sets

Create feature from any custom BED file:

```bash
# Your custom regions (e.g., enhancers, TADs, etc.)
# custom_regions.bed format:
# chr1    100000  102000  enhancer_1
# chr1    150000  153000  enhancer_2

bedtools intersect -a cpg_ref.bed.gz -b custom_regions.bed -loj -sorted | \
  bedtools groupby -g 1-3 -c 7 -o first | \
  cut -f4 | \
  yame pack -f2 - > custom_features.cm
```

---

## yame unpack

Convert CX binary format back to human-readable text.

### Basic Syntax

```bash
yame unpack [options] <cx_file>
```

### Common Options

| Option | Description | Example |
|--------|-------------|---------|
| `-a` | Unpack all samples | `yame unpack -a multi.cg` |
| `-f <N>` | Filter by min coverage | `yame unpack -f 5 data.cg` |
| `-s <name>` | Unpack specific sample | `yame unpack -s sample1 data.cg` |
| `-i <index>` | Unpack by sample index | `yame unpack -i 0 data.cg` |

---

## Unpacking Examples

### Format 3 - M/U Count Data

```bash
# Unpack all samples with coverage >= 1
yame unpack -a -f 1 sequencing.cg

# Output format (tab-separated):
# beta    coverage
# 0.75    4
# 0.33    3
# 0.88    8
# NA      0
# 0.50    2
```

**Columns:**
- **Beta**: Methylation fraction (M / (M+U))
- **Coverage**: Total reads (M + U)

### Format 0/b - Binary Data

```bash
# Unpack binary data
yame unpack binary_data.cg

# Output (one value per line):
# 1
# 0
# 1
# 1
# 0
```

### Format 2/s - State Data

```bash
# Unpack chromatin states
yame unpack ChromHMM.cm | head -10

# Output (state labels):
# 1_TssA
# 1_TssA
# 3_TxFlnk
# 4_Tx
# 4_Tx
# 7_Enh
# 15_Quies
```

### Format 7 - Reference Coordinates

```bash
# Unpack reference coordinates to BED
yame unpack cpg_nocontig.cr > cpg_coordinates.bed

# Output (BED format):
# chr1    10468   10469   chr1_10469
# chr1    10470   10471   chr1_10471
# chr1    10483   10484   chr1_10484
```

### Unpacking Specific Samples

```bash
# Check available samples first
yame info multi_sample.cg

# Unpack specific sample
yame unpack -s "sample_A" multi_sample.cg > sample_A.txt

# Unpack first sample (index 0)
yame unpack -i 0 multi_sample.cg > first_sample.txt
```

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

### Workflow 4: Quality Control Pipeline

Systematic QC for methylation data:

```bash
#!/bin/bash
# QC pipeline for CX files

CX_FILE=$1

echo "=== Quality Control Report for $CX_FILE ==="
echo

# 1. Basic file info
echo "File Information:"
yame info $CX_FILE
echo

# 2. Coverage statistics
echo "Coverage Statistics:"
yame summary $CX_FILE
echo

# 3. CpG island enrichment (should be enriched)
echo "CpG Island Enrichment:"
yame summary -m CpG_islands.cm $CX_FILE | awk 'NR>1 {print "Log2OR:", $9}'
echo

# 4. Check coverage distribution
echo "Coverage Distribution:"
yame unpack -a $CX_FILE | \
  awk '{cov=$2} cov>0 {if(cov<5) low++; else if(cov<10) med++; else high++} 
       END{print "Low(<5):", low, "Med(5-10):", med, "High(>10):", high}'
echo

# 5. Methylation distribution
echo "Methylation Distribution:"
yame unpack -a $CX_FILE | \
  awk '$1!="NA" {beta=$1; if(beta<0.3) hypo++; else if(beta<0.7) inter++; else hyper++}
       END{print "Hypo(<0.3):", hypo, "Inter(0.3-0.7):", inter, "Hyper(>0.7):", hyper}'

echo
echo "=== QC Complete ==="
```

---

## Best Practices

### For Packing

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

3. **Choose appropriate format**
   - Binary (0/1) → Format 0
   - States/categories → Format 2
   - M,U counts → Format 3
   - Fractions → Format 4

4. **Handle missing values properly**
   - Format 3: Use M=0, U=0 for no coverage
   - Format 4: Use "NA" for missing values
   - Format 2: Use "." for unassigned CpGs

5. **Test with small subset first**
   ```bash
   head -1000 data.txt | yame pack -fb - > test.cg
   yame info test.cg
   yame unpack test.cg | head
   ```

### For Feature Files

1. **Use -loj (left outer join) for completeness**
   ```bash
   bedtools intersect -a cpg_ref.bed.gz -b features.bed -loj -sorted
   ```

2. **Handle overlapping features explicitly**
   ```bash
   bedtools groupby -g 1-3 -c 7 -o first  # Take first
   # or
   bedtools groupby -g 1-3 -c 7 -o last   # Take last
   # or
   bedtools groupby -g 1-3 -c 7 -o collapse  # Keep all (comma-separated)
   ```

3. **Quality control your features**
   ```bash
   yame summary feature.cm
   # Check: N_mask (should be reasonable number of CpGs)
   ```

4. **Document feature creation**
   ```bash
   # Add comment in script
   # Feature: ChromHMM 15-state model
   # Source: ENCODE, hg38
   # Date created: 2025-01-15
   # Command: [paste command here]
   ```

### For Unpacking

1. **Filter by coverage when appropriate**
   ```bash
   yame unpack -f 5 data.cg  # Only CpGs with coverage ≥ 5
   ```

2. **Pipe to reduce storage**
   ```bash
   yame unpack data.cg | process_data.py  # Don't save intermediate
   ```

3. **Verify unpacked data**
   ```bash
   yame unpack data.cg | head  # Check format
   yame unpack data.cg | wc -l  # Check row count
   ```

---

## Format Selection Guide

Choose the optimal format for your data:

| Your Data Type | Format | Command Example |
|----------------|--------|-----------------|
| DMR sites (yes/no) | 0/b | `yame pack -fb` |
| Differential calls (hypo/no/hyper) | 5 | `yame pack -f5` |
| Chromatin states | 2/s | `yame pack -f2` |
| ChIP-seq peaks | 0/b | `yame pack -fb` |
| M and U counts | 3 | `yame pack -f3` |
| 450k array betas | 4 | `yame pack -f4` |
| WGBS betas | 4 | `yame pack -f4` |
| Gene annotations | 2/s | `yame pack -f2` |
| Genomic bins | 2/s | `yame pack -f2` |
| Query + background | 6 | `yame pack -f6` |
| CpG coordinates | 7 | `yame pack -f7` |

---

## Troubleshooting

### Common Issues

**Problem**: "Dimension mismatch" error

**Solution:**
```bash
# Check dimensions
zcat cpg_ref.bed.gz | wc -l  # Reference
wc -l your_data.txt            # Your data

# These must match! If not:
bedtools intersect -a cpg_ref.bed.gz -b your_data.bed -loj -sorted | \
  bedtools groupby -g 1-3 -c 7 -o first | \
  cut -f4 > aligned_data.txt
```

**Problem**: Packed file is unexpectedly large

**Solutions:**
- Check if you're using the right format (binary data shouldn't use format 3)
- Verify data is properly aligned
- Consider if data is actually compressible

**Problem**: Unpacked data shows all NAs

**Solutions:**
```bash
# Check file format
yame info yourfile.cg

# Verify data exists
yame summary yourfile.cg

# Try unpacking without filters
yame unpack yourfile.cg | head
```

**Problem**: Feature labels not showing in enrichment

**Solutions:**
```bash
# Verify format 2 was used
yame info feature.cm

# Check labels
yame unpack feature.cm | sort | uniq -c

# Recreate with proper format
yame pack -f2 states.txt > feature.cm
```

**Problem**: bedtools intersect produces wrong output

**Solutions:**
```bash
# Ensure proper sorting
bedtools sort -i file1.bed > file1_sorted.bed
bedtools sort -i file2.bed > file2_sorted.bed

# Use -sorted flag
bedtools intersect -a file1_sorted.bed -b file2_sorted.bed -sorted

# Check chromosome naming (chr1 vs 1)
head file1.bed
head file2.bed
```

**Problem**: Memory errors when packing large files

**Solutions:**
```bash
# Process in chunks
split -l 5000000 large_data.txt chunk_
for chunk in chunk_*; do
  yame pack -fb $chunk > ${chunk}.cg
done

# Combine (if needed)
cat chunk_*.cg > combined.cg
```

---

## Performance Optimization

### Speed Tips

1. **Use pipes to avoid I/O**
   ```bash
   # Slow (intermediate files)
   bedtools intersect ... > temp.bed
   cut -f4 temp.bed > temp2.txt
   yame pack -fb temp2.txt > output.cg
   
   # Fast (piped)
   bedtools intersect ... | cut -f4 | yame pack -fb - > output.cg
   ```

2. **Parallel processing for multiple files**
   ```bash
   # Process multiple samples in parallel
   parallel -j 8 'yame pack -f3 {} > {.}.cg' ::: *.txt
   ```

3. **Pre-sort large files**
   ```bash
   # Sort once, use many times
   bedtools sort -i large_file.bed > large_file_sorted.bed
   ```

### Storage Tips

1. **Don't double-compress**
   ```bash
   # CX files are already compressed, don't gzip them
   yame pack -fb data.txt > output.cg  # Good
   yame pack -fb data.txt | gzip > output.cg.gz  # Redundant
   ```

2. **Choose appropriate format**
   - Format 0 for binary: ~32x compression
   - Format 3 for M/U: ~10-20x compression
   - Format 2 for states: Varies by number of categories

3. **Remove intermediate text files**
   ```bash
   yame pack -fb data.txt > data.cg && rm data.txt
   ```

---

## Validation and Testing

Always validate your packed files:

```bash
#!/bin/bash
# Validation script

CX_FILE=$1

echo "Validating $CX_FILE..."

# 1. Check file is valid CX format
if yame info $CX_FILE > /dev/null 2>&1; then
  echo

```
