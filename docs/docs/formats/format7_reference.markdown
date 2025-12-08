---
title: Fmt 7 – Ref Coordinates
parent: 1. Storage & Format
nav_order: 7
---

# Format 7 – CpG Reference Coordinate Files (`.cr`)

Format 7 provides the **reference coordinate system** for all `.cx` files in YAME.

Every `.cg` / `.cm` file—regardless of format (0–6)—must match the row count and ordering of a **Format 7 reference file**.  
This makes Format 7 the backbone of the entire YAME infrastructure.

A `.cr` file contains:

- Chromosome  
- Start position  
- End position  
- CpG ID (e.g., `chr1_10469`)  
- Efficiently delta-encoded and BGZF-compressed for fast access  

---

## 1. Purpose of Format 7

Format 7 answers two fundamental questions:

1. **Where is CpG row *i* located in the genome?**  
2. **What row corresponds to genomic coordinate (chr, pos)?**

All downstream YAME operations assume that **every sample** is aligned to the **same** `.cr` file.

Typical uses:

- Establish genome-wise ordering of CpGs  
- Provide coordinate context to summary operations  
- Support indexing, subsetting, windowing, and masking  
- Serve as the base for feature creation (`format 2`)  
- Ensure consistent row alignment across multi-sample `.cx` files  

---

## 2. What Format 7 Contains

Each row represents **one CpG**:

```

chr1    10468   10469   chr1_10469
chr1    10470   10471   chr1_10471
chr1    10483   10484   chr1_10484
...

````

Columns:

1. **Chromosome**  
2. **Start** (0-based BED convention)  
3. **End** (start+1)  
4. **Name** (or ID), usually `chr_pos1`  

Internally, YAME stores:

- Delta-compressed positions  
- RLE encoding of chromosome boundaries  
- A name dictionary  

This makes `.cr` much smaller than a BED.

---

## 3. Packing to Format 7

### 3.1 From a BED file of CpG sites

```bash
yame pack -f7 cpg_coords.bed > cpg_reference.cr
````

Input BED must contain **at least** the first 3 columns; the 4th is optional but recommended.

Example input (`cpg_coords.bed`):

```text
chr1    10468   10469   chr1_10469
chr1    10470   10471   chr1_10471
chr1    10483   10484   chr1_10484
```

If no name is supplied, YAME generates row IDs automatically.

---

### 3.2 Converting an existing `.cr` file to BED

```bash
yame unpack cpg_reference.cr > cpg_ref.bed
```

This unpacking step is foundational for most workflows:

* Feature construction (Format 2)
* BED-based intersection
* Generating masks or windows

---

## 4. Using Format 7 as the Reference Genome for `.cx` Data

Every `.cg` / `.cm` file must:

* Contain one row per CpG
* Follow the same ordering as the `.cr` file
* Use the same genome build and coordinate convention

To verify alignment:

```bash
yame info sample.cg
```

This prints:

* Row count
* Format ID
* Validity checks
* Compatibility with `.cr` (if provided)

---

### 4.1 Aligning BED-like data to the reference CpGs

This is the standard workflow:

```bash
yame unpack cpg_reference.cr | gzip > cpg_ref.bed.gz

bedtools intersect -a cpg_ref.bed.gz -b input.bed -loj -sorted \
  | cut -f4 \
  | yame pack -fb - > binary.cg
```

This ensures:

* Same ordering as `.cr`
* One CpG per row
* Missing positions explicitly included (value 0)

---

### 4.2 Aligning M/U counts or fraction values

See examples in Format 3 / Format 4 documentation:

```bash
bedtools intersect -a cpg_ref.bed.gz -b mu_counts.bed -loj -sorted \
  | awk '{print $8"\t"$9}' \
  | yame pack -f3 - > sample.cg
```

---

## 5. Integration with YAME Commands

---

### 5.1 `yame rowsub`

Coordinate-based row selection requires `.cr`:

```bash
# Extract CpGs listed by genomic coordinate
yame rowsub -R cpg_reference.cr -L CpG_sites.txt sample.cg > subset.cg
```

### 5.2 `yame summary`

Format 7 itself is not summarized, but masks and queries rely on the order defined by `.cr`.

---

### 5.3 `yame chunk` / `yame chunkchar`

When chunking `.cg` files:

```bash
yame chunk sample.cg chunks/
```

the row boundaries respect `.cr` ordering.

---

### 5.4 Creating features (Format 2) requires `.cr`

Example:

```bash
yame unpack cpg_reference.cr | gzip > cpg_ref.bed.gz

bedtools intersect -a cpg_ref.bed.gz -b chromhmm.bed -loj -sorted \
  | bedtools groupby -g 1-3 -c 7 -o first \
  | cut -f4 \
  | yame pack -f2 - > ChromHMM.cm
```

---

## 6. Choosing or Building a CpG Reference

Most users rely on prebuilt `.cr` files, e.g.:

* **hg19**: CpG_nocontig CR file
* **hg38**: CpG_nocontig CR file
* **mm10**, **mm39**: Available via KYCGKB repositories

Custom genomes:

```bash
# Build your own CpG reference
grep -E -o 'CG' -b genome.fa \
  | awk '{pos=$1; ... construct BED rows ... }' \
  > cpg_coords.bed

yame pack -f7 cpg_coords.bed > cpg_reference.cr
```

---

## 7. When NOT to Use Format 7

Do NOT use Format 7 for:

* Storing methylation data → use Formats 0,1,3,4,5,6
* Storing mask or feature annotations → use Format 2
* Storing per-sample or bulk CpG calls → use Formats 0/3/4

Format 7 is strictly a **reference coordinate container**.

---

## 8. Minimal End-to-End Example

```bash
# 1. Unpack the reference CpGs
yame unpack cpg_nocontig.cr | gzip > cpg_ref.bed.gz

# 2. Align your BED data to reference CpGs
bedtools intersect -a cpg_ref.bed.gz -b peaks.bed -loj -sorted \
  | cut -f4 \
  | yame pack -fb - > peaks.cg

# 3. Create a feature file using ChromHMM
zcat ChromHMM.bed.gz \
  | bedtools intersect -a cpg_ref.bed.gz -b - -loj -sorted \
  | bedtools groupby -g 1-3 -c 7 -o first \
  | cut -f4 \
  | yame pack -f2 - > chromhmm.cm

# 4. Summarize enrichment
yame summary -m chromhmm.cm peaks.cg > enrich.txt
```

---

Format 7 is the **foundation** of all YAME workflows:
It defines CpG identity, ordering, genomic position, and compatibility across all `.cg` and `.cm` files.
A correct and consistent `.cr` file ensures that your entire methylation analysis pipeline remains coherent, efficient, and reproducible.

