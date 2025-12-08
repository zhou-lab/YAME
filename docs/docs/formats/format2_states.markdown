---
title: Fmt 2 - Categorical States
parent: 1. Storage & Format
nav_order: 2
---

# Format 2 – Categorical State Labels (Chromatin States / Genomic Features)

Format 2 stores **categorical state labels**, one per CpG, using run-length encoding optimized for repeated states across long genomic intervals.

It is ideal for:

- **ChromHMM 15-state / 25-state / full-stack segmentations**
- **Genomic annotations** (promoters, enhancers, gene bodies, exons)
- **Windowed features** (100 kb bins, GC-content partitions)
- **Cell-type or experiment-specific categorical tracks**

Format 2 is the **canonical** YAME format for feature annotation files (`.cm`).

---

## 1. Characteristics of Format 2

Format 2 represents:

- **One categorical label per CpG**
- Labels are stored as **strings**  
  (e.g., `TssA`, `TxFlnk`, `EnhW1`, `Promoter`, `Exon`)
- Internally compressed using run-length encoding:
  - Identical labels across adjacent CpGs stored as a *single block*
- Extremely compact for segment-based features (like ChromHMM)

**Typical extension:** `.cm`  
**Alternate shorthand:** `-fs` (for “state”)

---

## 2. Input Requirements

Input must:

- Have **one state label per row**
- Contain the **same number of rows** as the reference CpG list
- Match the **ordering** of the `.cr` reference file

Example input (`states.txt`):

```text
TssA
TssA
TxFlnk
Tx
Tx
EnhG
EnhW2
EnhW2
Quies
Quies
````

Any string is accepted:

* `Enhancer`, `Promoter`, `Intron`
* `1_TssA`, `2_TssFlnk`, `15_Quies`
* `chr1_100kb_bin42`
* `.` for missing / unassigned regions

(*Note: `.` is treated as a state label, not NA.*)

---

## 3. Packing to Format 2

### 3.1 From a text file

```bash
yame pack -f2 states.txt > chrom_states.cm
```

Shorthand:

```bash
yame pack -fs states.txt > chrom_states.cm
```

---

### 3.2 From ChromHMM BED segmentation (recommended workflow)

Given ChromHMM segments:

```
chr1    10000   10400   1_TssA
chr1    10400   10600   2_TssFlnk
chr1    10600   10800   8_EnhW1
...
```

Convert to Format 2:

```bash
yame unpack cpg_nocontig.cr | gzip > cpg_ref.bed.gz

zcat ChromHMM_segments.bed.gz \
  | bedtools sort \
  | bedtools intersect -a cpg_ref.bed.gz -b - -loj -sorted \
  | bedtools groupby -g 1-3 -c 7 -o first \
  | cut -f4 \
  | yame pack -f2 - > ChromHMM.cm
```

Explanation:

* **`-loj`** ensures all CpGs are retained, even in unannotated regions
  (unannotated sites get `"."`)
* **`groupby first`** resolves overlapping ChromHMM segments

This is the **standard method** used in all YAME examples.

---

### 3.3 From generic genomic features

Any BED file can be converted into categorical states:

```bash
# For example, promoters, enhancers, exons...

bedtools intersect -a cpg_ref.bed.gz -b promoters.bed -loj -sorted \
  | bedtools groupby -g 1-3 -c 7 -o first \
  | cut -f4 \
  | yame pack -f2 - > promoters.cm
```

Labels will be the region names (4th column).

---

### 3.4 From CPRA segments or window partitions

```bash
bedtools makewindows -g hg38.genome -w 100000 \
  | awk '{print $0"\tWin_"NR}' \
  | bedtools intersect -a cpg_ref.bed.gz -b - -loj -sorted \
  | bedtools groupby -g 1-3 -c 7 -o first \
  | cut -f4 \
  | yame pack -f2 - > Win100kb.cm
```

---

## 4. Unpacking Format 2

```bash
yame unpack chrom_states.cm | head
```

Example output:

```text
TssA
TssA
TxFlnk
Tx
Tx
EnhG
EnhW2
EnhW2
Quies
Quies
```

Unpacking simply returns the **string labels** originally packed (lossless).

---

## 5. Integration with Other YAME Commands

Format 2 plays a central role in **feature-level analysis**, especially enrichment testing.

---

### 5.1 `yame summary`

```bash
yame summary -m ChromHMM.cm sample.cg > enrichment.txt
```

For each unique state label, YAME reports:

* `N_mask` — how many CpGs fall in that state
* `N_overlap` — how many of those CpGs meet the query condition
* `Log2OddsRatio` — enrichment of the state for the query
* `Beta` — average methylation level (if query is Format 3)
* `Depth` — depth (if available)

Format 2 masks behave differently from binary masks:

* They **produce one summary row per label**, not a single aggregate.

---

### 5.2 `yame rowsub`

Extract all CpGs belonging to specific states:

```bash
# Extract only Enhancer-like states
echo -e "EnhW1\nEnhW2\nEnhG" > enhancer_states.txt

yame rowsub -R cpg_nocontig.cr -L enhancer_states.txt ChromHMM.cm > enhancer.cm
```

Or using Format 2’s built-in label mode:

```bash
# Add state name as first column in unpacked output
yame rowsub -T -R cpg_nocontig.cr ChromHMM.cm > labeled_output.cm
```

---

### 5.3 `yame mask`

Convert Format 2 to binary mask:

```bash
# Mask all sites that are NOT enhancers
yame mask ChromHMM.cm enhancer_mask.cx -v -o enhancer_only.cm
```

---

### 5.4 `yame dsample`

Downsampling Format 2 keeps **N non-“.” sites** per sample:

```bash
yame dsample -N 50000 ChromHMM.cm > ChromHMM_50k.cm
```

Unannotated sites (`"."`) are typically treated as non-eligible.

---

## 6. When Not To Use Format 2

Avoid Format 2 if:

* You need **numerical values** → use Format 1, 3, or 4
* You need binary 0/1 → use Format 0
* You need structured query/universe sets → use Format 6
* You need NA handling or floating-point → use Format 4

Format 2 is for **categorical labels only**.

---

## 7. Minimal End-to-End ChromHMM Example

```bash
# 1. Prepare reference CpGs
yame unpack cpg_nocontig.cr | gzip > cpg_ref.bed.gz

# 2. Convert ChromHMM segments to Format 2
zcat hg38_ChromHMM_25state.bed.gz \
  | bedtools sort \
  | bedtools intersect -a cpg_ref.bed.gz -b - -loj -sorted \
  | bedtools groupby -g 1-3 -c 7 -o first \
  | cut -f4 \
  | yame pack -f2 - > ChromHMM.cm

# 3. Enrichment of methylation sample
yame summary -m ChromHMM.cm sample.cg > chromhmm_enrichment.txt

# 4. Extract promoter-like states
grep -E "TssA|TssFlnk" ChromHMM.cm > promoters_of_interest.cm
```

---

Format 2 is the **workhorse format** for genomic annotations — especially ChromHMM — and is essential for all feature-based methylation analysis in YAME.

