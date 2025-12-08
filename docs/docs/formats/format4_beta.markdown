---
title: Fmt 4 – Fraction Values
parent: 1. Storage & Format
nav_order: 4
---

# Format 4 – Continuous Methylation Fractions (Beta Values)

Format 4 stores **continuous numeric values**, typically representing:

- **DNA methylation beta values** (0.0–1.0)
- **Methylation fractions** computed from M/U counts
- **Imputed methylation levels**
- **Normalized array intensities converted to β-values**
- **Any numeric score per CpG where NA is allowed**

Format 4 is the **right choice** when you need floating-point precision and optional NA handling.

---

## 1. Characteristics of Format 4

Format 4 provides:

- One **float** per CpG
- Optional **NA** values
- Compression using “NA + run-length encoding” (NA-RLE)
- Very efficient storage for array data and imputed datasets

**Typical extension:** `.cg`  
**Format flag:** `-f4`

---

## 2. Input Requirements

Input must have:

- **One row per CpG**
- A value in one of the following forms:
  - Float between 0 and 1 (inclusive)
  - `"NA"` meaning missing
  - Scientific notation allowed (e.g., `3e-2`)

Example input (`beta_values.txt`):

```text
0.75
0.33
0.88
NA
0.50
0.92
````

Valid input also includes:

```text
0
1
0.12345
0.9999
```

---

## 3. Packing to Format 4

### 3.1 Pack from a simple text vector

```bash
yame pack -f4 beta_values.txt > beta.cg
```

---

### 3.2 From array data (e.g., Illumina 450k, EPIC)

Suppose you have a `ProbeID → Beta` table:

```bash
# array_data.txt example:
# cg00001234    0.45
# cg00002456    0.71
# cg00003522    NA
```

Convert to Format 4:

```bash
# Intersect array probes with CpG reference
join -1 4 -2 1 -t$'\t' \
  <(sort -k4,4 cpg_ref_with_ids.bed) \
  <(sort -k1,1 array_data.txt) \
  | sort -k2,2 -k3,3n \
  | cut -f5 \
  | yame pack -f4 - > array_sample.cg
```

---

### 3.3 From M/U counts (Format 3 → Format 4)

To convert a Format 3 `.cg` into beta values:

```bash
yame unpack -a sample.cg > beta_cov.txt   # outputs beta and coverage

cut -f1 beta_cov.txt \
  | yame pack -f4 - > sample.beta.cg
```

Or directly:

```bash
yame unpack sample.cg | cut -f1 | yame pack -f4 - > sample.beta.cg
```

---

## 4. Unpacking Format 4

```bash
yame unpack beta.cg | head
```

Example output:

```text
0.75
0.33
0.88
NA
0.50
0.92
```

Format 4 unpacks to a **single float value per line**, identical to input (lossless).

---

## 5. Integration with Other YAME Commands

---

### 5.1 `yame summary`

Format 4 supports enrichment and window summarization:

```bash
yame summary -m features.cm beta.cg
```

Outputs include:

* `Beta` → **mean of numeric values inside the mask**
* `Depth` → always “NA” for Format 4
* `N_query` → number of **non-NA** CpGs
* `N_overlap` → number of CpGs inside mask with **non-NA** values

---

### 5.2 `yame rowop`

Useful numeric operations:

```bash
# Per-CpG mean across samples
yame rowop -o mean beta_multi.cg > mean.tsv

# Per-CpG standard deviation across samples
yame rowop -o std beta_multi.cg > std.tsv
```

Format 4 behaves identically to Format 3 for `mean`/`std`, except that NA is allowed.

---

### 5.3 `yame dsample`

Downsampling Format 4 is based on **non-NA** values:

```bash
yame dsample -N 50000 beta.cg > beta_50k.cg
```

* Eligible: CpGs where value is not NA
* Non-selected: set to **NA**

---

### 5.4 `yame rowsub`

```bash
# Subset by mask
yame rowsub -m promoter.cm beta.cg > beta.promoters.cg

# Subset by coordinate list
yame rowsub -R cpg_nocontig.cr -L CpG_sites.txt beta.cg > subset.cg
```

---

### 5.5 `yame mask`

Masking Format 4 replaces masked positions with **NA**:

```bash
yame mask beta.cg lowqual.cx -o beta.filtered.cg
```

---

## 6. When NOT to Use Format 4

Use a different format if:

* You need M/U counts → **Format 3**
* Your values are binary → **Format 0**
* Your labels are categorical → **Format 2**
* You need a structured universe/query → **Format 6**
* You need NA-free numeric integers → **Format 1**

Format 4 is the **only** format allowing floating-point values.

---

## 7. Minimal End-to-End Example

```bash
# 1. Prepare CPM reference
yame unpack cpg_nocontig.cr | gzip > cpg_ref.bed.gz

# 2. Convert array data to Format 4
join -1 4 -2 1 -t$'\t' \
  <(sort -k4,4 cpg_ref_with_ids.bed) \
  <(sort -k1,1 array_data.txt) \
  | sort -k2,2 -k3,3n \
  | cut -f5 \
  | yame pack -f4 - > array_sample.cg

# 3. Summarize across genomic features
yame summary -m ChromHMM.cm array_sample.cg > array_enrichment.txt

# 4. Extract CpGs for gene promoters
yame rowsub -m promoters.cm array_sample.cg > promoter_beta.cg

# 5. Replace low-confidence CpGs with NA using mask
yame mask array_sample.cg lowconf.cx -o array_filtered.cg
```

---

Format 4 is the **default for all continuous methylation values**, including array β-values and imputed sequencing fractions, with full support for NA, summarization, subsetting, and efficient compression.

