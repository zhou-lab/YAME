---
title: Fmt 3 – M & U Counts
parent: 1. Storage & Format
nav_order: 3
---

# Format 3 – Methylated & Unmethylated Counts

Format 3 stores **paired M/U (methylated/unmethylated) read counts per CpG** from bisulfite-based assays (e.g., WGBS, RRBS).

- **Typical extension:** `.cg`
- **Input:** two integers per CpG: `M` (methylated) and `U` (unmethylated)
- **Best for:** high-coverage sequencing where you want to preserve read counts and coverage

---

## 1. Input Requirements

Before packing into format 3, your data must be:

- Aligned to a **CpG reference coordinate file** (format 7, e.g. `cpg_nocontig.cr`)
- One row per CpG, in the **same order** and of the **same length** as the reference
- Each row: two non-negative integers: `M` and `U`

Example input (`mu_counts.txt`):

```text
12    3
8     8
0     15
20    2
0     0
5     10
````

Interpretation per line:

* Column 1: methylated count (M)
* Column 2: unmethylated count (U)
* Both must be ≥ 0; `M = U = 0` means “no coverage / missing”

---

## 2. Packing to Format 3

### 2.1 From an aligned M/U table

If you already have a table of `M U` pairs, aligned to the CpG reference:

```bash
yame pack -f3 mu_counts.txt > sample.cg
```

This creates a compressed `.cg` file in **format 3**.

---

### 2.2 From BED-like methylation calls

If your pipeline outputs per-CpG counts as BED:

```bash
# mu_calls.bed columns (example):
# chr  start  end  M  U

bedtools intersect -a cpg_ref.bed.gz -b mu_calls.bed -loj -sorted \
  | awk '{ if ($8 == ".") print "0\t0"; else print $8"\t"$9 }' \
  | yame pack -f3 - > sample.cg
```

Key points:

* `-loj` (left outer join) ensures **every reference CpG** appears exactly once.
* CpGs not covered in `mu_calls.bed` are assigned `0 0` (no coverage).
* The `awk` produces two columns (`M`, `U`) suitable for `yame pack -f3`.

---

## 3. Unpacking Format 3

By default, `yame unpack` converts M/U counts into **beta** and **coverage**:

```bash
yame unpack sample.cg | head
```

Example output:

```text
0.800   15
0.500   16
0.000   15
0.909   22
NA      0
0.333   15
```

Columns:

1. **Beta** – `M / (M + U)`; `NA` if `M + U == 0`
2. **Coverage** – `M + U`

You can filter by minimum coverage with `-f`:

```bash
# Only output CpGs with coverage ≥ 5
yame unpack -f 5 sample.cg > sample_cov5.txt
```

This is convenient for QC and downstream tools that expect beta + coverage.

---

## 4. Integration with Other YAME Commands

Format 3 works with most downstream operations:

### 4.1 `yame summary`

```bash
yame summary sample.cg
yame summary -m feature.cm sample.cg
```

`summary` will compute per-sample and per-feature:

* Number of CpGs in universe / query / mask
* Overlap counts and log2 odds ratio
* **Average beta** in feature regions
* **Average depth** in feature regions

---

### 4.2 `yame rowop` (row-wise operations)

Some common operations on format 3:

```bash
# Per-CpG mean beta across samples
yame rowop -o mean sample.cg > mean_beta.tsv

# Per-CpG standard deviation of beta across samples
yame rowop -o std sample.cg > beta_std.tsv

# Sum M and U across samples (true pseudobulk counts)
yame rowop -o musum sample.cg bulk.cg

# Binarize and sum (vote-based pseudobulk)
yame rowop -o binasum -c 3 sample.cg bulk_binasum.cg
```

Quick interpretations:

* `-o mean`
  For each CpG, computes the **mean beta** across samples that have coverage ≥ `-c` (default 1).
  Output: `beta_mean    n_samples_used`.

* `-o std`
  For each CpG, computes **standard deviation** of beta across samples with coverage ≥ `-c`.
  Output: `beta_sd    n_samples_used`.

* `-o musum`
  Sums M and U **directly** across samples, preserving counts.
  Output is a new format 3 `.cg` file.

* `-o binasum`
  Converts each sample to a binary call by comparing M vs U (and coverage ≥ `-c`), then counts methylated vs unmethylated “votes” across samples.

---

### 4.3 `yame dsample` (downsampling)

You can create downsampled versions of format 3 data to simulate reduced coverage or increased sparsity:

```bash
# Keep 50,000 covered CpGs per sample, fixed seed
yame dsample -N 50000 -s 1 sample.cg > sample_N50k.cg
```

Behavior for format 3:

* Eligible sites: CpGs with `M + U > 0`
* Randomly selects up to `N` eligible sites; the rest are **masked** by setting `M = U = 0`.

---

### 4.4 `yame rowsub` (row subsetting)

Subset CpGs (rows) from a format 3 file:

```bash
# By index list
yame rowsub -l row_ids.txt sample.cg > subset.cg

# By coordinate list and row coordinate file
yame rowsub -R cpg_nocontig.cr -L CpG_sites.txt sample.cg > subset.cg

# By binary mask
yame rowsub -m mask.cx sample.cg > masked_subset.cg
```

The output remains format 3.

---

## 5. When to Use Format 3 (vs Other Formats)

Choose **Format 3** when:

* You have **bisulfite sequencing data** with per-CpG counts.
* You care about **coverage** and want to model uncertainty (beta at low coverage vs high coverage).
* You want the option to derive beta, binary calls, or more complex statistics later.

Consider other formats when:

* You only need beta values: use **Format 4** for more compact storage.
* You want only 0/1 calls: use **Format 0**.
* You’re working with sparse single-cell binary data with an explicit universe: consider **Format 6**.

---

## 6. Minimal End-to-End Example

```bash
# 1. Prepare reference CpG coordinates (once per genome/build)
yame unpack cpg_nocontig.cr | gzip > cpg_ref.bed.gz

# 2. Align M/U calls to reference CpGs
bedtools intersect -a cpg_ref.bed.gz -b mu_calls.bed -loj -sorted \
  | awk '{ if ($8 == ".") print "0\t0"; else print $8"\t"$9 }' \
  | yame pack -f3 - > sample.cg

# 3. QC and summary
yame info sample.cg
yame summary sample.cg

# 4. Export beta and coverage for external tools (coverage ≥ 10)
yame unpack -f 10 sample.cg > sample_cov10.txt
```

This is the typical pattern for ingesting WGBS/RRBS data into YAME using Format 3.

