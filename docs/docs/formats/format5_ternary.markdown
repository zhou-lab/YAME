---
title: Fmt 5 – Differential Calls
parent: 1. Storage & Format
nav_order: 5
---

# Format 5 – Ternary Differential Methylation Calls

Format 5 stores a **single ternary value per CpG**, typically representing:

- **0** → hypomethylated  
- **1** → unchanged / not significant  
- **2** → hypermethylated  

This format is compact (2 bits per CpG), extremely fast to read, and ideal for storing *categorical differential methylation results*.

---

## 1. When to Use Format 5

Use Format 5 when your data:

- Represents **three discrete states**  
  (e.g., {hypo, unchanged, hyper})
- Comes from:
  - differential methylation analysis (DMRs)
  - EWAS with discretized outcomes
  - statistical calling pipelines
- Does **not** require:
  - continuous values (use Format 4)
  - counts (use Format 3)
  - categorical labels of arbitrary strings (use Format 2)

Format 5 is ideal for **large differential methylation compendia** because it’s tiny and fast.

---

## 2. Input Requirements

Input must include:

- **One line per CpG**
- A value from the set **{0, 1, 2}**

Example input (`diff_calls.txt`):

```text
0
1
2
1
0
2
````

Interpretation (default):

| Value | Meaning         |
| ----- | --------------- |
| 0     | Hypomethylated  |
| 1     | No change       |
| 2     | Hypermethylated |

(*These semantic labels are user-defined; YAME treats them as integers.*)

---

## 3. Packing to Format 5

### 3.1 Pack from text

```bash
yame pack -f5 diff_calls.txt > differential.cg
```

### 3.2 Generate Format 5 from an external DMR pipeline

If your DMR results are in BED or TSV, convert them into a CpG-level ternary track by intersecting with reference CpGs:

```bash
# Example: columns = chr, start, end, diff_call (0/1/2)

bedtools intersect -a cpg_ref.bed.gz -b dmr_calls.bed -loj -sorted \
  | awk '{ if ($8 == ".") print "1"; else print $8 }' \
  | yame pack -f5 - > diff.cg
```

Interpretation:

* Sites lacking a call become `1` (unchanged)
  (You can change this behavior depending on your workflow.)

---

## 4. Unpacking Format 5

```bash
yame unpack differential.cg | head
```

Output:

```text
0
1
2
1
0
2
```

Format 5 unpacks to integers **exactly as originally stored**.

---

## 5. Integration with YAME Commands

Format 5 integrates with major YAME tools.

---

### 5.1 `yame summary`

Use Format 5 to compute enrichment of hypo / hyper states across features:

```bash
yame summary -m ChromHMM.cm differential.cg > diff_summary.txt
```

Behavior:

* `N_query` = number of CpGs with **non-1** values (or all CpGs, depending on interpretation)
* `N_overlap` = number of CpGs in the mask with relevant values
* `Beta` = **mean ternary state** (0–2) across the mask
  (This is typically interpretable as “direction of change”)
* `Depth` = NA (unused)

If you want to treat only **hypo (0) vs hyper (2)** as binary, you may preprocess the data.

---

### 5.2 `yame rowop`

Format 5 supports row-wise aggregation across multiple samples.

#### Summing ternary states:

```bash
yame rowop -o binasum differential_multi.cg > diff_votes.cg
```

Interpretation:

* Converts:

  * 0 → unmethylated vote
  * 2 → methylated vote
  * 1 → ambiguous (ignored)
* Produces a Format 3 pseudobulk count file.

#### Mean or standard deviation:

```bash
yame rowop -o mean differential_multi.cg > mean_state.tsv
yame rowop -o std differential_multi.cg > state_variability.tsv
```

Output:

* **mean**: average state (0–2) across samples
* **std**: variability in differential calling

---

### 5.3 `yame dsample`

Downsampling Format 5 retains **N non-1 CpGs** per sample:

```bash
yame dsample -N 10000 diff.cg > diff_10k.cg
```

Rules:

* Eligible = CpGs where state ≠ 1 (i.e., changes)
* Non-selected = set to **1** (unchanged)

---

### 5.4 `yame rowsub`

```bash
# Extract CpGs that are differential
yame rowsub -m diff_mask.cx diff.cg > diff_only.cg

# Extract genomic subregions
yame rowsub -R cpg_nocontig.cr -L CpG_sites.txt diff.cg > subset.cg
```

---

### 5.5 `yame mask`

To mask out (remove) certain CpGs:

```bash
yame mask diff.cg lowqual_mask.cx -o filtered_diff.cg
```

Masked values become **1** (interpreted as “not differential”).
Or use `-v` to invert mask logic.

---

## 6. When Not to Use Format 5

Choose a different format if:

* You need floating-point values → **Format 4**
* You need per-CpG M/U counts → **Format 3**
* You need arbitrary string labels → **Format 2**
* You need a sparse binary representation with an explicit universe → **Format 6**

Format 5 is strictly for **three-state encoded data**.

---

## 7. Minimal End-to-End Example

```bash
# 1. Prepare reference CpGs
yame unpack cpg_nocontig.cr | gzip > cpg_ref.bed.gz

# 2. Convert DMR calls into ternary CpG-level states
bedtools intersect -a cpg_ref.bed.gz -b dmr_calls.bed -loj -sorted \
  | awk '{ if ($8 == ".") print "1"; else print $8 }' \
  | yame pack -f5 - > diff.cg

# 3. Evaluate enrichment of differential states across ChromHMM annotations
yame summary -m ChromHMM.cm diff.cg > diff_enrichment.txt

# 4. Extract hypermethylated CpGs (state = 2)
grep "^2$" diff.cg > hyper_sites.txt

# 5. Downsample for benchmarking
yame dsample -N 20000 diff.cg > diff_20k.cg
```

---

Format 5 is the compact, efficient format for **discrete differential methylation calls**, offering excellent compression and fast downstream feature-based analysis.

