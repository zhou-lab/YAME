---
title: Fmt 0 – Binary Data
parent: 1. Storage & Format
nav_order: 0
---

# Format 0 – Binary Presence/Absence Data

Format 0 stores **binary data**, one bit per CpG, representing presence/absence or a simple methylation state.

- **Typical extension:** `.cg`
- **Input:** a vector of **0/1 values**, one per CpG
- **Best for:**  
  - Binary methylation calls  
  - Presence/absence features  
  - DMR indicators  
  - Peak overlap masks  
  - Cell-level binary CpG accessibility  

Format 0 is the **smallest** and **fastest** CX format: 8 CpGs per byte (~32× compression over text).

---

## 1. Input Requirements

Your input must be:

- Aligned to a **CpG reference file** (`.cr`, format 7)
- Exactly **one line per CpG**
- Each line containing **0** or **1**

Example input (`binary_data.txt`):

```text
1
0
1
1
0
````

Interpretation:

* `1` → presence, methylated, peak overlap, or “TRUE”
* `0` → absence, unmethylated, no peak, or “FALSE”

There is **no NA value** in format 0. Use Format 4 or 6 if you need NA or universe masking.

---

## 2. Packing to Format 0

### 2.1 From a raw 0/1 text vector

```bash
yame pack -fb binary_data.txt > binary_output.cg
```

`-fb` selects Format 0 (`b`).

---

### 2.2 From BED features or presence/absence annotations

This is the **canonical way** to create a `.cm`/`.cg` mask from BED:

```bash
bedtools sort -i dmr_sites.bed \
  | bedtools intersect -a cpg_ref.bed.gz -b - -sorted -c \
  | cut -f4 \
  | yame pack -fb - > dmr_binary.cg
```

Explanation:

* `bedtools intersect -c` counts how many times each CpG overlaps the BED input
* `cut -f4` extracts that count (0 or >0)
* `yame pack -fb` converts the result into a compressed CX file

If a CpG overlaps ≥ 1 region, the value becomes `1`; otherwise `0`.

---

## 3. Unpacking Format 0

To get the 0/1 vector back:

```bash
yame unpack binary_output.cg | head
```

Example output:

```text
1
0
1
1
0
```

Format 0 does **not** include sample-level metadata unless packed from multiple samples; YAME will unpack each sample sequentially if multiple samples exist.

---

## 4. Integration with Other YAME Commands

Format 0 integrates well with many downstream YAME tools.

---

### 4.1 `yame summary`

Format 0 is ideal for:

* Feature masks
* Peak overlaps
* Region presence indicators

Example:

```bash
yame summary -m promoter.cm sample.cg
```

Outputs for each mask:

* `N_mask`
* `N_overlap`
* `log2OddsRatio`
* Binary fraction (`Beta`)
* Universe counts depending on the query

---

### 4.2 `yame rowop`

Useful operations on binary data:

```bash
# Sum across samples (pseudobulk for binary)
yame rowop -o binasum multi_sample.cg > pseudobulk.cg

# Convert to binary string representation
yame rowop -o binstring multi_sample.cg > patterns.txt
```

`binasum` → format 3 output, with M = #1 votes and U = #0 votes.

`binstring` → one binary string per row, e.g.:

```
01011001
11100011
...
```

---

### 4.3 `yame dsample`

Downsample a binary file to N “present” sites:

```bash
yame dsample -N 10000 -s 1 binary.cg > dsampled.cg
```

For Format 0:

* Eligible sites are those with value `1`
* Non-selected sites become `0`

---

### 4.4 `yame rowsub`

Subset rows from Format 0 the same as other formats:

```bash
# By row IDs
yame rowsub -l row_ids.txt binary.cg > subset.cg

# By coordinate list and reference
yame rowsub -R cpg_nocontig.cr -L CpG_sites.txt binary.cg > subset.cg

# By mask
yame rowsub -m promoter.cm binary.cg > subset.promoters.cg
```

---

### 4.5 `yame mask`

Format 0 plays especially well with `yame mask`:

```bash
# Mask out positions where mask == 1
yame mask binary.cg lowquality.cm -o masked_binary.cg
```

If you use `-c`, Format 0 can be **contextualized** into Format 6 to define a universe for sparse annotations.

---

## 5. When to Use Format 0 vs Other Formats

Use **Format 0** when:

* Your data is intrinsically binary (e.g., peak/no-peak)
* You want minimal storage footprint
* You want maximum speed for summarization / enrichment
* You are constructing feature files (promoters, enhancers, windows)

Consider alternatives if:

* You need M/U counts → **Format 3**
* You need NA handling or fractions → **Format 4**
* You need structured query + universe semantics → **Format 6**

---

## 6. Minimal End-to-End Example

```bash
# 1. Prepare CpG reference (only once)
yame unpack cpg_nocontig.cr | gzip > cpg_ref.bed.gz

# 2. Create binary mask from BED peaks
bedtools intersect -a cpg_ref.bed.gz -b H3K27ac.bed -sorted -c \
  | cut -f4 \
  | yame pack -fb - > H3K27ac.cm

# 3. Summarize enrichment of a methylation sample over H3K27ac peaks
yame summary -m H3K27ac.cm sample.cg > enrich.txt

# 4. Subset sample to only CpGs in promoter regions
yame rowsub -m promoters.cm sample.cg > sample.promoters.cg

# 5. Downsample the promoter mask to 5000 sites
yame dsample -N 5000 promoters.cm > promoters_5k.cm
```

