---
title: Fmt 6 – Set & Universe
parent: 1. Storage & Format
nav_order: 6
---

# Format 6 – Query & Universe (Boolean Sparse Representation)

Format 6 encodes two complementary boolean vectors for each CpG:

- **SET bit** → CpG is included in the *query*  
- **UNIVERSE bit** → CpG belongs to the background universe  

This format is essential for:

- **Enrichment testing** (query vs universe)
- **Sparse single-cell methylation** (accessible universe + binary call)
- **Binary presence/absence with explicit context**
- **Efficient representation of cell-specific CpG availability**

Format 6 allows YAME to know *both* the methylated/unmethylated call **and** whether the CpG is even *in scope* for that sample.

---

## 1. When to Use Format 6

Use Format 6 when:

- You want to run **enrichment tests** (`yame summary -m mask.cm sample.cx`)
- You want to represent **sparse CpG coverage** per sample (e.g., scWGBS)
- You want a **contextualized binary representation** in which:
  - CpG exists in the cell’s “universe”
  - CpG is methylated/unmethylated inside that universe

Format 6 is the only CX format that explicitly encodes **two bits per CpG**.

This makes it the natural fit for:

- Single-cell bisulfite data  
- Sparse pseudo-bulk models  
- Query/background analyses  
- Fine-grained masking workflows (e.g., `yame mask -c`)

---

## 2. Input Requirements

Input must contain **two columns**, each with 0 or 1:

```

<QueryBit>    <UniverseBit>

````

Example:

```text
1    1   # in query, in universe
0    1   # not in query, in universe
1    1   # in query, in universe
0    0   # outside the universe entirely
````

Meaning:

| Query | Universe | Interpretation                                      |
| ----- | -------- | --------------------------------------------------- |
| 1     | 1        | CpG is part of query set                            |
| 0     | 1        | CpG is in universe but not part of query            |
| 0     | 0        | CpG should be **ignored** entirely                  |
| 1     | 0        | **Invalid** (cannot be in query if not in universe) |

YAME enforces this logic automatically when packing/masking.

---

## 3. Packing to Format 6

### 3.1 Packing from a text file with two columns

```bash
yame pack -f6 query_universe.txt > sample.cx
```

---

### 3.2 Using BED inputs (query and universe tracks)

Generate a query mask:

```bash
bedtools intersect -a cpg_ref.bed.gz -b query_peaks.bed -sorted -c \
  | cut -f4 > query.txt
```

Generate a universe mask:

```bash
bedtools intersect -a cpg_ref.bed.gz -b regions_accessible.bed -sorted -c \
  | cut -f4 > universe.txt
```

Combine to two columns:

```bash
paste query.txt universe.txt | \
  awk '{print ($1>0?1:0) "\t" ($2>0?1:0)}' \
  | yame pack -f6 - > query.cx
```

---

### 3.3 From Format 0 + Universe mask using `yame mask -c` (recommended)

If you already have a binary vector (Format 0), and a universe mask:

```bash
# binary input: 1 = methylated, 0 = unmethylated
yame mask -c input_binary.cx universe_mask.cx -o contextual.cx
```

Rules:

* Universe = mask’s 1s
* Query = input’s binary values inside universe

This is common in single-cell pipelines.

---

## 4. Unpacking Format 6

```bash
yame unpack sample.cx | head
```

The output has two columns:

```text
<Query>    <Universe>
1          1
0          1
1          1
0          0
...
```

Matching the original input format.

---

## 5. Integration with YAME Commands

Format 6 integrates deeply with the core statistical logic of YAME.

---

### 5.1 `yame summary`

Format 6 is designed for enrichment analysis:

```bash
yame summary -m feature.cm sample.cx
```

Interpretation:

* `N_univ` = number of CpGs where `UniverseBit = 1`
* `N_query` = number of CpGs with both `UniverseBit = 1` **and** `QueryBit = 1`
* `N_mask` = size of feature
* `N_overlap` = how many CpGs are (in query) AND (in mask)
* `Log2OddsRatio` = enrichment of query inside mask
* `Beta` = fraction of `QueryBit = 1` among universe CpGs within mask
* `Depth` = NA (Format 6 has no coverage)

This is the recommended representation for:

* Comparing different CpG subsets
* Representing sparse single-cell methylomes
* Querying binary signatures within genomic features

---

### 5.2 `yame rowop`

Useful operations:

```bash
# Binary summation across samples
yame rowop -o binasum multi.cx > pseudobulk.cg

# Convert to binary string
yame rowop -o binstring multi.cx > patterns.txt
```

Rules:

* Query bit → methylated/unmethylated value
* Universe bit → determines coverage; if Universe = 0, skip site entirely

---

### 5.3 `yame dsample`

Downsampling Format 6:

```bash
yame dsample -N 10000 sample.cx > sample_10k.cx
```

Rules:

* Eligible sites = UniverseBit = 1
* Randomly keep **N** of them per sample
* Non-selected → UniverseBit set to **0**; QueryBit cleared

This is extremely useful for benchmarking single-cell sparsity.

---

### 5.4 `yame rowsub`

Format 6 supports all selection methods:

```bash
# Subset by mask
yame rowsub -m promoters.cm sample.cx > subset.cx

# Coordinate lists
yame rowsub -R cpg_nocontig.cr -L CpG_sites.txt sample.cx > subset.cx
```

---

### 5.5 `yame mask` (contextualization mode)

```bash
# Convert binary format 0 → format 6 with a universe mask
yame mask -c binary.cx universe.cx -o out.cx
```

* Universe bits taken from mask
* Query bits taken from binary vector
* Sites outside universe become NA-like (Universe = 0)

This workflow is central to building high-quality **single-cell methylome objects**.

---

## 6. When NOT to Use Format 6

Use a different format if:

* You need floating-point methylation values → **Format 4**
* You need M/U counts → **Format 3**
* You have categorical labels → **Format 2**
* You need simple 0/1 binary without context → **Format 0**

Format 6 is best when **context matters** (universe vs query).

---

## 7. Minimal End-to-End Example

```bash
# 1. Prepare reference CpGs
yame unpack cpg_nocontig.cr | gzip > cpg_ref.bed.gz

# 2. Build a universe mask (e.g., accessible CpGs)
bedtools intersect -a cpg_ref.bed.gz -b ATAC.bed -sorted -c \
  | cut -f4 \
  | yame pack -fb - > universe.cx

# 3. Build a binary methylation call track
bedtools intersect -a cpg_ref.bed.gz -b methylated_calls.bed -sorted -c \
  | cut -f4 \
  | yame pack -fb - > binary.cx

# 4. Contextualize into Format 6
yame mask -c binary.cx universe.cx -o cell.cx

# 5. Run enrichment against ChromHMM
yame summary -m ChromHMM.cm cell.cx > enrich.txt

# 6. Downsample to 50k accessible CpGs
yame dsample -N 50000 cell.cx > cell_50k.cx
```

---

Format 6 is the **essential representation for sparse binary methylation**, especially in *single-cell* analyses and *feature-level enrichment testing*.
It is compact, expressive, and foundational for high-performance YAME workflows.

```
