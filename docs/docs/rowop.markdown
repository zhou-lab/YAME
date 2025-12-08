---
title: 5. Aggregate Row-wise
nav_order: 5
---

# 5. Row-wise Aggregation of Methylation Sites

`yame rowop` provides fast, low-level operations applied **row-by-row** across all samples of a `.cx` file.  
These operations are frequently used in single-cell methylation analysis, pseudobulk aggregation, QC, and exploratory co-methylation analyses.

```

Usage: yame rowop -o <operation> <in.cx> <out>

````

If `<out>` is omitted, results are written to stdout.

---

# 5.1 Common Use Case: Summing Binary Data to Create Pseudobulks

A very common operation is to **sum binarized methylation calls** across cells belonging to the same cluster.

For example, to generate a pseudobulk for *cluster 1* from a single-cell `.cg` file:

```bash
yame subset -l cluster_1_id.txt single_cell.cg \
  | yame rowop -o binasum - > cluster1_pseudobulk.cg
````

What `binasum` does:

* For **binary formats (0/1)**:
  It counts how many cells have methylation = 1 vs 0.
* For **format 3 (M/U counts)**:
  It compares M vs U in each cell and votes for the majority state.
  Only rows with coverage ≥ `-c mincov` are counted.

The output is a **format 3** `.cg` file where each position contains:

* Total number of methylated votes (M)
* Total number of unmethylated votes (U)

This is ideal for downstream:

* Pseudobulk methylation beta estimates
* Differential methylation analysis
* Cluster-level QC

---

# 5.2 Available Operations

Below is a complete list of operations supported by `-o`.

---

## **1. `binasum` — Sum binarized methylation across samples**

```
yame rowop -o binasum input.cx output.cx
```

* Works on formats **0, 1, and 3**
* Produces a new **format 3** dataset
* For format 3 inputs, rows with depth < `-c` (default 1) are ignored
* Interprets each sample's methylation as a **binary vote**

Useful for:

* Pseudobulk construction
* Collapsing replicates
* Aggregating barcodes / technical partitions

---

## **2. `musum` — Sum M and U counts directly**

```
yame rowop -o musum input.cx output.cx
```

* Works only on **format 3**
* Directly adds up:

  * all M counts
  * all U counts
* Produces a new `.cx` file in format 3

Useful when:

* You want true aggregated read counts
* You trust the original M/U values (not just binary interpretation)

---

## **3. `mean` — Per-row methylation mean and count**

```
yame rowop -o mean input.cx > mean.tsv
```

Output (tab-delimited):

```
beta_mean    number_of_valid_cells
```

Only samples with coverage ≥ `-c` contribute.
Rows with no valid values output `NA    0`.

Useful for:

* QC summaries
* Computing single-cell methylation averages
* Identifying globally variable CpGs

---

## **4. `std` — Per-row methylation standard deviation**

```
yame rowop -o std input.cx > std.tsv
```

Output columns:

```
std_dev    counts
```

Uses the same filtering as `mean`.
Helpful for:

* Measuring variability
* Ranking variable CpGs
* Identifying DNA methylation landmarks

---

## **5. `binstring` — Convert methylation profiles to binary strings**

```
yame rowop -o binstring -b 0.6 input.cx > binstrings.txt
```

* Converts each sample to a binary string per row
* Threshold `-b` decides methylated vs unmethylated (default 0.5)
* Outputs one string **per row**

Example output:

```
0010110101
1100001110
...
```

Useful for:

* Clustering by Hamming distance
* Haplotype inference
* Co-methylation motif discovery
* Compression for machine learning models

---

## **6. `cometh` — Co-methylation of neighboring CpGs**

```
yame rowop -o cometh -w 5 input.cx > cometh.tsv
```

For each CpG *i* and neighbors *i+1 … i+w*, `cometh` outputs a vector encoding four categories:

* U0U1
* U0M1
* M0U1
* M0M1

Counts are stored in a compact 64-bit representation; use `-v` to print them explicitly:

```
i   U0U1-U0M1-M0U1-M0M1   U0U2-U0M2-M0U2-M0M2   ...
```

Additional behavior:

* Only format 3 is allowed
* Excludes intermediate methylation values (0.3–0.7)
* Requires both CpGs to have depth ≥ `-c mincov`
* `-w` sets window size (default: 5)

This is useful for:

* Detecting locally coordinated methylation
* Assessing haplotype methylation patterns
* Fragment-level epiallele analysis
* Identifying footprints or nucleosome-scale structure

---

# 5.3 Summary of Operations

| Operation   | Output Type      | Input Requirement | Purpose                       |
| ----------- | ---------------- | ----------------- | ----------------------------- |
| `binasum`   | `.cx` (format 3) | fmt 0/1/3         | Pseudobulk / aggregation      |
| `musum`     | `.cx` (format 3) | fmt 3             | True count summation          |
| `mean`      | text             | fmt 3             | Mean methylation per CpG      |
| `std`       | text             | fmt 3             | Standard deviation per CpG    |
| `binstring` | text             | fmt 3             | Binary methylation strings    |
| `cometh`    | text             | fmt 3             | Local co-methylation patterns |

---

# 5.4 Additional Notes

### Minimum coverage control (`-c`)

Many operations ignore rows where:

```
M + U < mincov
```

Default is 1.

### Verbose output (`-v`)

For `cometh`, verbose mode expands packed category counts into readable numbers.

### Binarization threshold (`-b`)

Relevant only for `binstring`.

---

# 5.5 Help and Subcommand Documentation

For detailed usage:

```bash
yame rowop -h
```

Subcommand documentation:

* [**rowop help page**]({% link docs/subcommands/YAME_rowop.markdown %})

