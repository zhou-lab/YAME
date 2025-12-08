---
title: 6. Combine, Split & Index
nav_order: 6
---

# 6. Combine, Split & Index

This section covers operations that manipulate **samples** (columns) in a `.cx` file:

- **`yame index`** — build or update a sample index  
- **`yame split`** — split a multi-sample `.cx` file into individual `.cx` files  
- **`yame subset`** — extract selected samples or states  
- **Combining files** — achieved using standard Unix `cat`  

These tools are essential for workflows that involve multiple `.cx` samples, such as merging epigenomic features, analyzing groups of samples, or reorganizing sample structure.

---

# 6.1 Generating and Updating Index Files (`yame index`)

A `.cx` file containing multiple samples stores its samples **sequentially**, and YAME uses an accompanying index file (`.cx.idx`) to record the byte offset of each sample.

You can generate an index using:

```bash
yame index yourfile.cx
````

This produces:

```
yourfile.cx.idx
```

with two columns:

```
sample_name    byte_offset
```

## Assigning Sample Names from a List

If the `.cx` has *N* samples but no index file, provide a sample-name list:

```bash
yame index -s sample_names.tsv yourfile.cx
```

`sample_names.tsv` contains names in its first column.

## Appending a New Sample (`-1`)

If you have added an extra `cdata` block to the end of a `.cx` file, you may append it to the existing index:

```bash
yame index -1 NewSampleName yourfile.cx
```

YAME locates the final block and records the offset of the newly appended sample.

## Output to Console

```bash
yame index -c yourfile.cx
```

This prints the index to stdout instead of writing `yourfile.cx.idx`.

---

# 6.2 Example: Merging Multiple `.cm` Feature Files

This example demonstrates how to:

1. Convert many BED files into `.cm` feature files
2. QC them
3. Merge them together
4. Generate a combined index

### **Step 1 — Prepare the table describing sample IDs and BED paths**

```
268     GSM648494       human_hm/268_sort_peaks.narrowPeak.bed
269     GSM648495       human_hm/269_sort_peaks.narrowPeak.bed
272     GSM575295       human_hm/272_b_sort_peaks.broadPeak.bed
...
```

### **Step 2 — Convert each BED into `.cm`**

```bash
cat controlfiles.tsv \
  | parallel --colsep '\t' -j 72 '
      id={1}; path={3};
      sortbed $path \
        | bedtools intersect -a cpg_nocontig.bed.gz -b - -sorted -c \
        | cut -f4 \
        | yame pack -f b - $id.cm
    '
```

Each `$id.cm` is a binary feature mask aligned to the CpG coordinate list.

### **Step 3 — QC each `.cm`**

```bash
awk '{print ""$1".cm", $2";"$4;}' controlfiles.tsv \
  | while read fn anno; do yame summary $fn; done \
  > qc.txt
```

Example filter: keep feature files with ≥ 5000 overlapping CpGs.

### **Step 4 — Merge and index**

```bash
awk '$1!~/QFile/ && $6>5000' qc.txt \
  | awk 'NR==FNR{a[$1]=1;}NR!=FNR&&($1".cm" in a){print $0;}' - controlfiles.tsv \
  | awk '{print ""$1".cm", $2";"$4;}' \
  | sort -k2,2 \
  | while read fn anno; do
        cat $fn >> merged.cm
        yame index -1 $anno merged.cm
    done
```

`merged.cm` is the concatenation of all retained `.cm` samples, with indexing updated each iteration.

---

# 6.3 Splitting Multi-Sample Files (`yame split`)

`yame split` takes a **multi-sample** `.cx` file and produces one `.cx` file per sample.

Basic usage:

```bash
yame split input.cx output_prefix
```

If sample names are present in the index, the output naming scheme becomes:

```
output_prefix<SampleName>.cx
```

Otherwise:

```
output_prefix_split_1.cx
output_prefix_split_2.cx
...
```

## Providing a Sample List

If the `.cx` has no index file but you know sample names:

```bash
yame split -s sample_list.txt input.cx prefix_
```

`sample_list.txt` should contain one name per line.

This preserves sample naming and ensures the `prefix_<sample>.cx` files correspond correctly.

For more help:

* [**split help page**]({% link docs/subcommands/YAME_split.markdown %})

---

# 6.4 Subsetting Samples (`yame subset`)

`yame subset` extracts a subset of samples from a multi-sample `.cx` file.
It uses the `.cx.idx` file to locate and extract the requested samples efficiently.

Basic syntax:

```bash
yame subset -l sample_list.txt input.cx > subset.cx
```

or:

```bash
yame subset input.cx SampleA SampleB SampleC > subset.cx
```

If you specify an output file via `-o`, YAME writes both:

* the new subset `.cx`
* a new index `.cx.idx`

Example:

```bash
yame subset -o cluster1.cx singlecell.cx Cell_01 Cell_07 Cell_33
```

### **Head / Tail extraction**

Useful for inspecting the first or last *N* samples:

```bash
yame subset -H 10 input.cx > first10.cx     # first 10 samples
yame subset -T 5 input.cx > last5.cx        # last 5 samples
```

### **Subsetting Format 2 states (`-s`)**

If the input is a **format 2** `.cx` file (categorical states), you may split states into binary masks:

```bash
yame subset -s -l state_list.txt -o states.cx chromatin_states.cx
```

This produces one binary vector per selected state.

For more help:

* [**subset help page**]({% link docs/subcommands/YAME_subset.markdown %})

---

# 6.5 Combining `.cx` Files

YAME does **not** provide a dedicated `combine` command because combining `.cx` files is equivalent to concatenation:

```bash
cat sample1.cx sample2.cx sample3.cx > combined.cx
yame index combined.cx
```

Rules:

* All `.cx` files must have the **same format** and **same row dimension**.
* After combining, run `yame index` to regenerate the sample index.

This pattern is used in the `.cm` merging example above.

---

# Summary of Commands

| Command       | Purpose                                  | Notes                                     |
| ------------- | ---------------------------------------- | ----------------------------------------- |
| `yame index`  | Build or update a `.cx.idx` sample index | Required for fast sample lookup           |
| `yame split`  | Produce one `.cx` per sample             | Naming uses index or user-supplied list   |
| `yame subset` | Extract selected samples or states       | Supports head/tail and format-2 filtering |
| `cat`         | Combine `.cx` files                      | Must re-index after concatenation         |

