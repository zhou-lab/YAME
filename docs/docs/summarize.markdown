---
title: 2. Summarize & Encode
nav_order: 2
---

# 2. Summarization of Packed `.cx` Files

`yame summary` and `yame info` provide quick, high-level summaries of `.cx` files created via `yame pack`. These tools help users inspect sample names, structural parameters, and feature-level statistics such as enrichment, methylation levels, and coverage.

To run a simple summary:

```bash
yame info yourfile.cg
yame summary yourfile.cg
````

Both single-sample and multi-sample `.cx` files are supported.

---

# Overview of What `yame summary` Computes

For each query sample, YAME reports:

* **N_univ** — total number of sites considered (universe)
* **N_query** — number of sites present or methylated in the query
* **N_mask** — number of sites in the mask
* **N_overlap** — intersection of query and mask
* **Log2OddsRatio** — an enrichment score
* **Beta** — average methylation or proportion (format dependent)
* **Depth** — approximate sequencing depth (only for formats with M/U counts)

Mask files may contain multiple masks; each is evaluated independently.

---

# Supported Input Formats

The query `.cx` file may be in any of the following internal YAME formats:

| Format    | Meaning                  | How Summary Interprets It              |
| --------- | ------------------------ | -------------------------------------- |
| **0 / 1** | Binary vectors           | Presence/absence                       |
| **2**     | State labels             | Multi-class state-specific summary     |
| **3**     | Methylation (M/U counts) | Beta, depth, overlap                   |
| **6**     | Binary with universe bit | Sparse methylation (e.g., single cell) |

Masks (`-m`) may also be in formats 0, 1, 2, or 6.

Format 7 (BED-like coordinates) is ignored by summary.

---

# Using Mask Feature Files

Masks allow aggregation of methylation or presence across functional features, e.g.:

* Fixed-size windows
* Chromatin states
* CpG islands, gene promoters
* Custom regions

Each mask is processed as a separate feature, producing a separate row in the summary output.

If no mask is provided, YAME produces a **global summary** of the entire dataset.

---

# Example: Summarization with a Window Mask File

Single-cell DNA methylomes are sparse, so window bins (e.g., 100 kb) are commonly used.

Download the mask file:

```bash
wget https://raw.githubusercontent.com/zhou-lab/KYCGKB_mm10/main/Win100k.20220228.cm
```

Run:

```bash
yame summary -m Win100k.20220228.cm single_cell.cg
```

Example output:

```
QFile            Query      MFile                   Mask        N_univ      N_query     N_mask      N_overlap   Log2OddsRatio   Beta    Depth
single_cell.cg   Sample1    Win100k.20220228.cm     chr1:30     21867837    1861715     589         48          -0.07           0.688   0.1
single_cell.cg   Sample1    Win100k.20220228.cm     chr1:31     21867837    1861715     574         36          -0.48           0.917   0.1
```

---

# Explanation of Output Columns

### **1. `QFile`**

The query `.cx` file.

### **2. `Query`**

Sample name from the sample index.
If missing, YAME assigns numerical IDs.

### **3. `MFile`**

The mask file used (or `"global"` when no mask is provided).

### **4. `Mask`**

The name of each mask:

* State label for format 2
* Sample name for binary masks
* Key names if `-T` is used

### **5. `N_univ`**

Number of valid universe sites, depending on format:

* Format 6 uses only universe-bit sites
* Other formats use total vector length

### **6. `N_query`**

Number of “present” positions in the query:

* Formats 0/1: value = 1
* Format 3: sites with nonzero M+U
* Format 6: universe-bit and set-bit sites

### **7. `N_mask`**

Number of positions included in the mask.

### **8. `N_overlap`**

Intersection of query and mask.

### **9. `Log2OddsRatio`**

Enrichment score computed from the 2×2 contingency table of query × mask membership.

### **10. `Beta`**

Average methylation or binary fraction inside the mask:

* Format 3: mean M/(M+U)
* Binary formats: fraction of “1”
* Format 6: fraction of SET sites within universe

### **11. `Depth`**

Average sequencing depth (M+U) across mask sites (format 3 only).

---

# Special Behaviors

### **Multiple Masks**

If the mask file contains multiple samples, each is processed independently.

### **State Masks (Format 2)**

Each state produces a distinct summary row.

### **Universe Subsetting (`-u`)**

Applies an additional universe mask to both query and features.

### **Memory Mode (`-M`)**

Loads all masks into RAM to minimize disk access.

### **Header Suppression (`-H`)**

Removes the header line for scripting convenience.

---

# Command Reference

```
Usage: yame summary [options] <query.cx>

Options:
  -m FILE     Mask feature (.cx) file. May contain multiple masks.
  -M          Load all masks into memory.
  -u FILE     Optional universe .cx file.
  -H          Suppress header output.
  -q NAME     Name to use when query file is '-'.
  -F          Use full file paths in output.
  -T          Always show section names (format 2).
  -s FILE     Sample list overriding the query index.
  -h          Display help.
```

---

# Additional Documentation

See also:

* [**Feature enrichment guide**]({% link docs/enrichment.markdown %})
* [**Subcommand documentation**]({% link docs/subcommands/YAME_summary.markdown %})
