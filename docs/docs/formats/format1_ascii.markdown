---
title: Format 1 – ASCII Values with RLE Compression
parent: 1. Storage & Format
nav_order: 1
---

# Format 1 – ASCII Values with RLE (Run-Length Encoding)

Format 1 stores **single-byte ASCII values** for each CpG, compressed using **run-length encoding (RLE)**.  
It is the most flexible CX format for encoding *arbitrary categorical or numeric symbols* that can be represented as **one ASCII character** (0–255).

Examples of valid values:

- ASCII digits: `0`, `1`, `2`, …  
- ASCII letters: `a`, `b`, `A`, `Z`  
- Punctuation or symbolic codes: `.`, `*`, `+`, `?`  
- Encoded or discretized values where each state is a **single-byte token**

This format is appropriate whenever the data is:

- **Discrete**, but not limited to binary (Format 0)  
- **Symbolic**, but not needing full string labels (Format 2)  
- **Piecewise constant**, allowing large compression via RLE  

---

## 1. When to Use Format 1

Use **Format 1** when your per-CpG track:

- Represents **one ASCII symbol per CpG**
- Uses **small alphabets**, leading to long repeated runs  
  (e.g., most values are `'0'`, `'1'`, `'2'`, or `'.'`)
- Must remain symbolic rather than numeric
- Must be efficiently stored and accessed  
- Does NOT require:
  - floating-point values (use Format 4)
  - multi-byte integers (Format 1 only stores **1 byte**)
  - arbitrary text labels (use Format 2)
  - pairs of values (use Format 3)

Typical applications:

- Discretized methylation states (`0`, `1`, `2`, `3`)
- Coverage tiers (`A`, `B`, `C`, `D`)
- Peak strength encoded as ASCII values  
- Custom one-character categorical encodings  
- Low-memory symbolic genome segmentation

---

## 2. Input Requirements

Format 1 accepts **any single ASCII character** per line.

Valid examples:

```text
0
1
2
a
b
A
.
*
?
````

Multi-character strings like `"10"` or `"state1"` are **not allowed** — use Format 2 if you need labels longer than one byte.

Missing values:

* There is **no NA in Format 1**
* You may encode missing as `'.'` or `'0'`, depending on your design

---

## 3. Packing to Format 1

### 3.1 Pack from a simple ASCII list

```bash
yame pack -f1 ascii_values.txt > ascii_track.cg
```

Where `ascii_values.txt` contains one character per line.

### 3.2 Convert numeric or discrete BED values to ASCII

Example: discretizing peak intensities into ASCII symbols:

```bash
bedtools intersect -a cpg_ref.bed.gz -b peaks.bed -loj -sorted \
  | awk '{ if ($8==".") print "."; else if ($8<5) print "a"; else print "b"; }' \
  | yame pack -f1 - > peak_states.cg
```

### 3.3 Overlap count → ASCII (integer → ASCII conversion)

```bash
bedtools intersect -a cpg_ref.bed.gz -b regions.bed -sorted -c \
  | cut -f4 \
  | awk '{ printf("%c\n", $1); }' \
  | yame pack -f1 - > counted.cg
```

Any number **0–255** may be converted to ASCII using `printf("%c")`.

---

## 4. Unpacking Format 1

To recover ASCII characters:

```bash
yame unpack ascii_track.cg | head
```

Output is exactly one character per line:

```text
0
0
0
1
1
.
.
a
a
2
```

It is **lossless**, even for non-printable bytes.

---

## 5. Integration with Other YAME Commands

---

### 5.1 `yame summary`

```bash
yame summary -m feature.cm ascii_track.cg
```

Interprets ASCII data as **integer values 0–255** internally when computing summary statistics:

* `Beta` = mean of ASCII byte values
* `N_query` = number of CpGs where value ≠ 0
* `Log2OddsRatio` treats nonzero values as "hits"

If your encoding uses characters where `'0'` is not the null value, you may want to preprocess.

---

### 5.2 `yame rowsub`

Supports all selection modes:

```bash
yame rowsub -m promoters.cm ascii_track.cg > promoters_ascii.cg
yame rowsub -R cpg.cr -L CpGs.txt ascii_track.cg > subset.cg
```

---

### 5.3 `yame dsample`

Downsampling uses the rule:

* **Non-zero ASCII values** = eligible
* Retain N
* Others become `'0'` (ASCII zero)

Example:

```bash
yame dsample -N 50000 ascii_track.cg > ascii_50k.cg
```

---

### 5.4 `yame mask`

Masking replaces values with `'0'`:

```bash
yame mask ascii_track.cg mask.cm -o masked.cg
```

---

## 6. When NOT to Use Format 1

Avoid Format 1 if:

| Requirement                | Use Instead  |
| -------------------------- | ------------ |
| Values longer than 1 byte  | **Format 2** |
| Floating-point, NA support | **Format 4** |
| M/U count pairs            | **Format 3** |
| Sparse universe + binary   | **Format 6** |
| Only 0/1 values            | **Format 0** |

Format 1 is best for **compact single-byte symbolic encodings**.

---

## 7. Minimal End-to-End Example

```bash
# 1. Prepare reference CpGs
yame unpack cpg_nocontig.cr | gzip > cpg_ref.bed.gz

# 2. Encode peak intensity categories as ASCII
bedtools intersect -a cpg_ref.bed.gz -b peaks.bed -loj -sorted \
  | awk '{ if ($8==".") print "."; else if ($8<5) print "a"; else print "b"; }' \
  | yame pack -f1 - > peaks_ascii.cg

# 3. Summarize across annotations
yame summary -m ChromHMM.cm peaks_ascii.cg > peak_state_enrichment.txt

# 4. Extract promoter subset
yame rowsub -m promoters.cm peaks_ascii.cg > promoter_ascii.cg

# 5. Downsample symbolic data
yame dsample -N 10000 peaks_ascii.cg > peaks_10k.cg
```

---

Format 1 provides an efficient, flexible container for **single-byte, symbolic per-CpG encodings**, bridging the gap between simple binary formats and complex state-based annotations.

