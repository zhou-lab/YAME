---
title: Fmt 1 - Integer Values
parent: 1. Storage & Format
nav_order: 1
---

# Format 1 – Integer Values with RLE Compression

Format 1 stores **integer per-CpG values** using **run-length encoding (RLE)** for efficient compression.  
It is ideal when your data consists of **small integers** repeated across long stretches of the genome.

Examples:

- Per-CpG **count values** (e.g., motif hits, coverage statistics)
- Fixed **categorical encodings** where numbers repeat often
- Processed or discretized tracks that map 1 integer → 1 CpG
- Anything that is not 0/1 (Format 0) and not M/U pairs (Format 3)

---

## 1. When to Use Format 1

Use **Format 1** when:

- You have **integer-valued** data (0, 1, 2, 3…)
- You expect **long runs of identical values**  
  (e.g., 0 for most CpGs, intermittent 1’s or 2’s)
- Storage size matters — Format 1 can give **10–100× compression** for RLE-friendly data
- You do *not* need:
  - multiple values per CpG (use Format 3)
  - floating-point values (use Format 4)
  - category labels (use Format 2)

---

## 2. Input Requirements

Input must:

- Contain **one integer per CpG**
- Be aligned to the reference coordinate file (`.cr`)
- Match the reference length exactly

Example input (`int_values.txt`):

```text
0
0
0
1
1
0
2
2
2
2
1
````

Valid values:

* Any non-negative integer (`0`, `1`, `2`, …)
* Larger integers are allowed but compress best when small and repeated

Missing values are **not supported** — use Format 4 if you need NA.

---

## 3. Packing to Format 1

### 3.1 Pack from a text file

```bash
yame pack -f1 int_values.txt > integer_track.cg
```

or using shorthand:

```bash
yame pack -f1 int_values.txt > mytrack.cg
```

---

### 3.2 Create integer tracks from BED (example: count overlapping regions)

```bash
bedtools intersect -a cpg_ref.bed.gz -b regions.bed -sorted -c \
  | cut -f4 \
  | yame pack -f1 - > regions_count.cg
```

Explanation:

* `bedtools intersect -c` → count overlaps for each CpG
* `cut -f4` → extract the count
* `yame pack -f1` → compress into format 1

If your BED contains multiple types of annotations (e.g., peak strengths), you can convert those into integers before packing.

---

## 4. Unpacking Format 1

To recover integer values:

```bash
yame unpack integer_track.cg | head
```

Output:

```text
0
0
0
1
1
0
2
2
2
2
1
```

Format 1 unpacks **exactly** what was packed (lossless RLE).

---

## 5. Integration with Other YAME Commands

Although Format 1 is less common than 0/2/3, it integrates with key YAME tools.

---

### 5.1 `yame summary`

You can compute summaries over masks:

```bash
yame summary -m features.cm integer_track.cg
```

Summary statistics for Format 1 include:

* `N_query` — number of CpGs with nonzero value
* `Beta` — **mean integer value** inside each mask
* `Depth` — unused (always NA)
* `Log2OddsRatio` — enrichment of nonzero entries in each mask

---

### 5.2 `yame rowsub`

Format 1 supports all row subsetting methods:

```bash
# By index
yame rowsub -l ids.txt integer_track.cg > subset.cg

# By coordinate list
yame rowsub -R cpg_nocontig.cr -L CpG_sites.txt integer_track.cg > subset.cg

# By mask
yame rowsub -m promoter_mask.cm integer_track.cg > promoter_subset.cg
```

---

### 5.3 `yame dsample`

Downsampling for Format 1 uses the rule:

* Eligible sites = *nonzero* integers
* Randomly keep `N` of them
* Remaining sites become `0`

Example:

```bash
yame dsample -N 20000 -s 123 int_track.cg > int_track_20k.cg
```

---

### 5.4 `yame mask`

You can mask out sites in Format 1 the same as Format 0 or 3:

```bash
yame mask int_track.cg low_quality.cm -o cleaned.cg
```

Masked CpGs become `0`.

---

## 6. When Format 1 Is Not Appropriate

Choose another format if:

* You have two values per CpG → use **Format 3**
* You need floating-point values → use **Format 4**
* You need NA handling → use **Format 4**
* You have categorical labels → use **Format 2**
* You need sparse universe-aware binary representation → use **Format 6**

---

## 7. Minimal End-to-End Example

```bash
# 1. Prepare reference CpGs (once per genome version)
yame unpack cpg_nocontig.cr | gzip > cpg_ref.bed.gz

# 2. Compute coverage (example: count overlaps)
bedtools intersect -a cpg_ref.bed.gz -b regions.bed -sorted -c \
  | cut -f4 \
  | yame pack -f1 - > regions_track.cg

# 3. Summarize over genomic features
yame summary -m ChromHMM.cm regions_track.cg > summary.txt

# 4. Subset to promoter CpGs
yame rowsub -m promoters.cm regions_track.cg > promoter_track.cg

# 5. Downsample high-intensity bins
yame dsample -N 5000 regions_track.cg > regions_5k.cg
```

---

Format 1 is a compact, efficient choice whenever your data is **integer-valued and run-length compressible**—a perfect middle ground between dense numeric matrices and coarse binary formats.

