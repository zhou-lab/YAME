---
title: 3. Visualization with hprint
nav_order: 3
---

# 3. Transposed Visualization with `yame hprint`

`yame hprint` is a powerful tool for visual inspection of methylation data. It prints the data "horizontally" (transposed), where samples are rows and genomic sites are columns. 

It supports two main modes:
- **Region View**: Detailed visualization of a specific genomic interval.
- **Whole-Genome View**: High-level summary of methylation across all chromosomes.

---

## 3.1 Basic Usage

```bash
yame hprint [options] -R <ref.cr> <in.cx>
```

- `-R <ref.cr>`: (Required for genomic views) The format 7 reference file used to create the data.
- `<in.cx>`: The input methylation file (supports Format 0, 3, 4, 6, and 7).

### 3.1.1 Color Output

By default, `yame hprint` uses ANSI escape codes to provide a rich color visualization. 
- To **disable** color (e.g., for piping to text-processing tools), use the `-c` flag.

### 3.1.2 Controlling Width and Windowing

If the region or genome being viewed is wider than the terminal (default: 80 columns), `yame hprint` automatically performs **on-the-fly averaging**. It calculates the mean beta value for each window of sites to fit the requested width.

- `-w <int>`: Set the maximum number of data columns (default: 80).
- `-l <int>`: Set the width of the sample label column (default: 20).
- `-t <int>`: Ruler tick interval in columns (default: 10).

---

## 3.2 View Modes

### 3.2.1 Region View (`-r`)

Visualize a specific chromosome or sub-region:

```bash
# View an entire chromosome
yame hprint -R hg38.cr -r chr16 data.cx

# View a specific coordinate range (1-based)
yame hprint -R hg38.cr -r chr16:10,000,000-10,100,000 data.cx
```

In region view, the ruler shows coordinates and the number of CpGs included. If windowing is active, it also shows the window size (`win=N`).

### 3.2.2 Whole-Genome View

If you provide `-R` but omit `-r`, `yame hprint` shows the entire dataset compressed into the requested width (default 80 columns). This is excellent for identifying large-scale methylation patterns or chromosome-wide changes.

---

## 3.3 Granular Decile Mode (`-g`)

By default, `yame hprint` uses three symbols to represent methylation levels:
- `H`: High beta (> 0.67)
- `M`: Medium beta (0.33 - 0.67)
- `L`: Low beta (< 0.33)
- `.`: Missing / No coverage

For more detail, use the **granular mode** (`-g`):

```bash
yame hprint -g -R hg38.cr -r chr16 data.cx
```

In granular mode:
- Beta values are mapped to **10 deciles** represented by digits `0` to `9`.
- **0**: 0.0 - 0.1
- **5**: 0.5 - 0.6
- **9**: 0.9 - 1.0

### 3.3.1 Blue-White-Red Spectrum

When granular mode and color are active, `yame hprint` uses a high-quality diverging color map:
- **Blue (0-4)**: Progresses from deep blue to light blue (unmethylated).
- **White (5)**: Intermediate methylation.
- **Red (6-9)**: Progresses from light red to deep red (methylated).

---

## 3.4 Technical Implementation

`yame hprint` is optimized for high performance and low memory footprint:

- **Two-Pass Scanning**: Genomic coordinates are gathered using a two-pass scan of the reference file, avoiding large memory allocations even for very large regions.
- **On-the-Fly Averaging**: Calculations are performed directly from source buffers using offsets.
- **Zero-Copy**: Data is processed in its original memory location without redundant slicing or copying.
- **Partial Decompression**: For Format 3 (`.cg`), only the required region is decompressed, making it efficient for large files.
