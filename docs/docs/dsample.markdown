---
title: 7. Mask Data
nav_order: 7
---

# 7. Downsampling & Masking Methylation Sites

YAME provides tools to **control sparsity** and **apply masks** to methylation data:

- `yame dsample` — randomly downsample sites to simulate lower coverage or sparsity.
- `yame mask` — apply a binary mask to zero out sites or convert a binary format into a contextualized format 6.

These functions are especially useful for benchmarking methods at different sparsity levels, building controlled simulation datasets, or restricting analyses to a specific universe of CpGs.

---

## 7.1 Random Downsampling with `yame dsample`

`yame dsample` randomly keeps a fixed number of non-NA sites per sample and masks out the rest.

Supported input formats:

- **Format 3** (`M/U` counts)
- **Format 6** (binary with universe bit; often used for single-cell sparse data)

Basic usage:

```bash
yame dsample -N 10000 -s 1 input.cg > downsampled.cg
````

This keeps **10,000** covered sites per sample (or all if fewer are available), using seed `1` for reproducibility.

### 7.1.1 What `dsample` Does (Format 3 vs Format 6)

* **Format 3 (`.cg`, M/U counts)**

  * Eligible sites are those with `M+U > 0`.
  * `dsample` randomly selects **N** such sites.
  * **Selected sites** keep their original `M` and `U`.
  * **Non-selected sites** are masked by setting `M=U=0` (treated as missing/NA).

* **Format 6 (universe-bit binary)**

  * Eligible sites are those in the **universe** (`FMT6_IN_UNI`).
  * `dsample` randomly selects **N** universe positions.
  * **Selected sites** remain unchanged.
  * **Non-selected sites** have their universe bit cleared (`FMT6_SET_NA`), effectively dropping them from the analyzable universe.

In both cases, if `N` is larger than the number of eligible sites, **all eligible sites are kept** (no error).

### 7.1.2 Key Options

```bash
yame dsample [options] <in.cx> [out.cx]
```

Options:

* `-N [int]` — number of eligible sites to keep **per sample** (default: `100`).
* `-s [int]` — random seed (default: current time). Use a fixed seed for reproducible downsampling.
* `-r [int]` — number of **independent replicates per sample** (default: `1`).

  * Each replicate is downsampled separately from the same input sample.
* `-h` — show help.

Output destination:

* If `out.cx` (positional) or `-o` is provided: write to that file and also write an index.
* If no output is given: write to **stdout** (no index).

### 7.1.3 Replicates and Index Naming

When `-r` is greater than 1, `dsample` creates multiple downsampled versions of each input sample.

* If the input `.cx` has an index:

  * The original sample names are used as a **base**.
  * Replicates are suffixed: `SampleA-0`, `SampleA-1`, ..., `SampleA-(r-1)`.
* If no input index exists:

  * Samples are named `0`, `1`, `2`, ... internally.
  * Replicates follow the same `base-rep` naming pattern.

Example: create **5** replicates with 50k sites each:

```bash
yame dsample -N 50000 -r 5 input.cg downsampled.cg
```

The resulting index in `downsampled.cg.idx` will contain entries like:

```
Sample1-0
Sample1-1
...
Sample1-4
Sample2-0
...
```

### 7.1.4 Typical Use Cases

* **Benchmarking methods at different sparsity levels**
  Run the same pipeline on `N = 1e5`, `N = 5e4`, `N = 1e4` to see robustness to coverage.

* **Generating multiple randomized sparsity replicates**
  For each sample, simulate multiple downsampling replicates with different seeds or with `-r`.

* **Single-cell simulations with format 6**
  Use `dsample` to progressively restrict the universe of accessible CpGs and observe performance changes.

For more help with `dsample`, run:

```bash
yame dsample -h
```

or see the
[**dsample help page**]({% link docs/subcommands/YAME_dsample.markdown %}).

---

## 7.2 Masking and Contextualization with `yame mask`

`yame mask` applies a **row-wise mask** to a `.cx` file and optionally converts binary data into **format 6** for contextualized single-cell usage.

Basic usage:

```bash
yame mask input.cg mask.cx -o masked.cg
```

Here:

* `input.cg` — query methylation file (format 0, 1, or 3).
* `mask.cx` — mask file (format 0, 1, or 3; internally converted to a binary format 0).
* Output `masked.cg` contains only the unmasked positions.

### 7.2.1 Supported Inputs and Mask Semantics

* **Mask file (`mask.cx`)**:

  * Can be format 0, 1, or 3.
  * Format 1 and 3 masks are **converted to format 0**:

    * For format 1: `1` is treated as masked, `0` unmasked.
    * For format 3: sites with `M+U > 0` become `1` (masked), zeros become `0`.
* **Query file (`input.cg`)**:

  * Format 3: M/U counts.
  * Format 0/1: binary.

The mask must have the **same row length** as the query; otherwise the command will abort with an error.

By default, bits that are `1` in the mask are **masked out** (removed).

If you set `-v`, the mask is **inverted**, so bits that are `0` in the original mask become masked.

### 7.2.2 Operations Without Contextualization (default)

```bash
yame mask input.cg mask.cx -o masked.cg
```

* If the query is **format 3**:

  * For every site where the mask bit is `1`, set `M=U=0`.
  * Effect: those sites are treated as missing.

* If the query is **format 0**:

  * Perform a binary AND with the complement of the mask (`c &= ~mask`).
  * Effect: all `1`s in the mask are forced to `0` in the query.

This is handy for:

* Removing blacklist CpGs from an existing `.cg`.
* Removing low-quality or low-coverage sites.
* Restricting analysis to a curated panel of CpGs.

### 7.2.3 Contextualizing to Format 6 (`-c`)

With `-c`, `yame mask` turns a **binary query** plus a **mask** into a **format 6** object:

```bash
yame mask -c input_binary.cx mask.cx -o contextualized.cx
```

Behavior:

* The **mask** defines the **universe** (sites that exist in the cell).
* The **binary values** in the query define whether each universe site is methylated (`1`) or unmethylated (`0`):

  * If `mask[i] = 1`:

    * If `input[i] = 1` ➜ set format 6 as methylated (`FMT6_SET1`).
    * If `input[i] = 0` ➜ set format 6 as unmethylated (`FMT6_SET0`).
  * If `mask[i] = 0`:

    * Site is outside the universe (no entry in the resulting format 6 vector).

With `-v`, the universe is effectively the complement of the mask (invert mask before contextualizing).

This is useful for:

* Defining **cell- or experiment-specific universes** while retaining 0/1 methylation calls.
* Converting feature masks into sparse, contextualized single-cell objects.

### 7.2.4 Command Summary

```bash
yame mask [options] <in.cx> <mask.cx>
```

Options:

* `-o [PATH]` — output `.cx` file name. If missing, write to stdout (no index).
* `-c` — contextualize to format 6 using `1`s in mask as the universe.
* `-v` — invert the mask (mask `0`s instead of `1`s).
* `-h` — help.

Example workflows:

**Mask out low-quality sites:**

```bash
yame mask highcov.cg lowqual_mask.cx -o highcov_masked.cg
```

**Restrict to a predefined universe and create fmt6:**

```bash
yame mask -c raw_binary.cx universe_mask.cx -o cell_fmt6.cx
```

For more help with `mask`, run:

```bash
yame mask -h
```

or see the
[**mask help page**]({% link docs/subcommands/YAME_mask.markdown %}).

