---
title: 4. Subset Rows
nav_order: 4
---

# 4. Subsetting Rows in a Packed `.cx` File

YAME provides several tools for extracting subsets of rows (genomic sites) or breaking a `.cx` file into smaller pieces. The main tools covered here are:

- `yame rowsub` — extract specific CpG rows or ranges  
- `yame chunk` — divide a `.cx` file into size-based chunks  
- `yame chunkchar` — divide large text files into size-based chunks  

This section describes how to use each tool and the most common workflows.

---

# 4.1 Subset CpG Rows with `yame rowsub`

`yame rowsub` allows you to extract specific genomic rows from a `.cx` file.  
This is useful when:

- You want methylation values only at specific CpG sites  
- You want to subset a contiguous genomic block  
- You want to filter using a binary mask (`fmt0`)  
- You want to extract rows by genomic coordinate rather than integer index  

The output is written to **stdout**, allowing direct piping or saving to file.

---

## Basic Usage

```bash
yame rowsub [options] <in.cx> > subset.cx
````

The command supports **three major ways** to select rows:

---

## **A. Select rows by integer index (`-l`)**

If you have a list of row indices (1-based), e.g.:

```
12
455
9012
...
```

Then run:

```bash
yame rowsub -l row_ids.txt yourfile.cx > subset.cx
```

No sorting is required; YAME internally maintains the given order.

---

## **B. Select rows by genomic coordinate labels (`-L` + `-R`)**

This is the most common method for extracting specific CpG sites.

You will need:

1. A **row coordinate file** (`.cr`) containing genome coordinates for each `.cx` row
   YAME provides:

   * [mm10 row coordinates](https://github.com/zhou-lab/KYCGKB_mm10)
   * [hg38 row coordinates](https://github.com/zhou-lab/KYCGKB_hg38)

2. A **list of site names** in the format `chr_beg1`, one per line, e.g.:

```
chr16_18300002
chr16_18300046
chr16_18300140
chr16_18300162
chr16_18300172
```

Then run:

```bash
yame rowsub -R cpg_nocontig.cr -L CpG_sites.tsv yourfile.cx > subset.cx
```

YAME internally:

1. Loads the `.cr` coordinate file
2. Maps each `chrX_pos` string to the correct row index
3. Extracts only those rows
4. Outputs a smaller `.cx`

You can also include coordinates as the first dataset in the output with:

```bash
yame rowsub -1 -R cpg_nocontig.cr -L CpG_sites.tsv yourfile.cx > subset.cx
```

---

## **C. Select rows using a binary mask (`-m`)**

Mask file must be **format 0 or 1**, representing a vector of 0/1 flags.

Example:

```bash
yame rowsub -m mask.cx yourfile.cx > subset.cx
```

Any row where the mask equals `1` is retained.

This acts similarly to feature masking in `yame summary`.

---

## **D. Select a contiguous block (`-B` or `-I`)**

### Select a row range (0-based):

```bash
yame rowsub -B 1000_2000 yourfile.cx > subset.cx
```

Extracts rows **1000–1999**.

If only one number is given:

```bash
yame rowsub -B 1000 yourfile.cx
```

Outputs a single row.

### Select a block by block index (`-I`)

Useful for chunked batch processing:

```bash
yame rowsub -I 5_1000000 yourfile.cx
```

This extracts:

* Block index = 5
* Block size = 1,000,000 rows
* Row range = 5,000,000 to 5,999,999

If block size is omitted, default = 1,000,000.

---

# Summary: Ways to Select Rows

| Method                   | Option                  | When to Use                        |
| ------------------------ | ----------------------- | ---------------------------------- |
| Integer row list         | `-l <file>`             | You know row indices               |
| Genomic coords           | `-L <file> -R <row.cr>` | You know genome coordinates        |
| Binary mask              | `-m <mask.cx>`          | Filtering by precomputed selection |
| Range of rows            | `-B <beg_end>`          | Simple slicing by index            |
| Block slicing            | `-I <block_blockSize>`  | Batch processing                   |
| Include coords in output | `-1`                    | To append `.cr` for clarity        |

---

# 4.1.1 `rowsub` Common Examples

### Extract a single genomic region

```bash
yame rowsub -B 500000_510000 input.cx > region.cx
```

### Extract CpGs in a BED region

(using row coordinate file + interval expansion)

```bash
awk '{print $1"_"$2+1}' region.bed > coords.txt
yame rowsub -R genome.cr -L coords.txt input.cx > subset.cx
```

### Apply a mask file

```bash
yame rowsub -m mymask.cx input.cx > masked.cx
```

---

# 4.2 Chunking `.cx` Files with `yame chunk`

`yame chunk` splits a packed `.cx` file into **multiple smaller `.cx` files**, each containing a fixed number of rows.

This is essential for:

* Distributed computing
* Parallel model training
* Memory-efficient processing
* Splitting extremely large `.cx` files into manageable pieces

---

## Basic Usage

```bash
yame chunk -s <chunkSize> input.cx output_dir/
```

If `output_dir` is not provided, a directory named:

```
input.cx_chunks/
```

is automatically created.

---

## Example

Split into chunks of 500,000 rows each:

```bash
yame chunk -s 500000 input.cx chunks/
```

This produces:

```
chunks/0.cx
chunks/1.cx
chunks/2.cx
...
```

Each chunk file maintains the **same number of samples** as the original.
All samples are split identically row-wise.

---

# 4.3 Chunking Text Files with `yame chunkchar`

`yame chunkchar` works like `chunk`, but for plain text files rather than `.cx` files.
This is useful for splitting:

* BED files
* FASTA headers
* List files
* Any long line-based text file

---

## Basic Usage

```bash
yame chunkchar -s <chunkSize> input.txt
```

By default, output is written to:

```
input.txt_chunks/
```

---

## Example

Split a large text file into 1M-line chunks:

```bash
yame chunkchar -s 1000000 sites.txt
```

Outputs:

```
sites.txt_chunks/0.txt
sites.txt_chunks/1.txt
sites.txt_chunks/2.txt
...
```

Each file contains up to `chunkSize` lines.

---

# 4.4 Help and Developer References

For additional details:

* Run with `-h`

  ```bash
  yame rowsub -h
  yame chunk -h
  yame chunkchar -h
  ```
  
* See full subcommand documentation:

  * [**rowsub help page**]({% link docs/subcommands/YAME_rowsub.markdown %})
  * [**chunk help page**]({% link docs/subcommands/YAME_chunk.markdown %})
  * [**chunkchar help page**]({% link docs/subcommands/YAME_chunkchar.markdown %})

---

# 4.5 Summary Table

| Command     | Input | Output          | Purpose                                        |
| ----------- | ----- | --------------- | ---------------------------------------------- |
| `rowsub`    | `.cx` | `.cx` to stdout | Fine-grained row selection                     |
| `chunk`     | `.cx` | multiple `.cx`  | Split methylation matrix into fixed-size parts |
| `chunkchar` | text  | multiple `.txt` | Split large text files                         |

---
