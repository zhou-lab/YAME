---
title: Methylation Data Summarization
nav_order: 4
---

# Summarization of the packed .cx file

`yame summary` and `yame info` provide a brief summarization and basic parameters of the input .cx file. Assume you have a `.cg` file from `yame pack`, then you simply need to run the following to obtain the information.  

```bash
yame info yourfile.cg
```

```bash
yame summary yourfile.cg
```

`yame summary` also allows user to input a `-m` feature file, this is particularly useful if you want to aggregate methylation sites information over certain features and conduct enrichment testings. 

## An example of yame summary with mask feature file

For example, in single cell DNA methylome analysis, a typical approach is to merge the CG sites into a fixed window size bin because of high sparsity levels of the data. To do that, we can download the `Win100k.20220228.cm` file from [KYCG_mm10](https://github.com/zhou-lab/KYCGKB_mm10), and then run the following

```bash
yame summary -m Win100k.20220228.cm single_cell.cg
```

The ouput will look something like this 

```
QFile            Query      MFile                   Mask        N_univ      N_query     N_mask      N_overlap   Log2OddsRatio   Beta    Depth
single_cell.cg   Sample1    Win100k.20220228.cm     chr1:30     21867837    1861715     589         48          -0.07           0.688   0.1
single_cell.cg   Sample1    Win100k.20220228.cm     chr1:31     21867837    1861715     574         36          -0.48           0.917   0.1
```

The output columns are:
1. The file name of the query .cx
2. The sample name in the query.cx
3. The file name of the feature file
4. The name of each mask in the feature file
5. The number of CG sites covered in the universe
6. The number of CG sites covered in the query file
7. The number of CG sites covered in each mask
8. The number of CG sites covered both in the query file and the mask
9. Log2 of the odds ratio for the enrichment overlap
10. Average methylation levels of the query file in the mask
11. An approximation of the depth of coverage in the mask

For a full guideline on how to test feature enrichment using yame summary, check out [enrichment help]({% link docs/enrichment.markdown %}).

For more help with `summary`, run `yame summary` in the terminal or check out the
[pack help page]({% link docs/subcommands/YAME_summary.markdown %}).