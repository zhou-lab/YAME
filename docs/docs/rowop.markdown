---
title: Row Subset/Operations
nav_order: 6
---

# Row operations of the methylation sites 

`yame rowop` provides a fast and convenient way to perform some frequently used row operations such as obtaining the mean and standard deviation statistics from `yame rowop -o mean` and `yame rowop -o std`. 

The most common use of `yame rowop` might be taking the sum of binary data to merge the pseudobulks. For example, we have a single cell methylation `.cg` file, and we want to merge cluster #1 to obatin the cluster #1 pseudobulks, we can combine this with [`yame subset`]({% link docs/subset.markdown %}) by 

```bash
yame subset -l cluster_1_id.txt single_cell.cg | yame rowop - -o binasum > single_cell_pseudobulk1.cg 
```

`yame rowop` also allows 
1. Binarization of the `.cg` file with the output of binary strings using `yame rowop -o binstring` 
2. Co-methylation patterns of neighboring CGs using `yame rowop -o  cometh`. 

For more help with `rowop`, run `yame rowop` in the terminal or check out the
[rowop help page]({% link docs/subcommands/YAME_rowop.markdown %}).

`yame chunk` and `yame chunkchar` breakdown `.cx` file into smaller and more manageable parts.

For more help with `chunk`, run `yame chunk` in the terminal or check out the
[chunk help page]({% link docs/subcommands/YAME_split.markdown %}).
