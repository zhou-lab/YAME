<p align="center">
  <a href="">
    <img alt="Logo" src="https://github.com/user-attachments/assets/9384208f-deb3-4b2e-a574-b397dbc83ca4" height="240" />
  </a>
</p>


# YAME — Yet Another Methylation Encoder

A methylation toolset designed for sequence-level DNA methylation data management. It is a command-line C program capable of performing sequence-level enrichment testing, row operations (such as merging pseudobulks), downsampling, and other related tasks with ultra fast speed.

**YAME** is a fast and lightweight toolkit for storing, manipulating, and analyzing large-scale DNA methylation data.  
It introduces a family of compact binary formats (**CX formats**) that represent methylation values, MU counts, categorical states, fraction data, masks, and genomic coordinates in a uniform compressed structure.

YAME provides command-line tools for:

- Efficiently **packing** text files into `.cx` formats  
- **Unpacking** and exporting `.cx` back to human-readable form  
- **Downsampling**, **subsetting**, and **row operations**  
- **Enrichment testing** and methylation **summarization**  
- Managing **multi-sample** CX files and sample indices  

YAME is designed for both bulk and single-cell DNA methylome workflows.

---

## 🌟 Key Features

- **High compression performance** on methylation matrices  
- Supports **MU counts**, **binary methylation**, **chromatin states**, **fractions**, **differential calls**, and **CpG coordinate streams**  
- Scales to **hundreds of thousands** of single cells  
- Consistent internal API: all data stored as `cdata_t` blocks inside BGZF frames  
- Integrates naturally with bedtools, KYCGKB, and other methylation workflows  
- Fully documented with structured guides in the `docs/` folder

---


## INSTALL

```
conda install yame -c bioconda
```

### Citing YAME

Goldberg*, Fu*, Atkins, Moyer, Lee, Deng, Zhou† ["KnowYourCG: Facilitating Base-level Sparse Methylome Interpretation"](https://www.science.org/doi/10.1126/sciadv.adw3027) Science Advances (2025)

## Usage and Documentation

A User Guide has been created to provide detailed documentation of YAME. The guide can be found at: https://zhou-lab.github.io/YAME/.
