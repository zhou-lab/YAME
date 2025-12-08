<p align="center">
  <a href="">
    <img alt="Logo" src="https://github.com/user-attachments/assets/9384208f-deb3-4b2e-a574-b397dbc83ca4" height="240" />
  </a>
</p>

# YAME — Yet Another Methylation Encoder

[![Install with Bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/recipes/yame/README.html)
[![Documentation](https://img.shields.io/badge/docs-online-blue.svg)](https://zhou-lab.github.io/YAME/)

A fast and lightweight toolkit for storing, manipulating, and analyzing large-scale DNA methylation data at the sequence level.

For detailed documentation, tutorials, and usage examples, visit the [YAME User Guide](https://zhou-lab.github.io/YAME/).

## Overview

YAME is designed for efficient sequence-level DNA methylation data management, capable of handling both bulk and single-cell DNA methylome workflows. It introduces a family of compact binary formats (**CX formats**) that represent methylation values, MU counts, categorical states, fraction data, masks, and genomic coordinates in a uniform compressed structure.

### 🌟 Key Features

- **Ultra-fast performance** with high compression for methylation matrices
- **Scalable** to hundreds of thousands of single cells
- **Versatile data support**: MU counts, binary methylation, chromatin states, fractions, differential calls, and CpG coordinate streams
- **Comprehensive toolkit**: packing, unpacking, downsampling, subsetting, row operations, enrichment testing, and summarization
- **Consistent internal API**: all data stored as `cdata_t` blocks inside BGZF frames
- **Integrates seamlessly** with bedtools, KYCGKB, and other methylation workflows

## Installation

Install YAME using conda from the bioconda channel:

```bash
conda install yame -c bioconda
```

## Citation

If you use YAME in your research, please cite:

Goldberg*, Fu*, Atkins, Moyer, Lee, Deng, Zhou† (2025). "KnowYourCG: Facilitating Base-level Sparse Methylome Interpretation." *Science Advances*. [https://doi.org/10.1126/sciadv.adw3027](https://www.science.org/doi/10.1126/sciadv.adw3027)

## Support

- **Documentation**: [https://zhou-lab.github.io/YAME/](https://zhou-lab.github.io/YAME/)
- **Issues**: Please report bugs and feature requests on the [GitHub Issues page](https://github.com/zhou-lab/YAME/issues)

## License

YAME is dual-licensed:

- **AGPL-3.0** for academic, educational, and non-profit research use
- **Commercial License** for commercial applications

### Academic & Non-Profit Use
YAME is free to use for academic research, educational purposes, and non-profit organizations under the [GNU Affero General Public License v3.0 (AGPL-3.0)](LICENSE).

### Commercial Use
If you wish to use YAME in commercial products or services, or if the AGPL-3.0 restrictions are not suitable for your use case, please contact us for a commercial license: [zhouw3@chop.edu]

---

Developed by the Zhou Lab
