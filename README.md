<p align="center">
  <a href="">
    <img alt="Logo" src="https://github.com/user-attachments/assets/9384208f-deb3-4b2e-a574-b397dbc83ca4" height="240" />
  </a>
</p>

# YAME — Yet Another Methylation Encoder

[![Install from zhou-lab](https://img.shields.io/badge/install-zhou--lab-brightgreen.svg)](https://anaconda.org/zhou-lab/yame)
[![Documentation](https://img.shields.io/badge/docs-online-blue.svg)](https://zhou-lab.github.io/YAME/)

A bit-packer for DNA methylation data — arrays and whole genomes, where analysis is bitwise, and so stays fast from 28K probes to 29M CpGs.

For detailed documentation, tutorials, and usage examples, visit the [YAME User Guide](https://zhou-lab.github.io/YAME/).

## Overview

YAME packs DNA methylation into bits. A family of compact binary formats (**CX formats**) holds methylation values, MU counts, categorical states, fractions, masks and coordinate streams — as little as one bit per CpG — inside a single record layout.

The payoff is that questions become bit operations. A data file stores no coordinates: row *i* means whatever row *i* of the reference means, so intersecting a methylome with a feature set is a bitwise AND rather than a genomic join. An Infinium manifest and a whole-genome CpG set are the same kind of object under that rule — a row space, differing only in length — which is why the same commands serve array and sequencing data.

### 🌟 Key Features

- **Bit-level packing**: 1 bit per CpG for binary calls, 2 for set/universe — a whole-genome hg38 track is 3.5 MB
- **Array and sequencing alike**: hg38, mm10, mm39, MSA, EPICv2, EPIC, HM450, HM27, MM285, Mammal40
- **Scalable** to hundreds of thousands of single cells
- **Versatile data support**: MU counts, binary methylation, chromatin states, fractions, differential calls, and CpG coordinate streams
- **Comprehensive toolkit**: packing, unpacking, downsampling, subsetting, row operations, enrichment testing, and summarization
- **Consistent internal API**: all data stored as `cdata_t` blocks inside BGZF frames
- **Integrates seamlessly** with bedtools, KYCGKB, and other methylation workflows

## Installation

```bash
conda install -c zhou-lab -c conda-forge yame
```

The `zhou-lab` channel is published by CI on every release tag, so it is
always current. The [bioconda recipe](https://bioconda.github.io/recipes/yame/README.html)
lags well behind and its build predates `yame fetch`, so it cannot download
the reference data the other commands resolve `-R` and `-m` against — prefer
the channel above until that catches up.

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
