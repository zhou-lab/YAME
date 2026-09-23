<div align="center">

# YAME — Yet Another Methylation Encoder

[![build](https://github.com/zhou-lab/YAME/actions/workflows/conda-build.yml/badge.svg)](https://github.com/zhou-lab/YAME/actions/workflows/conda-build.yml)
[![conda](https://img.shields.io/conda/vn/zhou-lab/yame?label=conda)](https://anaconda.org/zhou-lab/yame)
[![license](https://img.shields.io/badge/license-BSD--2--Clause%20(academic)%20%2F%20commercial-blue.svg)](LICENSE)
[![coverage](https://img.shields.io/endpoint?url=https%3A%2F%2Fzhou-lab.github.io%2FYAME%2Fcoverage.json)](https://github.com/zhou-lab/YAME/blob/main/test/run.sh)
[![docs](https://img.shields.io/badge/docs-online-blue.svg)](https://zhou-lab.github.io/YAME/)

**A bit-packer for DNA methylation data** — arrays and whole genomes, where
analysis is bitwise, and so stays fast from 28K probes to 29M CpGs.

📖 **[Documentation](https://zhou-lab.github.io/YAME/)** ·
🤖 [llms.txt](https://zhou-lab.github.io/YAME/llms.txt) ·
📦 [Install](https://anaconda.org/zhou-lab/yame)

</div>

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

Use of this software is available to academic and non-profit institutions for
research purposes under the [2-Clause BSD License](LICENSE). For use or
transfers of the software to commercial entities, please inquire with
Dr. Wanding Zhou at zhouw3@chop.edu. © 2021-present The Children's Hospital
of Philadelphia.

---

Developed by the Zhou Lab
