---
layout: article
---

An optimised re-implementation of the *index*, *call-methylation* and *eventalign* modules in [Nanopolish](https://github.com/jts/nanopolish). Given a set of basecalled Nanopore reads and the raw signals, *f5c call-methylation* detects the methylated cytosine and *f5c eventalign* aligns raw nanopore signals (events) to the reference k-mers. *f5c* can optionally utilise NVIDIA graphics cards for acceleration. For best performance and easy usability, it is recommended to use f5c on [BLOW5 format](https://www.nature.com/articles/s41587-021-01147-4). Use [slow5tools](https://github.com/hasindu2008/slow5tools) for FAST5->BLOW5 conversion and [blue-crab](https://github.com/Psy-Fer/blue-crab) for POD5->BLOW5 conversion.

First, the reads have to be indexed using `f5c index`. Then, invoke `f5c call-methylation` to detect methylated cytosine bases. Finally, you may use `f5c meth-freq` to obtain methylation frequencies. Alternatively, invoke `f5c eventalign` to perform event alignment. The results are almost the same as from nanopolish except a few differences due to floating point approximations.

- **f5c v1.2 onwards support nanopore R10.4.1 chemistry (must specify --pore r10 if FAST5 input, autodetected for S/BLOW5 input)**.
- **f5c v1.4 onwards support nanopore RNA004 chemistry (make specify --pore rna004 if FAST5 input, autodetected for S/BLOW5 input)**.

*Full Documentation* : [https://hasindu2008.github.io/f5c/docs/overview](https://hasindu2008.github.io/f5c/docs/overview)
*Latest release* : [https://github.com/hasindu2008/f5c/releases/latest](https://github.com/hasindu2008/f5c/releases/latest)
*Pre-print* : [https://doi.org/10.1101/756122](https://www.biorxiv.org/content/10.1101/756122v1)
*Publication* : [https://doi.org/10.1186/s12859-020-03697-x](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03697-x)
*Supplementary*: [nanopore_signal_alignment_supplementary_material.pdf](https://hasindu2008.github.io/f5c/nanopore_signal_alignment_supplementary_material.pdf)

[![GitHub Downloads](https://img.shields.io/github/downloads/hasindu2008/f5c/total?logo=GitHub)](https://github.com/hasindu2008/f5c/releases)
[![BioConda Install](https://img.shields.io/conda/dn/bioconda/f5c?label=BioConda)](https://anaconda.org/bioconda/f5c)
[![Build Status](https://travis-ci.org/hasindu2008/f5c.svg?branch=master)](https://travis-ci.org/hasindu2008/f5c)

Please cite the following when using *f5c* in your publications:

> Gamaarachchi, H., Lam, C.W., Jayatilaka, G. et al. GPU accelerated adaptive banded event alignment for rapid comparative nanopore signal analysis. BMC Bioinformatics 21, 343 (2020). https://doi.org/10.1186/s12859-020-03697-x

```
@article{gamaarachchi2020gpu,
  title={GPU accelerated adaptive banded event alignment for rapid comparative nanopore signal analysis},
  author={Gamaarachchi, Hasindu and Lam, Chun Wai and Jayatilaka, Gihan and Samarakoon, Hiruna and Simpson, Jared T and Smith, Martin A and Parameswaran, Sri},
  journal={BMC bioinformatics},
  volume={21},
  number={1},
  pages={1--13},
  year={2020},
  publisher={BioMed Central}
}
