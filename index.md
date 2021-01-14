---
layout: article
---

**f5c** is an optimised re-implementation of the *call-methylation* and *eventalign* modules in [Nanopolish](https://github.com/jts/nanopolish). Given a set of basecalled Nanopore reads and the raw signals, *f5c call-methylation* detects the methylated cytosine and *f5c eventalign* aligns raw nanopore DNA signals (events) to the reference k-mers. *f5c* can optionally utilise NVIDIA graphics cards for acceleration.

First, the reads have to be indexed using `f5c index`. Then, invoke `f5c call-methylation` to detect methylated cytosine bases. Finally, you may use `f5c meth-freq` to obtain methylation frequencies. Alternatively, invoke `f5c eventalign` to perform event alignment. The results are almost the same as from nanopolish except a few differences due to floating point approximations.


*Full Documentation* : [https://hasindu2008.github.io/f5c/docs/overview](https://hasindu2008.github.io/f5c/docs/overview)

*Latest release* : [https://github.com/hasindu2008/f5c/releases/latest](https://github.com/hasindu2008/f5c/releases/latest)

*Pre-print* : [https://doi.org/10.1101/756122](https://www.biorxiv.org/content/10.1101/756122v1)

*Publication* : [https://doi.org/10.1186/s12859-020-03697-x](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03697-x)

[![GitHub Downloads](https://img.shields.io/github/downloads/hasindu2008/f5c/total?logo=GitHub)](https://github.com/hasindu2008/f5c/releases)
[![BioConda Install](https://img.shields.io/conda/dn/bioconda/f5c?label=BioConda)](https://anaconda.org/bioconda/f5c)
[![Build Status](https://travis-ci.org/hasindu2008/f5c.svg?branch=master)](https://travis-ci.org/hasindu2008/f5c)



