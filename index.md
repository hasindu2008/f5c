---
layout: article
---

**f5c** is an optimised re-implementation of the call-methylation module in
[Nanopolish](https://github.com/jts/nanopolish). Given a set of basecalled
Nanopore reads and the raw signals, f5c detects the methylated cytosine bases.
f5c can optionally utilise NVIDIA graphics cards for acceleration.

First the reads have to be indexed using `f5c index` (or `nanopolish index` -
`f5c index` is the same code as `nanopolish index`). Then invoke
`f5c call-methylation` to detect methylated cytosine bases. The result is
almost the same as from nanopolish except a few differences due to floating
point approximations.

[Latest release](https://github.com/hasindu2008/f5c/releases/latest)
