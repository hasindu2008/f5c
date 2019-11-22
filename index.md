---
layout: article
---

**f5c** is an optimised re-implementation of the *call-methylation* and *eventalign* modules in [Nanopolish](https://github.com/jts/nanopolish). Given a set of basecalled Nanopore reads and the raw signals, *f5c call-methylation* detects the methylated cytosine and *f5c eventalign* aligns raw nanopore DNA signals (events) to the base-called read. *f5c* can optionally utilise NVIDIA graphics cards for acceleration.

First the reads have to be indexed using `f5c index`. Then invoke `f5c call-methylation` to detect methylated cytosine bases. Finally, you may use `f5c meth-freq` to obtain methylation frequencies. Alternatively, invoke `f5c eventalign` to perform event alignment. The results are almost the same as from nanopolish except a few differences due to floating point approximations.


[Go to Full Documentation](https://hasindu2008.github.io/f5c/docs/overview)

[Go to Latest release](https://github.com/hasindu2008/f5c/releases/latest)

[Go to pre-print](https://www.biorxiv.org/content/10.1101/756122v1)


