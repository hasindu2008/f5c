# f5c

An optimised re-implementation of the call-methylation module in [Nanopolish](https://github.com/jts/nanopolish). Given a set of basecalled Nanopore reads and the raw signals, f5c detects the methylated cytosine bases. f5c can optionally utilise NVIDIA graphics cards for acceleration.

First the reads have to be indexed using `f5c index` (or `nanopolish index` - f5c index is the same code as nanopolish index). Then invoke `f5c call-methylation` to detect methylated cytosine bases. The result is almost the same as from nanopolish except a few differences due to floating point approximations.

![#f03c15](https://placehold.it/15/f03c15/000000?text=+) `This project is under active development.`

[![Build Status](https://travis-ci.org/hasindu2008/f5c.svg?branch=master)](https://travis-ci.org/hasindu2008/f5c)

## Quick start

If you are a Linux user and want to quickly try out download the compiled binaries from the [latest release](https://github.com/hasindu2008/f5c/releases). For example:
```sh
wget "https://github.com/hasindu2008/f5c/releases/download/v0.0-alpha/f5c-v0.0-alpha-binaries.tar.gz"
tar xvf f5c-v0.0-alpha-binaries.tar.gz
cd f5c-v0.0-alpha/
./f5c_x86_64_linux        # CPU version
./f5c_x86_64_linux_cuda   # cuda supported version
```
Binaries should work on most Linux distributions and the only dependency is `zlib` which is available by default on most distros.

To build from source, check out the [documentation](https://hasindu2008.github.io/f5c/).

## Acknowledgement
This extensively reuses code and methods from [Nanopolish](https://github.com/jts/nanopolish).
The event detection code is from Oxford Nanopore's [Scrappie basecaller](https://github.com/nanoporetech/scrappie).
Some code snippets have been taken from [Minimap2](https://github.com/lh3/minimap2) and [Samtools](http://samtools.sourceforge.net/).
