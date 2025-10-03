---
title: Quick Start
---
If you are a Linux user and want to quickly try out download the compiled
binaries from the [latest release](https://github.com/hasindu2008/f5c/releases/latest).
Note that currently the binary only runs on *x86_64* platform.

For example:

```sh
VERSION=v1.6
wget "https://github.com/hasindu2008/f5c/releases/download/$VERSION/f5c-$VERSION-binaries.tar.gz"
tar xvf f5c-$VERSION-binaries.tar.gz
cd f5c-$VERSION/
./f5c_x86_64_linux        # CPU version
./f5c_x86_64_linux_cuda   # cuda supported version
```

Binaries should work on most Linux distributions as the only dependency is `zlib` which is available by default on most distributions. For compiled binaries to work, your processor must support SSSE3 instructions or higher (processors after 2007 have these) and your operating system must have GLIBC 2.17 or higher (Linux distributions from 2014 onwards typically have this).

You can also use conda to install *f5c* as `conda install f5c -c bioconda -c conda-forge`.

From f5c v1.6 onwards, experimental binaries for AMD GPUs are also provided under [releases](https://github.com/hasindu2008/f5c/releases).
