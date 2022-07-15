---
title: Quick Start
---
If you are a Linux user and want to quickly try out download the compiled
binaries from the [latest release](https://github.com/hasindu2008/f5c/releases/latest).
Note that currently the binary only runs on *x86_64* platform.

For example:

```sh
VERSION=v1.0
wget "https://github.com/hasindu2008/f5c/releases/download/$VERSION/f5c-$VERSION-binaries.tar.gz"
tar xvf f5c-$VERSION-binaries.tar.gz
cd f5c-$VERSION/
./f5c_x86_64_linux        # CPU version
./f5c_x86_64_linux_cuda   # cuda supported version
```

Binaries should work on most Linux distributions and the only dependency is
`zlib` which is available by default on most distros.
