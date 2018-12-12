---
title: Quick Start
---
If you are a Linux user and want to quickly try out download the compiled
binaries from the [latest release](https://github.com/swaywm/sway/releases/latest).
Note that currently the binary only runs on *x86_64* platform.

For example:

```sh
wget "https://github.com/hasindu2008/f5c/releases/download/v0.0-alpha/f5c-v0.0-alpha-binaries.tar.gz"
tar -xvf f5c-v0.0-alpha-binaries.tar.gz
cd f5c-v0.0-alpha/
./f5c_x86_64_linux        # CPU version
./f5c_x86_64_linux_cuda   # cuda supported version
```

Binaries should work on most Linux distributions and the only dependency is
`zlib` which is available by default on most distros.
