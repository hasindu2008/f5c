---
title: Compile from source
---
Instead of using the pre-compiled binary, you can build the binary from source.

While we have tried hard to avoid the dependency hell, three dependencies
(zlib, HDF5 and htslib) could not be avoided:
- zlib (`zlib1g-dev` on Debian/Ubuntu, `zlib-devel` on Fedora/CentOS, `zlib` on Arch Linux)
- htslib (`libhts1` on Debian/Ubuntu, and install `htslib` from AUR if you are on Arch Linux)
- hdf5 (`libhdf5-dev` on Debian/Ubuntu, `hdf5-devel` on Fedora/CentOS, `hdf5` on Arch Linux)

## Build from tar release
f5c provides a minimal tar ball for users to build a binary on their computer.

Download the [latest release](https://github.com/hasindu2008/f5c/releases/latest)
from GitHub, then extract and `cd` to the working directory.
```sh
scripts/install-hts.sh  # download and compile the htslib
./configure
make
```

## Build from git source

#### Method 1 (recommended)

Install HDF5 (and zlib development libraries) listed at the [top](http://127.0.0.1:4000/docs/compile-from-source).
This compiles only htslib locally and uses the system-wide HDF5 installation.

Now build f5c
```sh
git clone https://github.com/hasindug/f5c
cd f5c
autoreconf
scripts/install-hts.sh
./configure
make
```

#### Method 2 (time consuming)

This compiles all the libraries locally and statically link them with the
binary, which could take up to half an hour. Make sure you have zlib installed
on your computer.

Now build f5c
```sh
git clone https://github.com/hasindug/f5c
cd f5c
make
```

#### Method 3 (not recommended)

This method is not recommended as the pre-compiled htslib on Debian/Ubuntu is
not up to date and f5c depends on a newer version of htslib. Building on Arch
Linux should be fine though.

First install HDF5 and htslib using your package manager.

Then build f5c
```sh
git clone https://github.com/hasindug/f5c
cd f5c
autoreconf
./configure --enable-systemhts
make
```

## NVIDIA CUDA Support
To build for the GPU, you need to have the CUDA toolkit properly installed. Make sure you have added the nvcc (NVIDIA C Compiler) to your PATH.  

The building instructions are the same as above except that you should call make as:
```sh
make cuda=1
```
Optionally you can provide the CUDA architecture as:
```sh
make cuda=1 CUDA_ARCH=-arch=sm_xy
```

If you get an error that `/usr/local/cuda/` does not exist, for the moment you
will have to create a symbolic link to direct to the correct CUDA installation
or else edit the line `CUDALIB=-L/usr/local/cuda/lib64/ -lcudart_static -lrt -ldl`
in the Makefile. We will make our installation scripts more intelligent in the future releases.
