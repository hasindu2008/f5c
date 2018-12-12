---
title: Compile from source
---
Instead of using the pre-compiled binary, you can build the binary from source.

f5c depends on the following libraries:
- zlib
- htslib
- hdf5

## Build from tar release
f5c provides a minimal tar ball for users to build a binary on their computer.
```sh
scripts/install-hts.sh  # download and compile the htslib
./configure
make
```

## Build from git source

#### Method 1 (recommended)

Dependencies : Install the HDF5 (and zlib development libraries)
```sh
On Debian/Ubuntu : sudo apt-get install libhdf5-dev zlib1g-dev
On Fedora/CentOS : sudo dnf/yum install hdf5-devel zlib-devel
On Arch Linux: sudo pacman -S hdf5
On OS X : brew install hdf5
```

Now build f5c
```sh
git clone https://github.com/hasindug/f5c
cd f5c
autoreconf #not required if a release
scripts/install-hts.sh
./configure
make
```

#### Method 2 (time consuming)

Dependencies : Install the HDF5 and zlib development libraries
```
On Debian/Ubuntu : sudo apt-get install zlib1g-dev
On Fedora/CentOS : sudo dnf/yum install zlib-devel
```

Now build f5c
```sh
git clone https://github.com/hasindug/f5c
cd f5c
autoreconf #optional
scripts/install-hts.sh #optional
scripts/install-hdf5.sh #optional
./configure --enable-localhdf5 #optional
make
```

#### Method 3 (not recommended)

Dependencies : Install HDF5 and hts
```sh
On Debian/Ubuntu : sudo apt-get install libhdf5-dev zlib1g-dev libhts1
```

Now build f5c
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
