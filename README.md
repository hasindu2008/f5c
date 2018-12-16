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

## Building
Users are recommended to build from the [latest release](https://github.com/hasindu2008/f5c/releases/latest) tar ball. You need a compiler that supports C++11. Quick example for Ubuntu :
```
sudo apt-get install libhdf5-dev zlib1g-dev   # install HDF5 and zlib development library
wget "https://github.com/hasindu2008/f5c/releases/download/v0.0-alpha/f5c-v0.0-alpha-release.tar.gz"
tar xvf f5c-v0.0-alpha-release.tar.gz
cd f5c-v0.0-alpha/
scripts/install-hts.sh  # download and compile the htslib
./configure             
make                    # make cuda=1 to enable CUDA support
```
The commands to install HDF5 (and zlib development libraries) on some popular distributions :
```
On Debian/Ubuntu : sudo apt-get install libhdf5-dev zlib1g-dev
On Fedora/CentOS : sudo dnf/yum install hdf5-devel zlib-devel
On Arch Linux: sudo pacman -S hdf5
On OS X : brew install hdf5
```
If you cannot install HDF5 system wide, you can build it locally by skipping `scripts/install-hts.sh` and `./configure`. However, building HDF5 takes ages.
Building from the Github repository additionally requires `autoreconf` which can be installed on ubuntu using `sudo apt-get install autoconf automake`.
Other building options are detailed [here](building.md).
### NVIDIA CUDA support
To build for the GPU, you need to have the CUDA toolkit properly installed. Make sure you have added the nvcc (NVIDIA C Compiler) to your PATH.
The building instructions are the same as above except that you should call make as :
```
make cuda=1
```
Optionally you can provide the CUDA architecture as :
```
make cuda=1 CUDA_ARCH=-arch=sm_xy
```
If you get an error that /usr/local/cuda/ does not exist, for the moment you will have to create a symbolic link to direct to the correct CUDA installation or else edit the line `CUDALIB = -L/usr/local/cuda/lib64/ -lcudart_static -lrt -ldl` in the Makefile. We will make our installation scripts more intelligent in the future releases.

## Example usage
Follow the same steps as in [Nanopolish tutorial](https://nanopolish.readthedocs.io/en/latest/quickstart_call_methylation.html) while replacing `nanopolish` with `f5c`. If you only want to perform a quick test of f5c without aligning reads :
```
#download and extract the dataset including sorted alignments
wget -O f5c_na12878_test.tgz "http://genome.cse.unsw.edu.au/tmp/f5c_na12878_test.tgz"
tar xf f5c_na12878_test.tgz
#index and call methylation
f5c index -d chr22_meth_example/fast5_files chr22_meth_example//reads.fastq
f5c call-methylation -b chr22_meth_example/reads.sorted.bam -g chr22_meth_example/humangenome.fa -r chr22_meth_example/reads.fastq > chr22_meth_example/result.tsv
```

## Acknowledgement
This extensively reuses code and methods from [Nanopolish](https://github.com/jts/nanopolish).
The event detection code is from Oxford Nanopore's [Scrappie basecaller](https://github.com/nanoporetech/scrappie).
Some code snippets have been taken from [Minimap2](https://github.com/lh3/minimap2) and [Samtools](http://samtools.sourceforge.net/).
