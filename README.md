# f5c

An optimised re-implementation of the *index*, *call-methylation* and *eventalign* modules in [Nanopolish](https://github.com/jts/nanopolish). Given a set of basecalled Nanopore reads and the raw signals, *f5c call-methylation* detects the methylated cytosine and *f5c eventalign* aligns raw nanopore signals (events) to the reference k-mers. *f5c* can optionally utilise NVIDIA graphics cards for acceleration. For best performance and easy usability, it is recommended to use f5c on [BLOW5 format](https://www.nature.com/articles/s41587-021-01147-4). Use [slow5tools](https://github.com/hasindu2008/slow5tools) for FAST5->BLOW5 conversion and [blue-crab](https://github.com/Psy-Fer/blue-crab) for POD5->BLOW5 conversion.

First, the reads have to be indexed using `f5c index`. Then, invoke `f5c call-methylation` to detect methylated cytosine bases. Finally, you may use `f5c meth-freq` to obtain methylation frequencies. Alternatively, invoke `f5c eventalign` to perform event alignment. The results are almost the same as from nanopolish except a few differences due to floating point approximations.

- **f5c v1.2 onwards support nanopore R10.4.1 chemistry (must specify --pore r10 if FAST5 input, autodetected for S/BLOW5 input)**.
- **f5c v1.4 onwards support nanopore RNA004 chemistry (make specify --pore rna004 if FAST5 input, autodetected for S/BLOW5 input)**.

*Full Documentation* : [https://hasindu2008.github.io/f5c/docs/overview](https://hasindu2008.github.io/f5c/docs/overview)<br/>
*Latest release* : [https://github.com/hasindu2008/f5c/releases/latest](https://github.com/hasindu2008/f5c/releases/latest)<br/>
*Pre-print* : [https://doi.org/10.1101/756122](https://www.biorxiv.org/content/10.1101/756122v1)<br/>
*Publication* : [https://doi.org/10.1186/s12859-020-03697-x](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03697-x)<br/>
*Supplementary*: [nanopore_signal_alignment_supplementary_material.pdf](https://hasindu2008.github.io/f5c/nanopore_signal_alignment_supplementary_material.pdf)<br/>
*Talk Video* : [https://www.nvidia.com/en-us/on-demand/session/gtcspring21-s31391](https://www.nvidia.com/en-us/on-demand/session/gtcspring21-s31391)<br/>

[![GitHub Downloads](https://img.shields.io/github/downloads/hasindu2008/f5c/total?logo=GitHub)](https://github.com/hasindu2008/f5c/releases)
[![BioConda Install](https://img.shields.io/conda/dn/bioconda/f5c?label=BioConda)](https://anaconda.org/bioconda/f5c)
[![x86_64](https://github.com/hasindu2008/f5c/actions/workflows/f5c-x86_64.yml/badge.svg)](https://github.com/hasindu2008/f5c/actions/workflows/f5c-x86_64.yml)
[![AARCH64](https://www.travis-ci.com/hasindu2008/f5c.svg?branch=master)](https://app.travis-ci.com/hasindu2008/f5c)
[![CodeFactor](https://www.codefactor.io/repository/github/hasindu2008/f5c/badge/master)](https://www.codefactor.io/repository/github/hasindu2008/f5c/overview/master)

Please cite the following when using *f5c* in your publications:

> Gamaarachchi, H., Lam, C.W., Jayatilaka, G. et al. GPU accelerated adaptive banded event alignment for rapid comparative nanopore signal analysis. BMC Bioinformatics 21, 343 (2020). https://doi.org/10.1186/s12859-020-03697-x

```
@article{gamaarachchi2020gpu,
  title={GPU accelerated adaptive banded event alignment for rapid comparative nanopore signal analysis},
  author={Gamaarachchi, Hasindu and Lam, Chun Wai and Jayatilaka, Gihan and Samarakoon, Hiruna and Simpson, Jared T and Smith, Martin A and Parameswaran, Sri},
  journal={BMC bioinformatics},
  volume={21},
  number={1},
  pages={1--13},
  year={2020},
  publisher={BioMed Central}
}
```

## Quick start

If you are a Linux user and want to quickly try out download the compiled binaries from the [latest release](https://github.com/hasindu2008/f5c/releases). For example:
```sh
VERSION=v1.6
wget "https://github.com/hasindu2008/f5c/releases/download/$VERSION/f5c-$VERSION-binaries.tar.gz" && tar xvf f5c-$VERSION-binaries.tar.gz && cd f5c-$VERSION/
./f5c_x86_64_linux        # CPU version
./f5c_x86_64_linux_cuda   # cuda supported version
```
Binaries should work on most Linux distributions as the only dependency is `zlib` which is available by default on most distributions. For compiled binaries to work, your processor must support SSSE3 instructions or higher (processors after 2007 have these) and your operating system must have GLIBC 2.17 or higher (Linux distributions from 2014 onwards typically have this).

You can also use conda to install *f5c* as `conda install f5c -c bioconda -c conda-forge`.

From f5c v1.6 onwards, experimental binaries for AMD GPUs are also provided under [releases](https://github.com/hasindu2008/f5c/releases).

## Building from source

### Building a release

Users are recommended to build from the  [latest release](https://github.com/hasindu2008/f5c/releases) tar ball. You need a compiler that supports C++11. Quick example for Ubuntu :
```sh
sudo apt-get install libhdf5-dev zlib1g-dev   #install HDF5 and zlib development libraries
VERSION=v1.6
wget "https://github.com/hasindu2008/f5c/releases/download/$VERSION/f5c-$VERSION-release.tar.gz" && tar xvf f5c-$VERSION-release.tar.gz && cd f5c-$VERSION/
scripts/install-hts.sh  # download and compile the htslib
./configure
make                    # make cuda=1 to enable CUDA support
```
The commands to install hdf5 (and zlib) __development libraries__ on some popular distributions :
```sh
On Debian/Ubuntu : sudo apt-get install libhdf5-dev zlib1g-dev
On Fedora/CentOS : sudo dnf/yum install hdf5-devel zlib-devel
On Arch Linux: sudo pacman -S hdf5
On OS X : brew install hdf5
```

### Other building options

- Building from the Github repository additionally requires invoking `autoreconf --install` to generate the *configure* script. `autoreconf` can be installed on Ubuntu using `sudo apt-get install autoconf automake`.
- If you want only S/BLOW5 support you can disable HDF5 by invoking `./configure --disable-hdf5 && make`.
- If you skip `scripts/install-hts.sh` and `./configure`, hdf5 will be compiled locally. It is a good option if you cannot install hdf5 library system wide. However, building hdf5 takes ages.
- *f5c* from version 0.8.0 onwards by default requires vector instructions (SSSE3 or higher for Intel/AMD and neon for ARM) for builtin *slow5lib*. If your processor is an ancient processor with no such vector instructions, invoke make as `make no_simd=1`.
- You can optionally enable *zstd* support for builtin *slow5lib* when building *f5c* by invoking `make zstd=1`. This requires __zstd 1.3 development libraries__ installed on your system (*libzstd1-dev* package for *apt*, *libzstd-devel* for *yum/dnf* and *zstd* for *homebrew*).
- On Mac M1 or in any system if `./configure` cannot find the hdf5 libraries installed through the package manager, you can specify the location as *LDFLAGS=-L/path/to/shared/lib/ CPPFLAGS=-I/path/to/headers/*. For example on Mac M1:
	```
	./configure LDFLAGS=-L/opt/homebrew/lib/ CPPFLAGS=-I/opt/homebrew/include/
	make
	```
- Instructions to build a docker image and conda installation are detailed [here](https://hasindu2008.github.io/f5c/docs/misc-install).
- Other uncommon building options are detailed [here](https://hasindu2008.github.io/f5c/docs/building).
- An SIMD accelerated version contributed by [@dkhyland](https://github.com/dkhyland) is available in the [*simd* branch](https://github.com/hasindu2008/f5c/tree/simd).

### NVIDIA CUDA support

To build for the GPU, you need to have the CUDA toolkit installed. Make sure nvcc (NVIDIA C Compiler) is in your PATH.

The building instructions are the same as above except that you should call make as :
```
make cuda=1
```
Optionally you can provide the CUDA architecture as :
```
make cuda=1 CUDA_ARCH=-arch=sm_xy
```
If your CUDA library is not in the default location /usr/local/cuda/lib64, point to the correct location as:
```
make cuda=1 CUDA_LIB=/path/to/cuda/library/
```
Visit [here](https://hasindu2008.github.io/f5c/docs/cuda-troubleshoot) for troubleshooting CUDA related problems.

### AMD ROCM support

From f5c v1.6, AMD GPUs are supported. To build for such GPUs, you need to have the ROCM toolkit installed.

The building instructions are the same as above for the CPU, except that you should call make as :
```
make rocm=1
```
Optionally you can provide the ROCM architecture as :
```
make rocm=1 ROCM_ARCH=--offload-arch=gfxnnn
```
If your ROCM library is not in the default location /opt/rocm, point to the correct location as:
```
make rocm=1 ROCM_LIB=/path/to/rocm/library/
```

## Usage

```sh
# indexing #

f5c index -d [fast5_folder] [read.fastq|fasta] 		# for FAST5
f5c index --slow5 [slow5_file] [read.fastq|fasta]	# for S/BLOW5

# methylation calling #

# for S/BLOW5 (R9.4 or R10.4 DNA data)
f5c call-methylation -b [reads.sorted.bam] -g [ref.fa] -r [reads.fastq|fasta] --slow5 [slow5_file] > [meth.tsv]
# for FAST5, R9.4 DNA data
f5c call-methylation -b [reads.sorted.bam] -g [ref.fa] -r [reads.fastq|fasta] > [meth.tsv]
# for FAST5, R10.4 DNA data
f5c call-methylation -b [reads.sorted.bam] -g [ref.fa] -r [reads.fastq|fasta] --pore r10 > [meth.tsv]
# methylation frequency
f5c meth-freq -i [meth.tsv] > [freq.tsv]

# event align #

# for S/BLOW5 (R9.4 DNA/RNA or R10.4 DNA data)
f5c eventalign -b [reads.sorted.bam] -g [ref.fa] -r [reads.fastq|fasta] --slow5 [slow5_file] > [events.tsv]
# for FAST5 (R9.4 DNA data)
f5c eventalign -b [reads.sorted.bam] -g [ref.fa] -r [reads.fastq|fasta]  > [events.tsv]
# for FAST5 (R9.4 RNA data)
f5c eventalign -b [reads.sorted.bam] -g [ref.fa] -r [reads.fastq|fasta] --rna > [events.tsv]
# for FAST5 (R10.4 DNA data)
f5c eventalign -b [reads.sorted.bam] -g [ref.fa] -r [reads.fastq|fasta] --pore r10 > [events.tsv]
# for FAST5 (RNA004 RNA data)
f5c eventalign -b [reads.sorted.bam] -g [ref.fa] -r [reads.fastq|fasta] --pore rna004 --rna > [events.tsv]
```

Visit the [man page](https://hasindu2008.github.io/f5c/docs/commands) for all the commands and options. See [here](https://hasindu2008.github.io/f5c/docs/output) for description of output formats.

### Example

Follow the same steps as in [Nanopolish tutorial](https://nanopolish.readthedocs.io/en/latest/quickstart_call_methylation.html) while replacing `nanopolish` with `f5c`. If you only want to perform a quick test of f5c :
```sh
#download and extract the dataset including sorted alignments
wget -O f5c_na12878_test.tgz "https://f5c.bioinf.science/f5c_na12878_test"
tar xf f5c_na12878_test.tgz

###### Using S/BLOW5 as input (recommended) ######

#index, call methylation and get methylation frequencies
f5c index --slow5 chr22_meth_example/reads.blow5 chr22_meth_example/reads.fastq
f5c call-methylation --slow5 chr22_meth_example/reads.blow5 -b chr22_meth_example/reads.sorted.bam -g chr22_meth_example/humangenome.fa -r chr22_meth_example/reads.fastq > chr22_meth_example/result.tsv
f5c meth-freq -i chr22_meth_example/result.tsv > chr22_meth_example/freq.tsv
#event alignment
f5c eventalign --slow5 chr22_meth_example/reads.blow5 -b chr22_meth_example/reads.sorted.bam -g chr22_meth_example/humangenome.fa -r chr22_meth_example/reads.fastq > chr22_meth_example/events.tsv


###### Using FAST5 as input ######

#index, call methylation and get methylation frequencies
f5c index -d chr22_meth_example/fast5_files chr22_meth_example/reads.fastq
f5c call-methylation -b chr22_meth_example/reads.sorted.bam -g chr22_meth_example/humangenome.fa -r chr22_meth_example/reads.fastq > chr22_meth_example/result.tsv
f5c meth-freq -i chr22_meth_example/result.tsv > chr22_meth_example/freq.tsv
#event alignment
f5c eventalign -b chr22_meth_example/reads.sorted.bam -g chr22_meth_example/humangenome.fa -r chr22_meth_example/reads.fastq > chr22_meth_example/events.tsv

#### NOTE: If you are using FAST5 format, make sure to specify --pore r10 for R10.4.1 data and --rna for RNA data. These are autodetected for S/BLOW5 in latest f5c.

```

More examples can be found [here](https://hasindu2008.github.io/f5c/docs/example-usage).

## Acknowledgement
This reuses code and methods from [Nanopolish](https://github.com/jts/nanopolish).
The event detection code is from Oxford Nanopore's [Scrappie basecaller](https://github.com/nanoporetech/scrappie).
Some code snippets have been taken from [Minimap2](https://github.com/lh3/minimap2) and [Samtools](http://samtools.sourceforge.net/).
