# f5c

An optimised re-implementation of the call-methylation module in [Nanopolish](https://github.com/jts/nanopolish). Given a set of basecalled Nanopore reads and the raw signals, f5c detects the methylated cytosine bases. f5c can optionally utilise NVIDIA graphics cards for acceleration.

First the reads have to be indexed using `f5c index` (or `nanopolish index` - f5c index is the same code as nanopolish index). Then invoke `f5c call-methylation` to detect methylated cytosine bases. The result is almost the same as from nanopolish except a few differences due to floating point approximations.

![#f03c15](https://placehold.it/15/f03c15/000000?text=+) `This project is under active development.`

[![Build Status](https://travis-ci.org/hasindu2008/f5c.svg?branch=master)](https://travis-ci.org/hasindu2008/f5c)

## Quick start


## Building

Compiled binaries for Linux for x86-64 are available in the releases.


### Building for CPU only

You need a compiler that supports c++11.
While we have tried hard to avoid the dependency hell, three dependencies (zlib, HDF5 and HTS) could not be avoided.

Currently 3 building methods are supported.
1. Locally compiled HTS library and system wide HDF5 library (recommended)
2. Locally compiled HTS and HDF5 libraries (HDF5 local compilation - takes a bit of time)
3. System wide HTS and HDF5 libraries (not recommended as HTS versions can be old)

#### Method 1 (recommended)

Dependencies : Install the HDF5 (and zlib development libraries)
```
On Debian/Ubuntu : sudo apt-get install libhdf5-dev zlib1g-dev
On Fedora/CentOS : sudo dnf/yum install hdf5-devel zlib-devel
On Arch Linux: sudo pacman -S hdf5
On OS X : brew install hdf5
```

Now build f5c
```
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
```
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
```
On Debian/Ubuntu : sudo apt-get install libhdf5-dev zlib1g-dev libhts1
```

Now build f5c
```
git clone https://github.com/hasindug/f5c
cd f5c
autoreconf
./configure --enable-systemhts
make
```

### building with NVIDIA CUDA support

To build for the GPU, you need to have the CUDA toolkit installed. Make sure you have added the nvcc (NVIDIA C Compiler) to your PATH.
The building instructions are the same except that you should call make as :
```
make cuda=1
```
Optionally you can provide the CUDA architecture as :
```
make cuda=1 CUDA_ARCH=-arch=sm_xy
```



## Example

Follow the same steps as in [Nanopolish tutorial](https://nanopolish.readthedocs.io/en/latest/quickstart_call_methylation.html) while replacing `nanopolish` with `f5c`. If you only want to perform a quick test of f5c without aligning reads :
```
#download and extract the dataset including sorted alignments
wget -O f5c_na12878_test.tgz "http://genome.cse.unsw.edu.au/tmp/f5c_na12878_test.tgz"
tar xf f5c_na12878_test.tgz

#index and call methylation
f5c -d index chr22_meth_example/fast5_files chr22_meth_example//reads.fastq
f5c call-methylation -b chr22_meth_example/reads.sorted.bam -g chr22_meth_example/humangenome.fa -r chr22_meth_example/reads.fastq > chr22_meth_example/result.tsv
```




## Usage

### Indexing

```
Usage: f5c index [OPTIONS] -d nanopore_raw_file_directory reads.fastq
Build an index mapping from basecalled reads to the signals measured by the sequencer
f5c index is equivalent to nanopolish index by Jared Simpson

  -h, --help                           display this help and exit
  -v, --verbose                        display verbose output
  -d, --directory                      path to the directory containing the raw ONT signal files. This option can be given multiple times.
  -s, --sequencing-summary             the sequencing summary file from albacore, providing this option will make indexing much faster
  -f, --summary-fofn                   file containing the paths to the sequencing summary files (one per line)
```

### Calling methylation

```
Usage: f5c call-methylation [OPTIONS] -r reads.fa -b alignments.bam -g genome.fa
   -r FILE                 fastq/fasta read file
   -b FILE                 sorted bam file
   -g FILE                 reference genome
   -t INT                  number of threads [8]
   -K INT                  batch size (number of reads loaded at once) [512]
   -h                      help
   --min-mapq INT          minimum mapping quality [30]
   --secondary             consider secondary mappings or not [no]
   --skip-unreadable       skip any unreadable fast5 or terminate program [yes]
   --verbose INT           verbosity level [0]
   --version               print version
   --disable-cuda          disable running on CUDA [no] (only if compiled for CUDA)
debug options:
   --kmer-model FILE       custom k-mer model file (used for debugging)
   --print-events          prints the event table (used for debugging)
   --print-banded-aln      prints the event alignment (used for debugging)
   --print-scaling         prints the estimated scalings (used for debugging)
   --print-raw             prints the raw signal (used for debugging)
   --debug-break           break after processing the first batch (used for debugging)
```


## Docker image

To build a docker image
```
git clone https://github.com/hasindug/f5c
cd f5c
docker build .
```

Note down the image uuid and run f5c as :
```
docker run -v /path/to/local/data/data/:/data/ -it :image_id  ./f5c call-methylation -r /data/reads.fa -b /data/alignments.sorted.bam -g /data/ref.fa
```



## Acknowledgement
This extensively reuses code and methods from [Nanopolish](https://github.com/jts/nanopolish).
The event detection code is from Oxford Nanopore's [Scrappie basecaller](https://github.com/nanoporetech/scrappie).
Some code snippets have been taken from [Minimap2](https://github.com/lh3/minimap2) and [Samtools](http://samtools.sourceforge.net/).
