# f5c

**An attempt to re-implement Nanopolish call-methylation in C**

![#f03c15](https://placehold.it/15/f03c15/000000?text=+) `This project is under active development and is not complete.`

[![Build Status](https://travis-ci.org/hasindu2008/f5c.svg?branch=master)](https://travis-ci.org/hasindu2008/f5c)

## Building

Currently 3 methods are supported.
1. Locally compiled HTS library and system wide HDF5 library (recommended if you have root access)
2. Locally compiled HTS and HDF5 libraries (If you do not have root access. However compiling HDF5 takes a bit of time)
3. System wide HTS and HDF5 libraries (only works on Ubuntu 16 or higher)

### Method 1 (recommended if you have root access)

Dependencies : Install the HDF5 (and zlib development libraries)
``` 
On Debian/Ubuntu : sudo apt-get install libhdf5-dev zlib1g-dev 
On Fedora/CentOS : sudo dnf/yum install hdf5-devel zlib-devel
```

Now build f5c

```
git clone https://github.com/hasindug/f5c
cd f5c
autoreconf
scripts/install-htslib.sh
./configure
make
```

### Method 2 (if you have no root acesss)


```
git clone https://github.com/hasindug/f5c
cd f5c
autoreconf
scripts/install-htslib.sh
scripts/install-hdf5.sh
./configure --enable-localhdf5
make
```

### Method 3 

Dependencies : Install HDF5 and hts
``` 
On Debian/Ubuntu : sudo apt-get install libhdf5-dev zlib1g-dev libhts1
```

```
git clone https://github.com/hasindug/f5c
cd f5c
autoreconf
./configure --enable-systemhts
make
```

## Running

```
Usage: ./f5c [OPTIONS] -r reads.fa -b alignments.bam -g genome.fa
```

## Credits
This contains code and methods from [Nanopolish](https://github.com/jts/nanopolish).
The event detection code is from Oxford Nanopore's [Scrappie basecaller](https://github.com/nanoporetech/scrappie).
Some code snippets have been taken from [Minimap2](https://github.com/lh3/minimap2) and [Samtools](http://samtools.sourceforge.net/).
