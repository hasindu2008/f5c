# f5c

**An attempt to re-implement Nanopolish call-methylation in C**

- ![#f03c15](https://placehold.it/15/f03c15/000000?text=+) `This project is under active development and is not complete.`

[![Build Status](https://travis-ci.com/hasindug/f5c.svg?token=pN7xnsxgLrRxbAn8WLVQ&branch=master)(https://travis-ci.com/hasindug/f5c)

Install the HDF5 (and zlib development libraries)
``` 
On Debian/Ubuntu : sudo apt-get install libhdf5-dev zlib1g-dev 
On Fedora/CentOS : sudo dnf/yum install hdf5-devel zlib-devel
```
## Building

```
git clone https://github.com/hasindug/f5c
cd f5c
autoreconf
src/install-htslib.sh
./configure
make
```

## Running

```

```

## Credits
This contains code and methods from [Nanopolish](https://github.com/jts/nanopolish).
The event detection code is from Oxford Nanopore's [Scrappie basecaller](https://github.com/nanoporetech/scrappie).
Some code snippets have been taken from [Minimap2](https://github.com/lh3/minimap2) and [Samtools](http://samtools.sourceforge.net/).