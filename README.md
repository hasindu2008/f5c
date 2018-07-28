# f5c

**Nanopolish call-methylation in C**


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
Usage: ./f5c [OPTIONS] -r reads.fa -b alignments.bam -g genome.fa
```

