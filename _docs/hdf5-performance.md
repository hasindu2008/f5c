---
title: HDF5 Performance Issue
author: Thomas Lam
---

Oxford Nanopore reads are stored in `.fast5` files, based on the [HDF5](https://www.hdfgroup.org/HDF5/)
file format - one fast5 file for one read. The HDF5 file format is complicated and hence we have to rely on a HDF5 library to read those files. At the
moment there is only one HDF5 implementation, which is the [official one](https://www.hdfgroup.org/downloads/hdf5)
that includes libraries and utilities for handling HDF5 files.

After implementing IO interleaving[^1] we observed that on certain systems (mostly >16 CPU cores with mechanical HD RAID)
the processing time is much faster than fast5 file access. Running f5c on the small chr22_meth_example dataset (~150 Mbases) on an HPC (72 Intel Xeon Gold 6154 threads and RAID 6 setup of 12 HD) gave us:

```
[meth_main] Data loading time: 197.460 sec
[meth_main]     - bam load time: 3.209 sec
[meth_main]     - fasta load time: 26.853 sec
[meth_main]     - fast5 load time: 167.311 sec
[meth_main]         - fast5 open time: 94.801 sec
[meth_main]         - fast5 read time: 65.457 sec
[meth_main] Data processing time: 60.996 sec
```
Data loading takes three times the processing time! It was worse on a ~60Gbase NA12878 dataset. Data loading took ~69 hours and processing took only ~3 hours.

```
[meth_main] Data loading time: 249389.512 sec = 69.27 hours
[meth_main]     - bam load time: 2099.134 sec
[meth_main]     - fasta load time: 21418.816 sec
[meth_main]     - fast5 load time: 225797.798 sec
[meth_main]         - fast5 open time: 147881.362 sec
[meth_main]         - fast5 read time: 69539.236 sec
[meth_main] Data processing time: 10523.977 sec = 2.92 hours
```

The breakdown of the data loading time above shows that majority of the time is for fast5 loading. BAM file access consumes comparatively very little time - the access pattern is sequential. Fasta file access (which includes random access to fasta files containing the reference and the reads - uses faidx) takes ~6 hours. Fast5 access takes a massive amount of time : ~62 hours.


These times were measured when loading data into the data batch struct (`db_t` defined in [`f5c.h`](https://github.com/hasindu2008/f5c/blob/master/src/f5c.h#L199-L247)) which contains bam records, fasta cache and fast5 file) by embedding a time counter in our code (commit [`1357c4`](https://github.com/hasindu2008/f5c/commit/1357c403b2f60580055a31c24d672c9016001d93)).

```c
typedef struct {
    bam1_t **bam_rec;
    char **fasta_cache;
    fast5_t **f5;
    // and other struct members...
} db_t;
```

The header only file `fast5lite` uses the HDF5 C library heavily. If we trace down the call tree further, we can see that `fast5_open()` calls
`H5Fopen()`:

```c
static inline hid_t fast5_open(char *filename) {
    hid_t hdf5file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    return hdf5file;
}
```

And in `fast5_read()`, there are multiple calls to the HDF5 library (code
stripped to show only HDF5 library functions usage):

```c
static inline int32_t fast5_read(hid_t hdf5_file, fast5_t* f5) {
    hid_t space;
    hsize_t nsample;
    herr_t status;

    hid_t dset = H5Dopen(hdf5_file, signal_path, H5P_DEFAULT);
    space = H5Dget_space(dset);
    status = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     f5->rawptr);

    H5Sclose(space);
    H5Dclose(dset);

    // get channel parameters
    const char* scaling_path = "/UniqueGlobalKey/channel_id";

    hid_t scaling_group = H5Gopen(hdf5_file, scaling_path, H5P_DEFAULT);
    H5Gclose(scaling_group);
    free(signal_path);

    return 0;
}
```

In the chr22_meth_example test set, there are 2.3GB of fast5 files, which takes the HDF5 library a total of 167 seconds to open and read. This means that the average load speed is 14.103 MB/s - probably a lot of random accesses to the disk.

What if the fast5 files could be sequentially accessed? To find out the file reading speed of the disk, we used `dd` to generate a big
random file and `cat` it to find out the read time:

```console
$ dd if=/dev/urandom of=test.txt iflag=fullblock bs=64M count=16
16+0 records in
16+0 records out
1073741824 bytes (1.1 GB, 1.0 GiB) copied, 5.93124 s, 181 MB/s
$ sync; echo 3 | tee /proc/sys/vm/drop_caches # cleaning disk cache
$ time cat test.txt > /dev/null

real    0m0.738s
user    0m0.001s
sys     0m0.210s
```

So the disk is able to sequentially read at a rate of 1.355 GB/s, while libhdf5 can only load the files at 14.103 MB/s, whoppingly slower than a normal sequential disk read by 98.4 times.


At this point of the development, the data processing time is almost
always faster than the data loading time, meaning that regardless of the number of the optimisation techniques we put in, the fast5 file I/O will always be the
bottleneck.

[^1]: brief explanation on how IO interleaving works in f5c can be found [here](https://github.com/hasindu2008/f5c/blob/master/src/meth_main.c#L12-L23)
