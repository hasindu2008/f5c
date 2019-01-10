---
title: HDF5 Performance Issue
author: Thomas Lam
---

The Oxford Nanopore reads are stored in `.fast5` files, based on the [HDF5](https://www.hdfgroup.org/HDF5/)
file format. So in order to process them we obviously need to use a HDF5 library
which provides us a way to read the file and process them accordingly. At the
moment there is only one HDF5 implementation, which is the [official one](https://www.hdfgroup.org/downloads/hdf5)
that includes libraries and utilities for dealing with HDF5 files.

When implementing IO interleaving[^1] we discovered that on certain platforms
the processing time is much faster than the file loading time. This is
interesting since theoretically the disks should be able to transfer files at a
higher rate than that. So we went on investigating which part of the code is
causing a bottleneck when reading in fast5 files.

We are interested to see the time breakdown of IO operations when initializing
the data batch struct (`db_t` defined in [`f5c.h`](https://github.com/hasindu2008/f5c/blob/master/src/f5c.h#L199-L247)),
which contains bam records, fasta cache and fast5 file.

```c
typedef struct {
    bam1_t **bam_rec;
    char **fasta_cache;
    fast5_t **f5;
    // and other struct members...
} db_t;
```

So we embedded a time counter in our code (commit [`1357c4`](https://github.com/hasindu2008/f5c/commit/1357c403b2f60580055a31c24d672c9016001d93))
to display a time breakdown for loading different types of files, which should
hopefully show us some clue.

Running f5c with chr22_meth_example dataset on an HPC shows us some really
interesting results:

```
[meth_main] Data loading time: 197.460 sec
[meth_main]     - bam load time: 3.209 sec
[meth_main]     - fasta load time: 26.853 sec
[meth_main]     - fast5 load time: 167.311 sec
[meth_main]         - fast5 open time: 94.801 sec
[meth_main]         - fast5 read time: 65.457 sec
[meth_main] Data processing time: 60.996 sec
```

As you can see bam filles and fasta files used up only 15% of the data load
time. The majority of the time used is for loading fast5 files. Now we can
confirm that the culprit of slow file IO time is the fast5 files operations. The
header file only `fast5lite` library uses the HDF5 C library heavily. If we
trace down the call tree further, we can see that `fast5_open()` calls
`H5Fopen()` for us:

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

This shows pretty clearly that libhdf5 is a bottleneck here. But we still have
to show that the disk is capable of reading a file at a higher speed to prove
that the HDF5 library has a slow file IO operation procedure. In the chr22_meth_example
test set, there are 2.3GB of fast5 files, which takes the HDF5 library a total of
167 seconds to open and read. This means that the load speed is 14.103 MB/s. To
find out the file reading speed of the disk, we used `dd` to generate a big
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

So the disk is able to read at a rate of 1.355 GB/s, while libhdf5 can only
process the file at 14.103 MB/s, whoppingly slower than a normal disk read by
98.4 times. At this point of the development, the data processing time is almost
always faster than the data load time, meaning that regardless of the number of
the optimization technique we put in, the fast5 files IO will always be the
bottleneck. Unless we have a better (faster) HDF5 implementation, f5c can only
perform as good as the HDF5 library.

[^1]: brief explanation on how IO interleaving works in f5c can be found [here](https://github.com/hasindu2008/f5c/blob/master/src/meth_main.c#L12-L23)
