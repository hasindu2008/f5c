---
title: Portable binaries
chart: true
---

As we value ease of use and portability of f5c, we are providing pre-compiled binaries in our releases. At the moment we provide the binaries for Linux on x86_64 architecture . Those pre-compiled binaries are portable; i.e. they work on most of the Linux distributions out there today. If you download and extract a release binary tar ball, you would see two binaries : one for the CPU only and the other with additional CUDA support.

To produce the portable binaries we follow a hybrid linking strategy where some libraries are statically linked while the others are dynamically linked. More reasons to use static linking can be found [here](http://lh3.github.io/2014/07/12/about-static-linking). The f5c release binary is statically-linked with [HTSlib](https://github.com/samtools/htslib), [HDF5](https://www.hdfgroup.org/solutions/hdf5). Additionally the NVIDIA CUDA runtime library is statically linked in the cuda binary. We only link libc, zlib and librt dynamically since we do not have a hard dependency on a particular version and they are in standard paths (`/usr/lib/` or `/lib/`) anyway.

You can produce the same binary with `make cuda=1` (without `./configure` before running this command) assuming that you have the CUDA runtime library installed[^1]. This will download HTSlib and HDF5 source code to a local working directory. They are then compiled to produce archive libraries (.a) and be linked with the f5c binary.

The release binaries are compiled on an Ubuntu 14 operating system with GCC 4.8 and CUDA toolkit 6.5 installed. Older versions are used as the glibc is dynamically linked. The performance of a release binary is on par with the one that is fresh compiled - compiled using newer compilers. The graph below shows that there is no significant performance difference between the two. The Dell XPS15 laptop used for testing had Ubuntu 18 as the operating system with GCC 7 and CUDA Toolkit 9.5.

```chart
{
  "type": "bar",
  "data": {
  "labels": [
    "CPU",
    "GPU"
  ],
  "datasets": [
    {
    "labels" : "pre-compiled release-binary",  
    "data": [ 178.2, 97.6 ],
    "backgroundColor": "rgba(255, 99, 132, 0.2)",
    "borderColor": "rgba(255,99,132, 1)",
    "borderWidth": 1
    },
    {
    "labels" : "fresh compiled binary",
    "data": [164.6,90.8],
    "backgroundColor": "rgba(54, 162, 235, 0.2)",
    "borderColor": "rgba(255,99,132, 1)",
    "borderWidth": 1      
    }
  ],
  "options": {
    "scales": {
      "yAxes": [{
        "scaleLabel": {
          "display": true,
          "labelString": "time/s"
        },
        "ticks": {
          "beginAtZero": true
        }
      }]
    }
  }
}
}
```

[^1]: See <https://github.com/hasindu2008/f5c/blob/master/.travis.yml> for the commands used to compile the statically-linked binary.
