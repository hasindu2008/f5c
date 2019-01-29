---
title: Static Compilation
chart: true
---

As we value ease of use and portability of f5c, we are providing statically-linked binary in our releases. More reasons to use static linking can be found [here](http://lh3.github.io/2014/07/12/about-static-linking).

The f5c release binary is statically-linked with [HTSlib](https://github.com/samtools/htslib), [HDF5](https://www.hdfgroup.org/solutions/hdf5) and the NVIDIA CUDA runtime library. You can produce the same binary with `make cuda=1` (without `./configure` before running this command) assuming that you have the CUDA runtime library installed[^1]. This will download HTSlib and HDF5 source code to a local working directory. They are then compiled to produce archive libraries (.a) and be linked with the f5c binary. We only link libc, zlib and librt dynamically since we do not have a hard dependency on a particular version and they are in standard paths (`/usr/lib/` or `/lib/`) anyway.

The performance of the statically linked binary is on par with the one that is dynamically-linked. The graph below shows that there is no significant performance difference between a statically-linked binary and a dynamically-linked one.

```chart
{
  "type": "bar",
  "data": {
  "labels": [
    "CPU (dynamically-linked)",
    "CPU (statically-linked)",
    "GPU (dynamically-linked)",
    "GPU (statically-linked)"
  ],
  "datasets": [
    {
    "label": "Run time (lower is better)",
    "data": [
      116.375,
      112.269,
      106.688,
      113.487
    ],
    "backgroundColor": [
      "rgba(255, 99, 132, 0.2)",
      "rgba(54, 162, 235, 0.2)",
      "rgba(255, 99, 132, 0.2)",
      "rgba(54, 162, 235, 0.2)"
    ],
    "borderColor": [
      "rgba(255,99,132, 1)",
      "rgba(54, 162, 235, 1)",
      "rgba(255,99,132, 1)",
      "rgba(54, 162, 235, 1)"
    ],
    "borderWidth": 1
    }
  ]
  },
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
```

[^1]: See <https://github.com/hasindu2008/f5c/blob/master/.travis.yml> for the commands used to compile the statically-linked binary.
