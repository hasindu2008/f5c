---
title: Parameter Profiles
---

For achieving maximal runtime performance in f5c, performance governing parameters must be selected by the user as explained in [f5c-perf-hints](https://hasindu2008.github.io/f5c/docs/f5c-perf-hints). From f5c-v0.5 onwards, preset parameter profiles are provided for a wide range of computer systems from embedded systems to servers. The user can simply specify the parameter profile to f5c with *-x* option followed by the name of the profile (See tables below) that closely matches the specification of the user's computer system.  Performance governing parameters currently set by a profile are *t*, *K*, *B*, *max-lf*, *avg-epk* and *max-epk*. Note that the user can also provide any of these parameters explicitly in addition to a profile. In such a circumstance, parameters set by *-x* will be overridden by the user provided parameters.

In the tables below, the number of cores includes hyper-threads in Intel processors. The GPU specification can be ignored if you are using the CPU only version of f5c.

## Generic profiles

| Profile name                       | Description                                                                  |
|------------------------------------|------------------------------------------------------------------------------|
| *laptop-low*                       | A low-end laptop with a 4-core CPU, 4 GB RAM and a GPU with 2 GB memory      |
| *laptop-mid* (or simply *laptop*)  | A mid-range laptop with an 8-core CPU, 8 GB RAM and a GPU with 3 GB memory    |
| *laptop-high*                      | A high-end laptop with a 12-core CPU, 16 GB RAM and a GPU with 4 GB memory   |
| *desktop-low*                      | A low-end workstation with an 8-core CPU, 32 GB RAM and a GPU with 8 GB memory      |
| *desktop-mid* (or simply *desktop*)| A mid-range workstation with a 12-core CPU, 32 GB RAM and a GPU with 10 GB memory    |
| *desktop-high*                     | A high-end workstation with a 16-core CPU, 64 GB RAM and a GPU with 12 GB memory   |
| *hpc-low*                      | A low-end server with a 32-core CPU, 128 GB RAM and a GPU with 16 GB memory      |
| *hpc-mid* (or simply *hpc*)| A mid-range server with a 48-core CPU, 256 GB RAM and a GPU with 32 GB memory    |
| *hpc-high*                     | A high-end server with a 64-core CPU, 384 GB RAM and a GPU with 40 GB memory   |

## Specific profiles

| Profile name                       | Description                                                                  |
|------------------------------------|------------------------------------------------------------------------------|
| *jetson-nano*                       | An NVIDIA Jetson Nano with a 4-core CPU and 4 GB integrated memory      |
| *jetson-tx2*                      | An NVIDIA Jetson TX2 with a 6-core CPU and 8 GB integrated memory    |
| *jetson-xavier*                      | An NVIDIA Jetson Xavier with an 8-core CPU and 16 GB integrated memory  |
| *nci-gadi*                      | A GPU node in the [NCI Gadi supercomputer](https://nci.org.au/our-systems/hpc-systems) equipped with Tesla V100      |



### Notes:
- The above profiles are for a typical nanopore dataset and if your dataset properties (e.g., read lengths, events per k-mer) are quite different, the user may still have to follow steps in [f5c-perf-hints](https://hasindu2008.github.io/f5c/docs/f5c-perf-hints) for best performance.
- The above list of computer systems is not exhaustive. I will keep adding new profiles whenever I get hold of a new computer system. If you wish a new profile to be added for a particular computer system, please open an [issue](https://github.com/hasindu2008/f5c/issues).
- Some of the above profiles are based on values obtained by running [f5c-tools](https://github.com/dkhyland/f5c-tools). As this method requires a number of iterations and thus very time consuming, some profiles contain interpolated values.
- The parameter values set by a particular profile can be found [here](). 

### Acknowledgement:
This parameters profile feature in f5c was developed by [David Hyland](https://github.com/dkhyland).
