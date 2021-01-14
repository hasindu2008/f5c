---
title: F5C performance hints
chart: true
---

Due to the vast differences among Nanopore datasets and the computing systems, having a universal set of parameters that gives the best runtime performance (execution time) in each case is impossible. If the default parameters are causing a significant impact on execution time, *f5c* will suggest user to tune the parameters. If preferred, users can terminate and relaunch *f5c* with better parameters. However, note that there may be false positives and false negatives.

## Table of Contents
- [Improving I/O performance](#io)
- [Improving CPU-GPU performance](#gpu)
  - [CPU-GPU Load balance](#loadbalance)
  - [GPU memory utilisation](#gpumem)

## <a name="io"></a> Improving I/O performance

#### Fast5 I/O taking more time than processing

If your fast5 files reside on a multi-disk system (a RAID array), you can launch with multiple I/O processes using *--iop* option. As a rule of thumb, if you have *n* disks, use *n* I/O processes.

If increased *iop* does not help, it means that your disk system cannot supply random disk requests fast enough to the processor. If your fast5s are located in a network mount storage, consider copying them to a fast local disk (ideally SSD).


## <a name="gpu"></a> Improving CPU-GPU performance

The tunable parameters in f5c that affects CPU-GPU performance are as follows.

| **parameter**  | **description**                                                                                           |
| :------------- | :-------------------------------------------------------------------------------------------------------- |
| *max-lf*       | reads with length <= *max-lf* x *average\_read\_length* are assigned to GPU and rest to CPU |
| *avg-epk*      | average number of events per base used for allocating GPU arrays       |
| *max-epk*      | reads with events per base <= *cuda-max-epk* are assigned to GPU, rest to CPU                        |
| *K*         | upper limit of the batch size with respect to the number of reads                                         |
| *B*          | upper limit of the batch size with respect to the number of bases                                         |
| *t*         | number of CPU threads                                                                                     |
| *ultra-thresh* | threshold to skip exceptionally long reads

For a comprehensive explanation of these parameters read our [pre-print](https://www.biorxiv.org/content/10.1101/756122v1.full).

### <a name="loadbalance"></a> CPU-GPU Load balance

When CPU and GPU are simultaneously processing parts of the batch of reads, either CPU or GPU taking too much time than the other, results in degraded overall performance. Hints to fix such load imbalance are as follows.

#### CPU time >> GPU time

If CPU is taking significantly more time than the GPU, the causes and the remedies can be :


- CPU processing too many reads due to a large batch size causing the GPU memory to overflow: decreasing -K may help
- CPU assigned with too many *very-long-reads* : increasing *--cuda-max-lf* may help
- CPU assigned with too many over-segmented reads : increasing *--cuda-max-epk* may help
- *Ultra-long-reads* which the CPU takes too much time for a given read: using *--skip-ultra* option may help. If still a problem decreasing *--ultra-thresh* would help.
- f5c was run with very few threads than the hardware CPU threads: increasing *-t* would help
- If nothing can fix, it may be that the CPU is too weak compared to the GPU and nothing can be done except upgrading the CPU


#### GPU time >> CPU time

If GPU is taking significantly more time than the CPU, the causes and the remedies can be : 

- GPU assigned with too many *very-long-reads* : decreasing *--cuda-max-lf* may help
- GPU assigned with too many over-segmented reads : decreasing *--cuda-max-epk* help
- CPU do not have much work to do probably due to a too low *Ultra-long-reads* threshold: increasing *--ultra-thresh* may help
- If nothing can fix, it may be that the CPU is too powerful compared to the GPU and nothing can be done except upgrading the GPU


### <a name="gpumem"></a> GPU memory utilisation

There are two types of arrays allocated on the GPU global memory, *read arrays* and *event arrays*. If either of them or both of them are under-utilised, the performance will not be ideal. Hints to maximise the GPU array utilisation is as follows.

#### all arrays under-utilised

If both *read arrays* and *event arrays* are under-utilised, it may be due to a very small batch size which is inadequate to fit the GPU memory. Try the following.

- increasing *-B*
- increasing *-K*

#### read arrays under-utilised

If only *read arrays* are under-utilised, this is mostly due to too many over-segmented reads assigned to the GPU or the default *average events per k-mer* value not adequate for this dataset. Try the following.

- decreasing *--max-epk*
- increasing *--avg-epk*


#### event arrays under-utilised


This is the inverse of the above and so try the reverse.

- increasing *--max-epk* 
- decreasing *--avg-epk*


