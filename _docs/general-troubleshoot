---
title: General Troubleshooting
---

### Error: Failed to allocate memory

If f5c exited with a "Failed to allocate memory" error, it means that the free RAM available on the system was not adequate. The peak RAM usage in *f5c* can be reduced by lowering the batch size parameters (-K and B). Reducing the number of threads (-t) also helps to reduce the peak RAM usage. 

Another factor that affects the RAM usage is the abundance of ultra-long reads. This is particularly the case if the error occurred in *align.c*. Peak RAM usage due to ultra-long reads can be reduced by skipping the ultra-long reads (--skip-ultra tmp.bam). If *f5c* still exits with the same error, you can reduce the threshold for ultra-long reads (e.g., --skip-ultra tmp.bam --ultra-thresh 50000).

### *f5c* crashes / get killed / segmentation fault

Linux allows memory overcommitment, i.e., programs are allowed to allocate more than the physical RAM available on the system, with the assumption that not all programmes will use the allocated memory at once.  However, when the system is running low on physical memory, the Out Of Memory (OOM) killer in Linux kills memory-hogging programmes. If *f5c* abruptly crashed (sometimes with the message "killed", but may also be a "segmentation fault"), it is likely to be due to high peak memory. You can verify this if you run the *dmesg* command on the terminal and observe a message like “Out of memory: Kill process X or sacrifice child”. Use the same techniques discussed above (Error: Failed to allocate memory topic) to reduce the peak memory usage of *f5c*.
