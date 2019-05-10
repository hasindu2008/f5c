#!/bin/ash

LD_PRELOAD=/home/hasindu/installs/gperftools/.libs/libtcmalloc.so ./f5c -b test/chr22_meth_example/reads.sorted.bam -g test/chr22_meth_example/humangenome.fa -r test/chr22_meth_example/reads.fastq -t64 --secondary=yes --min-mapq=0 --print-scaling=yes -K4096 > result.txt
