---
title: Example Usage
---

## Simple example 

Follow the same steps as in [Nanopolish tutorial](https://nanopolish.readthedocs.io/en/latest/quickstart_call_methylation.html) while replacing `nanopolish` with `f5c` and `scripts/calculate_methylation_frequency.py` with `f5c meth-freq` in the commands. If you only want to perform a quick test of f5c without aligning reads :
```sh
#download and extract the dataset including sorted alignments
wget -O f5c_na12878_test.tgz "https://f5c.page.link/f5c_na12878_test"
tar xf f5c_na12878_test.tgz

#index and call methylation
f5c index -d chr22_meth_example/fast5_files chr22_meth_example/reads.fastq
f5c call-methylation -b chr22_meth_example/reads.sorted.bam -g chr22_meth_example/humangenome.fa -r chr22_meth_example/reads.fastq > chr22_meth_example/result.tsv

#calculate methylation frequency
f5c meth-freq -i chr22_meth_example/result.tsv > chr22_meth_example/freq.tsv

#perform event alignment
f5c eventalign -b chr22_meth_example/reads.sorted.bam -g chr22_meth_example/humangenome.fa -r chr22_meth_example/reads.fastq > chr22_meth_example/events.tsv
```

## Workflow for HPC

todo
  

## Methylation calling Workflow for a distributed system

todo: meth-freq stuff
  
