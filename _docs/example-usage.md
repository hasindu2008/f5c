---
title: Example Usage
---
Follow the same steps in [Nanopolish tutorial](https://nanopolish.readthedocs.io/en/latest/quickstart_call_methylation.html)
while replacing `nanopolish` with `f5c`. If you only want to perform a quick test
of `f5c` without aligning reads, you can download a test dataset from
[here](http://genome.cse.unsw.edu.au/tmp/f5c_na12878_test.tgz).
```sh
# Extract dataset
tar xf f5c_na12878_test.tgz
# indexing and call methylation
f5c index -d chr22_meth_example/fast5_files chr22_meth_example//reads.fastq
f5c call-methylation -b chr22_meth_example/reads.sorted.bam -g chr22_meth_example/humangenome.fa -r chr22_meth_example/reads.fastq > chr22_meth_example/result.tsv
```
