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

## Generic methylation calling workflow

In the following example, the system is assumed to have an 8-core CPU. You can change the number of threads depending on your CPU. Increase -K and -B for better multi-threaded performance at cost of higher peak RAM.

```sh
#align reads using minimap2 and sort using samtools
minimap2 -x map-ont -a -t8 --secondary=no ref.fa reads.fq > reads.sam
samtools sort -@8 reads.sam > reads.bam
samtools index reads.bam

#index fast5s, call methylation and count frequencies
f5c index --iop 8 -t 8 -d fast5/ reads.fq 
f5c call-methylation -t 8 -r reads.fq -g ref.fa -b reads.bam -K 512 -B 2M > meth.tsv
f5c meth-freq -i meth.tsv -s > meth-freq.tsv
```

## Resource efficient methylation calling workflow for a dataset with many ultra-long reads

In the following example, the system is assumed to have an 32-core CPU and a GPU with 16 GB memory. You can change the number of threads depending on your CPU. Set -K and -B depending on the available GPU memory. 

```sh
#align reads using minimap2 and sort using samtools
minimap2 -x map-ont -a -t32 --secondary=no ref.fa reads.fq > reads.sam
samtools sort -@32 reads.sam > reads.bam
samtools index reads.bam

#index fast5s
f5c index --iop 64 -t 32 -d fast5/ reads.fq

#call methylation while skipping the ultra-long reads longer than 100,000 and writing those entries to tmp.bam
f5c call-methylation -t 32 -r reads.fq -g ref.fa -b reads.bam -K 1024 -B 10M --skip-ultra tmp.bam --ultra-thresh 100000 > meth.tsv

#call methylation on the ultra-long reads separately (tmp.bam) on the CPU with a lower batch size to cap the peak RAM
f5c call-methylation -t 32 -r reads.fq -g ref.fa -b tmp.bam -K 512 -B 5M --disable-cuda=yes > meth-ultra.tsv

#call methylation frequencies and merge them
f5c meth-freq -i meth.tsv -s > meth-freq.tsv
f5c meth-freq -i meth-ultra.tsv -s > meth-ultra-freq.tsv
f5c freq-merge meth-freq.tsv meth-ultra-freq.tsv > meth-freq-combined.tsv
```
  

## Methylation calling workflow for a dataset containing independent batches

Assume the dataset is composed of several batches which is usually the case if real-time base-calling was performed. You can perform an individual workflow independently on each batch and finally combine the methylation frequency counts. 


```bash
#perform methylation calling on each batch (10 batches in this example)
for i in {1..10}
do
  minimap2 -x map-ont -a -t8 --secondary=no ref.fa reads_i.fq > reads_i.sam
  samtools sort -@8 reads_i.sam > reads_i.bam
  samtools index reads_i.bam
  f5c index --iop 8 -t 8 -d fast5_i/ reads_i.fq 
  f5c call-methylation -t 8 -r reads_i.fq -g ref.fa -b reads_i.bam -K 512 -B 2M > meth_i.tsv
done

#merge the frequency counts
f5c freq-merge reads_{1..10}.sam > meth-freq.tsv
```

Frequency count merging can also be useful to utilise a distributed system (e.g. array job in an HPC environment) or when performing real-time methylation calling.

