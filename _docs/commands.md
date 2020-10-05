---
title: Commands and Options
---

## NAME

f5c(1) - Ultra-fast methylation calling and event alignment tool for nanopore sequencing data with optional GPU (CUDA) acceleration

## SYNOPSIS

* indexing:
  ```
  f5c index -d [fast5_folder] [read.fastq|fasta]
  ```
* methylation calling:
  ```
  f5c call-methylation -b [reads.sorted.bam] -g [ref.fa] -r [reads.fastq|fasta] > [meth.tsv]
  f5c meth-freq -i [meth.tsv] > [freq.tsv]
  ```
* event alignment:
  ```
  f5c eventalign -b [reads.sorted.bam] -g [ref.fa] -r [reads.fastq|fasta] > [events.tsv]
  ```


## DESCRIPTION

Given a set of base-called nanopore reads and associated raw signals, f5c call-methylation detects the methylated cytosine at genomic CpG cites and f5c eventalign aligns raw nanopore DNA signals (events) to the base-called read. f5c can optionally utilise CUDA enabled NVIDIA graphics cards for acceleration. f5c is a heavily re-engineered and optimised implementation of the call-methylation and eventalign modules in Nanopolish.

## COMMANDS 

* `index`:               
         Build an index map from basecalled reads to the signals measured by the sequencer (same as nanopolish index).
* `call-methylation`:    
         Classify nucleotides as methylated or not at genomic CpG sites (optimised nanopolish call-methylation).
* `meth-freq`:           
         Calculate methylation frequency at genomic CpG sites (optimised nanopolish calculate_methylation_frequency.py).
* `eventalign`:          
         Align nanopore events to reference k-mers (optimised nanopolish eventalign).


## OPTIONS

### index

`f5c index [OPTIONS] -d nanopore_raw_file_directory reads.fastq`

Build an index map from basecalled reads to the signals measured by the sequencer. f5c index is equivalent to nanopolish index by Jared Simpson.

*  `-h`, `--help`:                           
         display this help and exit
*  `-v`, `--verbose`:                        
         display verbose output
*  `-d`, `--directory`:                      
         path to the directory containing the raw ONT signal files. This option can be given multiple times.
*  `-s`, `--sequencing-summary`:             
         the sequencing summary file from albacore, providing this option will make indexing much faster
*  `-f`, `--summary-fofn`:                   
         file containing the paths to the sequencing summary files (one per line)


### call-methylation

`f5c call-methylation [OPTIONS] -r reads.fa -b alignments.bam -g genome.fa`

Classify nucleotides as methylated or not at genomic CpG cites (optimised nanopolish call-methylation). Note that the list below contains the options for both CPU-only and CPU-GPU versions of f5c.  Options related to the GPU (CUDA) do NOT apply to the CPU-only version.

#### basic options

* `-r FILE`:                     
  The file containing the base-called reads in FASTQ or FASTA format. Can be gzip compressed files.
* `-b FILE`:                    
  The file contaning the alignment records sorted based on genomic coordinates in BAM format.
* `-g FILE`:                    
  The file containing the reference genome in FASTA format.
* `-w STR`:      
  Only process the specified genomic region STR. STR should be in the format *chr:start-end*. Currently, multiple region strings are not supported. If this option is not specified, the whole genome will be processed.
* `-t INT`:                     
  Number of processing threads [default value: 8].
* `-K INT`:                     
  Maximum number of reads loaded at once to the memory [default value: 512]. A larger value maximises multithreading performance at cost of increased peak RAM.
* `-B FLOAT[K/M/G]`:            
  Maximum number of bases loaded at once to the memory [default value: 2.0M]. A larger value maximises multithreading performance at cost of increased peak RAM.
* `-h`:                         
  Print the help to the standard out.
* `-o FILE`:                    
  The file to write the output. If this option is not specified, the output will be written to the standard out. 
* `-x STR`:                     
  Parameter profile to be used for maximising the performance to a particular computer system. The profile parameters are always applied before other options, i.e., the user can override these parameters explicitly. Some example profiles are laptop, desktop, hpc. See [profiles](https://f5c.page.link/profiles) for the full list and details.
* `--iop INT`:                  
  Number of I/O processes to read FAST5 files [default value: 1]. Increase this value if reading FAST5 limits the overall performance. A higher value is always preferred for systems with multiple disks (RAID) and network file systems.
* `--min-mapq INT`:             
  Minimum mapping quality of an alignment (MAPQ in the BAM record) to be considered for methylation calling [default value: 20].
* `--secondary=yes|no`:         
  Whether secondary alignments are considered or not for methylation calling [default value: no].
* `--verbose INT`:              
  Verbosity level for the log messages [default value: 0].
* `--version`:                  
  Print the version number to the standard out.
* `--disable-cuda=yes|no`:      
  Disable running on the GPU or not [default value: no]. If this option is set to yes, GPU acceleration is disabled. 
* `--cuda-dev-id INT`:          
  CUDA device identifier to run GPU kernels on [default value: 0]. The device identifier of the first GPU is 0, the second GPU is 1 and so on. This can be found by invoking the  `nvidia-smi` command. Currently, only a single GPU can be specified. To utilise multiple GPUs, you have to manually invoke multiple f5c commands on different datasets with a different device identifier.  
* `--cuda-max-lf FLOAT`:        
  Process reads with read-length less than or equal to the product of *cuda-max-lf* and the average read length in the current batch on GPU. The rest is processed on CPU [default value: 3.0]. Useful for tuning the CPU-GPU load balance for atypical datasets. Refer to [performance guidelines](https://hasindu2008.github.io/f5c/docs/f5c-perf-hints) for details.
* `--cuda-avg-epk FLOAT`:       
  The average number of events-per-kmer used for allocating the arrays in GPU memory [default value: 2.0]. Useful for tuning the CPU-GPU load balance for atypical datasets. Refer to [performance guidelines](https://hasindu2008.github.io/f5c/docs/f5c-perf-hints) for details.
* `--cuda-max-epk FLOAT`:       
  process the reads with events-per-kmer less than or equal to *cuda_max_epk* on GPU. The rest is processed on CPU [default value: 5.0]. Useful for tuning the CPU-GPU load balance for atypical datasets. Refer to [performance guidelines](https://hasindu2008.github.io/f5c/docs/f5c-perf-hints) for details.

#### advanced options

* `--skip-ultra FILE`:    
  skip ultra long reads and write those entries to the bam file provided as the argument
* `--ultra-thresh INT`:  
  threshold to skip ultra long reads [100000]
* `--skip-unreadable=yes|no`:  
  skip any unreadable fast5 or terminate program [yes]
* `--kmer-model FILE`:  
  custom nucleotide 6-mer model file (format similar to test/r9-models/r9.4_450bps.nucleotide.6mer.template.model)
* `--meth-model FILE`:  
  custom methylation 6-mer model file (format similar to test/r9-models/r9.4_450bps.cpg.6mer.template.model)
* `--meth-out-version INT`:  
  methylation tsv output version (set 2 to print the strand column) [1]
* `--cuda-mem-frac FLOAT`:  
  Fraction of free GPU memory to allocate [0.9 (0.7 for tegra)]


#### developer options

* `--print-events=yes|no`:  
  prints the event table
* `--print-banded-aln=yes|no`:  
  prints the event alignment
* `--print-scaling=yes|no`:  
  prints the estimated scalings
* `--print-raw=yes|no`:  
  prints the raw signal
* `--debug-break INT`:  
  break after processing the specified no. of batches
* `--profile-cpu=yes|no`:  
  process section by section (used for profiling on CPU)
* `--write-dump=yes|no`:  
  write the fast5 dump to a file or not
* `--read-dump=yes|no`:  
  read from a fast5 dump file or not
         
### meth-freq

`meth-freq [options...]`

```
  -c [float]        Call threshold. Default is 2.5.
  -i [file]         Input file. Read from stdin if not specified.
  -o [file]         Output file. Write to stdout if not specified.
  -s                Split groups
  ```


### eventalign

`f5c eventalign [OPTIONS] -r reads.fa -b alignments.bam -g genome.fa`

```
   -r FILE                    fastq/fasta read file
   -b FILE                    sorted bam file
   -g FILE                    reference genome
   -w STR[chr:start-end]      limit processing to genomic region STR
   -t INT                     number of threads [8]
   -K INT                     batch size (max number of reads loaded at once) [512]
   -B FLOAT[K/M/G]            max number of bases loaded at once [2.0M]
   -h                         help
   -o FILE                    output to file [stdout]
   --iop INT                  number of I/O processes to read fast5 files [1]
   --min-mapq INT             minimum mapping quality [30]
   --secondary=yes|no         consider secondary mappings or not [no]
   --verbose INT              verbosity level [0]
   --version                  print version
   --disable-cuda=yes|no      disable running on CUDA [no]
   --cuda-dev-id INT          CUDA device ID to run kernels on [0]
   --cuda-max-lf FLOAT        reads with length <= cuda-max-lf*avg_readlen on GPU, rest on CPU [3.0]
   --cuda-avg-epk FLOAT       average number of events per kmer - for allocating GPU arrays [2.0]
   --cuda-max-epk FLOAT       reads with events per kmer <= cuda_max_epk on GPU, rest on CPU [5.0]
   -x STRING                  profile to be used for optimal CUDA parameter selection. user-specified parameters will override profile values
advanced options:
   --kmer-model FILE          custom k-mer model file
   --skip-unreadable=yes|no   skip any unreadable fast5 or terminate program [yes]
   --print-events=yes|no      prints the event table
   --print-banded-aln=yes|no  prints the event alignment
   --print-scaling=yes|no     prints the estimated scalings
   --print-raw=yes|no         prints the raw signal
   --debug-break [INT]        break after processing the specified batch
   --profile-cpu=yes|no       process section by section (used for profiling on CPU)
   --skip-ultra FILE          skip ultra long reads and write those entries to the bam file provided as the argument
   --ultra-thresh [INT]       threshold to skip ultra long reads [100000]
   --write-dump=yes|no        write the fast5 dump to a file or not
   --read-dump=yes|no         read from a fast5 dump file or not
   --summary FILE             summarise the alignment of each read/strand in FILE
   --sam                      write output in SAM format
   --print-read-names         print read names instead of indexes
   --scale-events             scale events to the model, rather than vice-versa
   --samples                  write the raw samples for the event to the tsv output
   --cuda-mem-frac FLOAT      Fraction of free GPU memory to allocate [0.9 (0.7 for tegra)]
```

## EXAMPLES


* download and extract the dataset including sorted alignments:
  ```
  wget -O f5c_na12878_test.tgz "https://f5c.page.link/f5c_na12878_test"
  tar xf f5c_na12878_test.tgz
  ```
* index, call methylation and get methylation frequencies:
   ```
   f5c index -d chr22_meth_example/fast5_files chr22_meth_example/reads.fastq
   f5c call-methylation -b chr22_meth_example/reads.sorted.bam -g chr22_meth_example/humangenome.fa -r chr22_meth_example/reads.fastq > chr22_meth_example/result.tsv
   f5c meth-freq -i chr22_meth_example/result.tsv > chr22_meth_example/freq.tsv
   ```
* event alignment:
  ```
  f5c eventalign -b chr22_meth_example/reads.sorted.bam -g chr22_meth_example/humangenome.fa -r chr22_meth_example/reads.fastq > chr22_meth_example/events.tsv
  ```

## AUTHOR

Hasindu Gamaarachchi wrote the framework of f5c, CUDA code and integrated with adapted components from Jared T. Simpson's Nanopolish [https://github.com/jts/nanopolish], with tremendous support from Chun Wai Lam, Gihan Jayatilaka and Hiruna Samarakoon. 

## LICENSE

f5c is licensed under the MIT License. f5c reuses code and methods from [Nanopolish](https://github.com/jts/nanopolish) which is also under the MIT License. The event detection code in f5c is from Oxford Nanopore's [Scrappie basecaller](https://github.com/nanoporetech/scrappie) which is under Mozilla Public License 2.0. Some code snippets have been taken from [Minimap2](https://github.com/lh3/minimap2) and [Samtools](https://github.com/samtools/samtools) that are under the MIT License.

If you use f5c, please cite Gamaarachchi, H., Lam, C.W., Jayatilaka, G. et al. GPU accelerated adaptive banded event alignment for rapid comparative nanopore signal analysis. BMC Bioinformatics 21, 343 (2020). https://doi.org/10.1186/s12859-020-03697-x

## SEE ALSO

Full documentation: https://hasindu2008.github.io/f5c/docs/overview

Source code: https://github.com/hasindu2008/f5c/

Publication: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03697-x
