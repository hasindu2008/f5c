---
title: Commands and Options
---

## NAME

f5c(1) - Ultra-fast methylation calling and event alignment tool for nanopore sequencing data with optional CUDA acceleration

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

Given a set of base-called nanopore reads and associated raw signals, f5c call-methylation detects the methylated cytosine and f5c eventalign aligns raw nanopore DNA signals (events) to the base-called read. f5c can optionally utilise CUDA enabled NVIDIA graphics cards for acceleration. f5c is a heavily re-engineered and optimised implementation of the call-methylation and eventalign modules in Nanopolish.

## COMMANDS 

* `index`:               
         Build an index mapping from basecalled reads to the signals measured by the sequencer (same as nanopolish index)
* `call-methylation`:    
         Classify nucleotides as methylated or not (optimised nanopolish call-methylation)
* `meth-freq`:           
         Calculate methylation frequency at genomic CpG sites (optimised nanopolish calculate_methylation_frequency.py)
* `eventalign`:          
         Align nanopore events to reference k-mers (optimised nanopolish eventalign)


## OPTIONS

### index

`f5c index [OPTIONS] -d nanopore_raw_file_directory reads.fastq`

Build an index mapping from basecalled reads to the signals measured by the sequencer. f5c index is equivalent to nanopolish index by Jared Simpson.


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
   --meth-out-version [INT]   methylation tsv output version (set 2 to print the strand column) [1]
   --cuda-mem-frac FLOAT      Fraction of free GPU memory to allocate [0.9 (0.7 for tegra)]
```

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

f5c is licensed under the MIT License. f5c reuses code and methods from Nanopolish [https://github.com/jts/nanopolish] which is also under the MIT License. The event detection code in f5c is from Oxford Nanopore's Scrappie basecaller [https://github.com/nanoporetech/scrappie] which is under Mozilla Public License 2.0. Some code snippets have been taken from Minimap2 [https://github.com/lh3/minimap2] and Samtools [https://github.com/samtools/samtools] that are under the MIT License.

If you use f5c, please cite Gamaarachchi, H., Lam, C.W., Jayatilaka, G. et al. GPU accelerated adaptive banded event alignment for rapid comparative nanopore signal analysis. BMC Bioinformatics 21, 343 (2020). https://doi.org/10.1186/s12859-020-03697-x

## SEE ALSO

Full documentation: https://hasindu2008.github.io/f5c/docs/overview

Source code: https://github.com/hasindu2008/f5c/

Publication: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03697-x
