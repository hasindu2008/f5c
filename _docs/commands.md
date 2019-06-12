---
title: Commands and options
---

## Available f5c tools

```
Usage: f5c <command> [options]

command:
         index               Build an index mapping from basecalled reads to the signals measured by the sequencer (same as nanopolish index)
         call-methylation    Classify nucleotides as methylated or not (optimised nanopolish call-methylation)
         meth-freq           Calculate methylation frequency at genomic CpG sites (optimised nanopolish calculate_methylation_frequency.py)
         eventalign          Align nanopore events to reference k-mers (optimised nanopolish eventalign)
```

### Indexing

```
Usage: f5c index [OPTIONS] -d nanopore_raw_file_directory reads.fastq
Build an index mapping from basecalled reads to the signals measured by the sequencer
f5c index is equivalent to nanopolish index by Jared Simpson

  -h, --help                           display this help and exit
  -v, --verbose                        display verbose output
  -d, --directory                      path to the directory containing the raw ONT signal files. This option can be given multiple times.
  -s, --sequencing-summary             the sequencing summary file from albacore, providing this option will make indexing much faster
  -f, --summary-fofn                   file containing the paths to the sequencing summary files (one per line)
```

### Calling methylation

```
Usage: f5c call-methylation [OPTIONS] -r reads.fa -b alignments.bam -g genome.fa
   -r FILE                    fastq/fasta read file
   -b FILE                    sorted bam file
   -g FILE                    reference genome
   -t INT                     number of threads [8]
   -K INT                     batch size (max number of reads loaded at once) [512]
   -B FLOAT[K/M/G]            max number of bases loaded at once [2.0M]
   -h                         help
   --min-mapq INT             minimum mapping quality [30]
   --secondary=yes|no         consider secondary mappings or not [no]
   --skip-unreadable=yes|no   skip any unreadable fast5 or terminate program [yes]
   --verbose INT              verbosity level [0]
   --version                  print version
   --disable-cuda=yes|no      disable running on CUDA [no]
   - cuda-dev-id INT          CUDA device ID to run kernels on [0]
   --cuda-max-lf FLOAT        reads with length <= cuda-max-lf*avg_readlen on GPU, rest on CPU [3.0]
   --cuda-avg-epk FLOAT       average number of events per kmer - for allocating GPU arrays [2.0]
   --cuda-max-epk FLOAT       reads with events per kmer <= cuda_max_epk on GPU, rest on CPU [5.0]
advanced options:
   --kmer-model FILE          custom k-mer model file
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
   --cuda-mem-frac FLOAT      fraction of free GPU memory to allocate [0.9 (0.7 for tegra)]
```

### Calculate methylation frequency (experimental, not thoroughly tested)
```
Usage: meth-freq [options...]

  -c [float]        Call threshold. Default is 2.5.
  -i [file]         Input file. Read from stdin if not specified.
  -s                Split groups
  ```


### Aligning events (experimental, not thoroughly tested)

```
Usage: f5c eventalign [OPTIONS] -r reads.fa -b alignments.bam -g genome.fa
   -r FILE                    fastq/fasta read file
   -b FILE                    sorted bam file
   -g FILE                    reference genome
   -t INT                     number of threads [8]
   -K INT                     batch size (max number of reads loaded at once) [512]
   -B FLOAT[K/M/G]            max number of bases loaded at once [2.0M]
   -h                         help
   --min-mapq INT             minimum mapping quality [30]
   --secondary=yes|no         consider secondary mappings or not [no]
   --skip-unreadable=yes|no   skip any unreadable fast5 or terminate program [yes]
   --verbose INT              verbosity level [0]
   --version                  print version
   --disable-cuda=yes|no      disable running on CUDA [no]
   - cuda-dev-id INT          CUDA device ID to run kernels on [0]
   --cuda-max-lf FLOAT        reads with length <= cuda-max-lf*avg_readlen on GPU, rest on CPU [3.0]
   --cuda-avg-epk FLOAT       average number of events per kmer - for allocating GPU arrays [2.0]
   --cuda-max-epk FLOAT       reads with events per kmer <= cuda_max_epk on GPU, rest on CPU [5.0]
advanced options:
   --kmer-model FILE          custom k-mer model file
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
   --cuda-mem-frac FLOAT      fraction of free GPU memory to allocate [0.9 (0.7 for tegra)]
```
