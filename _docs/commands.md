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
         freq-merge          Merge calculated methylation frequency tsv files (new feature)
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

### Calculate methylation frequency
```
Usage: meth-freq [options...]

  -c [float]        Call threshold. Default is 2.5.
  -i [file]         Input file. Read from stdin if not specified.
  -o [file]         Output file. Write to stdout if not specified.
  -s                Split groups
  ```


### Aligning events

```
Usage: f5c eventalign [OPTIONS] -r reads.fa -b alignments.bam -g genome.fa
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

### Merge methylation frequency calculated tsv files
```
Usage: freq-merge [options...]
To merge multiple methylation frequency files to a single file.
For each methylation calling output (.tsv) file, perform methylation frequency calculation separately (no concatenation required). 
Then feed those output (.tsv) files to this tool, to obtain the final methylation frequency calculated file. 

  -o [float]        Output file.
  -n [INT]          Number of methylation frequency .tsv files to be merged
  -f                n number of input filepaths should be followed
  e.g. freq-merge -o merged_freq.tsv -n 2 -f data1_freq.tsv data2_freq.tsv 
  ```
