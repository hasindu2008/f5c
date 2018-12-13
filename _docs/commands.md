---
title: Commands and options
---

## Indexing

```sh
Usage: f5c index [OPTIONS] -d nanopore_raw_file_directory reads.fastq
Build an index mapping from basecalled reads to the signals measured by the sequencer
f5c index is equivalent to nanopolish index by Jared Simpson

  -h, --help                           display this help and exit
  -v, --verbose                        display verbose output
  -d, --directory                      path to the directory containing the raw ONT signal files. This option can be given multiple times.
  -s, --sequencing-summary             the sequencing summary file from albacore, providing this option will make indexing much faster
  -f, --summary-fofn                   file containing the paths to the sequencing summary files (one per line)
```

## Calling methylation

```sh
Usage: f5c call-methylation [OPTIONS] -r reads.fa -b alignments.bam -g genome.fa
   -r FILE                 fastq/fasta read file
   -b FILE                 sorted bam file
   -g FILE                 reference genome
   -t INT                  number of threads [8]
   -K INT                  batch size (number of reads loaded at once) [512]
   -h                      help
   --min-mapq INT          minimum mapping quality [30]
   --secondary             consider secondary mappings or not [no]
   --skip-unreadable       skip any unreadable fast5 or terminate program [yes]
   --verbose INT           verbosity level [0]
   --version               print version
   --disable-cuda          disable running on CUDA [no] (only if compiled for CUDA)
debug options:
   --kmer-model FILE       custom k-mer model file (used for debugging)
   --print-events          prints the event table (used for debugging)
   --print-banded-aln      prints the event alignment (used for debugging)
   --print-scaling         prints the estimated scalings (used for debugging)
   --print-raw             prints the raw signal (used for debugging)
   --debug-break           break after processing the first batch (used for debugging)
```
