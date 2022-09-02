# Output

f5c call-methylation, eventalign and meth-freq have the same output format as nanopolish and thus I have not put much effort to document those formats.

## resqiiggle

f5c resquiggle which is under development is explained below. Note that as this is under development, the output format is a draft and will change in the future versions. When it is stable, this notice will be removed.

The default output is an intuitive TSV format with the following columns.

|Col|Type  |Name            |Description                                                            |
|--:|:----: |:------:        |:-----------------------------------------                             |
|1  |read_id|read_id         |Read identifier name                                                   |
|2  |int    |kmer_idx  |k-mer index on the basecalled read (0-based)                               |
|3  |int    |start_raw_idx     |Raw signal start index for the corresponding k-mer (0-based; BED-like; closed)                                  |
|4  |int    |end_raw_idx       |Raw signal end index for the corresponding k-mer (0-based; BED-like; open)                                   |

If a corresponding base has no corresponding signal samples, a '.' will be printed.


The above output is bulky. Specifying `-c` will produce a condensed output in PAF-like format (inspired by UNCALLED) with the following columns:
The query is the raw-signal and the target is the basecalled-read.

|Col|Type  |Description                               |
|--:|:----:|:-----------------------------------------|
|1  |string|Read identifier name                       |
|2  |int   |Raw signal length (number of samples)                    |
|3  |int   |Raw signal start index  (0-based; BED-like; closed)   |
|4  |int   |Raw signal end index (0-based; BED-like; open)       |
|5  |char  |Relative strand: "+" or "-"               |
|6  |string|Same as colum 1                     |
|7  |int   |base-called sequence length                    |
|8  |int   |start on basecalled sequence (0-based; BED-like; closed)  |
|9  |int   |end on basecalled sequence (0-based; BED-like; open)   |
|10 |int   |Number of residue matches                 |
|11 |int   |Alignment block length                    |
|12 |int   |Mapping quality (0-255; 255 for missing)  |

10,11 and 12 are to be decided.

Following optional tags are present:

|Tag|Type  |Description                               |
|--:|:----:|:-----------------------------------------|
|sc  |f| Post alignment  recalibrated scale parameter                     |
|sh  |f   |Post alignment reclibrated shift parameter                      |
|ss  |Z   |signal alignment string in format described below   |

*ss* string is a custom encoding that compacts the signal-base alignment. It can be thought of as an extended version CIGAR string that accommodates for signal alignment needs.

Consider the following example:

8,5,4,8I4,3D4,5

This means 8 signal samples map to the starting base of the sequence; the next 5 samples to the next base, the next 4 samples to the next base; 8 next samples are missing a mapping in the basecalled read (insertion to reference); 4 samples map to the next base; 3 bases in the basecalled read have no corresponding signal samples (deletion); 4 bases map to the next base; and 5 bases map to the next base.

Note that the start indexes of read and signal are the absolute values in columns 8 and 3 above. the ss string is relative to this.
, after a number means step one base in the basecall reference, while step number of samples preceding the , in the raw signal.


TODO: 
- For RNA, column 3 is large than 4. That is because the base left-end of the read map to the right-end of the signal. ss string is relative to this right-end coordinate and is computed leftwards.

- - start and -end in paf

