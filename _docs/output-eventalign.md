# eventalign output

The default output is as below:

|Col|Type  |Name |Description                               |
|--:|:----:|:----|:-----------------------------------------|
|1  |string|contig|Contig on the reference which the read maps to                       |
|2  |int   |position|Start index on the contig (0-based; BED-like; closed)                    |
|3  |string   |reference_kmer|the k-mer on the reference   |
|4  |int   |read_index|Index in the BAM file for the coresponding read (0-based)       |
|5  |char  |strand|Legacy, ignore               |
|6  |int|event_index|Index of the event on the event table (0-based; BED-like; closed)                    |
|7  |float   |event_level_mean|Mean level of the current values of the event                   |
|8  |float   |event_stdv| Standard deviation of the current values of the event  |
|9  |float   |event_length| Length of the event (in seconds)   |
|10 |int   |model_kmer|The k-mers on the pore-model which this event matched                   |
|11 |int   |model_mean|Scaled mean level on the pore-model for the matched k-mer (scaling.scale * level_mean + scaling.shift)                    |
|12 |int   |model_stdv|Scaled standard deviation on the pore-model for the matched k-mer (level_stdv * scaling.var)  |
|13 |int   |standardized_level|(event_level_mean - model_mean) / (sqrt(scalings.var) * model_stdv)  |

Following options can be used to modify default columns or print additional columns:

- `--print-read-name`: Column 4 will become read_name that prints the read ID.
- `--signal-index`: two additional columns will be printed, namely start_idx and end_idx. </br> 
   start_idx is the starting index on the raw signal which the corresponding k-mer maps to (0-based; BED-like; closed). </br> 
   end_idx is the ending index on the raw signal which the corresponding k-mer maps to (0-based; BED-like; open).
- `--samples`:   prints the comma separated signal samples corresponding to the mapped k-mer (scaled pA current values) </br>
    `scaled pA current values = (pA - scaling.shift) / scaling.scale` where `pA = (raw_signal + offset) * range / digitisation`
- `--scale-events`:   intead of scaling the model to the events, now events will be scaled to the model
    column 7 becomes `(event_level_mean-scaling.shift)/scaling.scale`
    column 11 becomes level_mean and column 12 becomes model_stdv

