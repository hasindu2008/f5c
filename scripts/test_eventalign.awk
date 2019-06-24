#!/bin/awk -f
function abs(x) {
	return (((x) < 0.0) ? -(x) : (x))
}

FNR == 1 {next}

BEGIN {status=0}
{	if ($1!=$14\
	|| abs($2-$15) > 2 \
	|| $3!=$16 \
	|| $4!=$17 \
	|| $5!=$18 \
	|| $6!=$19 \
	|| $7!=$20 \
	|| abs($8-$21) > 1 \
	|| abs($11-$23) > 0.02 \
	|| $12!=$24 \
	|| abs($13-$25) > 0.02) 
	{print "f5c - nanopolish mismatch at index " $1; status=1}}
END {if (status > 0) {exit 1}}


#check contig	[tig00000001]
#check ref_position	[15]  
#check ref_kmer	AACGCA
#check strand	[t]
#check event_index	[10711]
#check event_level_mean	[104.41]
#check event_stdv	[0.598]
#check event_length	[4]
#check event_duration -- not checked	[0.00100]
#check model_kmer --not checked	[AACGCA]
#check model_mean	[100.57]
#check model_stdv	[2.83]
#check standardized_level	[1.15]
