#!/bin/awk -f
function abs(x) {
	return (((x) < 0.0) ? -(x) : (x))
}

FNR == 1 {next}

BEGIN {status=0}
{	if ($1!=$(1+15)\
	|| abs($2-$(2+15)) > 2 \
	|| $3!=$(3+15) \
	|| $4!=$(4+15) \
	|| $5!=$(5+15) \
	|| $6!=$(6+15) \
	|| $7!=$(7+15) \
	|| $8!=$(8+15) \
	|| $9!=$(9+15) \
	|| $10!=$(10+15) \
	|| $11!=$(11+15) \
	|| $12!=$(12+15) \
	|| abs($13-$(13+15)) > 0.02 \
	|| $14!=$(14+15) \
	|| $15!=$(15+15)) \
	{print "f5c - nanopolish mismatch at line " NR; status=1}}
END {if (status > 0) {exit 1}}


#1 check contig	[tig00000001]
#2 check ref_position	[15]
#3 check ref_kmer	AACGCA
#4 check read_index [0]
#5 check strand	[t]
#6 check event_index	[10711]
#7 check event_level_mean	[104.41]
#8 check event_stdv	[0.598]
#9 check event_length [0.001]
#10 check model_kmer [AACGCA]
#11 check model_mean	[100.57]
#12 check model_stdv	[2.83]
#13 check standardized_level	[1.15]
#14 start_idx [250]
#15 end_idx [254]
