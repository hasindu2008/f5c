#!/bin/awk -f
function abs(x) {
	return (((x) < 0.0) ? -(x) : (x))
}

FNR == 1 {next}

BEGIN {status=0}
{	if ($1!=$12\
	|| abs($2-$13) > 2 \
	|| $3!=$14 \
	|| $4!=$15 \
	|| $5!=$16 \
	|| $6!=$17 \
	|| $7!=$18 \
	|| $9!=$20 \
	|| $10!=$21 \
	|| abs($11-$22) > 0.02) 
	{print "f5c - nanopolish mismatch at line " NR; status=1}}
END {if (status > 0) {exit 1}}


#check contig	[tig00000001]
#check ref_position	[15]  
#check ref_kmer	AACGCA
#check strand	[t]
#check event_index	[10711]
#check event_level_mean	[104.41]
#check event_stdv	[0.598]
#check model_kmer --not checked	[AACGCA]
#check model_mean	[100.57]
#check model_stdv	[2.83]
#check standardized_level	[1.15]
