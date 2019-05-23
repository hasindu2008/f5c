#!/bin/awk -f
function abs(x) {
	return (((x) < 0.0) ? -(x) : (x))
}

FNR == 1 {next}

BEGIN {status=0}
{	if ($2!=$13 \ 
	|| $3!=$14 \ 
	|| $4!=$15 \ 
	|| $5!=$16 \ 
	|| $6!=$17 \ 
	|| $7!=$18 \ 
	|| abs($8-$19) > 0.01\ 
	|| abs($9-$20) > 0.01\ 
	|| abs($10-$21) > 0.01 \ 
	|| $11!=$22 \ 
	|| abs($12-$23) > 0.01 )
	{print "f5c - nanopolish mismatch at index " $1; status=1}}
END {if (status > 0) {exit 1}}

#check uuids
#check template 
#check num_events
#check num_steps
#check num_skips
#check num_stays
#check sum_duration
#check scalings shift
#check scalings scale
#check scalings drift
#check scalings var