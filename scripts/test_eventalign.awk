#!/bin/awk -f
function abs(x) {
	return (((x) < 0.0) ? -(x) : (x))
}

FNR == 1 {next}

BEGIN {status=0}
{	if ($2!=$11 \
	|| $3!=$12 \
	|| $4!=$13 \
	|| $5!=$14 \
	|| abs($6-$15) > 0.01\
	|| abs($7-$16) > 0.01\
	|| abs($8-$17) > 0.01 \
	|| $9!=$18 \
	|| abs($10-$19) > 0.01 )
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
