#!/bin/awk -f
function abs(x) {
	return (((x) < 0.0) ? -(x) : (x))
}

FNR == 1 {next}

# to be used with join
BEGIN {status=0}
{	if ($2!=$(2+13) \
	|| $3!=$(3+13) \
	|| $4!=$(4+13) \
	|| $5!=$(5+13) \
	|| $6!=$(6+13) \
	|| $7!=$(7+13) \
	|| $8!=$(8+13) \
	|| $9!=$(9+13) \
	|| abs($10-$(10+13)) > 0.01 \
	|| abs($11-$(11+13)) > 0.01 \
	|| abs($12-$(12+13)) > 0.01 \
	|| $13!=$(13+13) \
	|| abs($14-$(14+13)) > 0.01 )
	{print "f5c - nanopolish mismatch at index " $1; status=1}}
END {if (status > 0) {exit 1}}

#1 check read_index
#2 check read_name
#3 check fast5_path
#4 check model_name
#5 check strand
#6 check num_events
#7 check num_steps
#8 check num_skips
#9 check num_stays
#10 check sum_duration
#11 check scalings shift
#12 check scalings scale
#13 check scalings drift
#14 check scalings var
