#!/bin/bash

# terminate script
die() {
	echo "$1" >&2
	echo
	exit 1
}

test(){
    testdir=$1

    echo "Test for $testdir"

	if [ $testdir = test/chr22_meth_example ]; then
        cat $testdir/meth.exp > $testdir/meth.txt
    else
        echo -e "chromosome\tstart\tend\tread_name\tlog_lik_ratio\tlog_lik_methylated\tlog_lik_unmethylated\tnum_calling_strands\tnum_cpgs\tsequence" > $testdir/meth.txt
        cat $testdir/meth.exp >> $testdir/meth.txt

        echo "valgrind check"
        valgrind --leak-check=full --error-exitcode=1 ./f5c meth-freq -i $testdir/meth.txt > $testdir/freq.txt || die "valgrind verification failed"

        echo ""
        echo "testing with cpg group split -s"
        ./f5c meth-freq -i $testdir/meth.txt -s > $testdir/freq.txt
        diff $testdir/freq-s.exp $testdir/freq.txt || die "cpg group split - verification failed"
    fi
    echo ""
    echo "testing with default options"
    ./f5c meth-freq -i $testdir/meth.txt > $testdir/freq.txt
    diff $testdir/freq.exp $testdir/freq.txt || die "verification failed"

    echo "______________________________________________________________"
    echo ""
}

help_msg() {
	echo "Test script for f5c."
	echo "Usage: f5c_dir/script/test_methfreq.sh [-c]"
	echo "-h                   Show this help message."
}

testdir=test/ecoli_2kb_region

# parse options
while getopts ch opt
do
	case $opt in
		c) testdir=test/chr22_meth_example;;
        h) help_msg
		   exit 0;;
		?) printf "Usage: %s [-c]\n" "$0"
		   exit 2;;
	esac
done

test "$testdir"
