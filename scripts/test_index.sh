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

	# if [ $testdir = test/chr22_meth_example ]; then
    #     cat $testdir/meth.exp > $testdir/meth.txt
    # else

    # fi
    echo ""
    echo "testing with iop options"
    ./f5c index -d $testdir/fast5_files $testdir/reads.fasta
    mv $testdir/reads.fasta.index.readdb $testdir/reads.fasta.index.readdb2
    ./f5c index -d $testdir/fast5_files $testdir/reads.fasta --iop 8
    diff $testdir/reads.fasta.index.readdb $testdir/reads.fasta.index.readdb2  || die "verification failed"

    rm $testdir/reads.fasta.index.readdb2
    echo "______________________________________________________________"
    echo ""
    echo "valgrind check"
    valgrind --leak-check=full --error-exitcode=1 ./f5c index -d $testdir/fast5_files $testdir/reads.fasta --iop 8 || die "valgrind verification failed"



}

help_msg() {
	echo "Test script for f5c."
	echo "Usage: f5c_dir/script/test_index.sh [-c]"
	echo "-h                   Show this help message."
}

testdir=test/ecoli_2kb_region

# parse options
while getopts ch opt
do
	case $opt in
		#c) testdir=test/chr22_meth_example;;
        h) help_msg
		   exit 0;;
		?) printf "Usage: %s [-c]\n" "$0"
		   exit 2;;
	esac
done

test "$testdir"
