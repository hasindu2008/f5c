#!/bin/bash

test(){

testdir=$1

echo -e "chromosome\tstart\tend\tread_name\tlog_lik_ratio\tlog_lik_methylated\tlog_lik_unmethylated\tnum_calling_strands\tnum_cpgs\tsequence" > $testdir/meth.txt
cat $testdir/meth.exp >> $testdir/meth.txt

#~/nanopolish/scripts/calculate_methylation_frequency.py -i $testdir/meth.txt >  $testdir/freq.exp

./f5c meth-freq -i $testdir/meth.txt > $testdir/freq.txt
diff $testdir/freq.exp $testdir/freq.txt

}

valgrind(){

testdir=$1
valgrind --leak-check=full --error-exitcode=1 ./f5c meth-freq -i $testdir/meth.txt > /dev/null
}

testdir=test/chr22_meth_example/
test($testdir)


testdir=test/ecoli_2kb_region/
test($testdir)
valgrind($testdir)