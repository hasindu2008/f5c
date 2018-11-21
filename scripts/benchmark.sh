#!/usr/bin/env bash

die() {
    local msg="${1}"
    echo "run.sh: ${msg}" >&2
    exit 1
}

set -x

#export LD_PRELOAD=/home/hasindu/installs/gperftools/.libs/libtcmalloc.so

f5cpath="./f5c"
nanopath="/home/hasindu/hasindu2008.git/nanopolish_bench/nanopolish"
testdir="test/chr22_meth_example/"

test -d $testdir || die "$testdir does not exist"

bamfile=$testdir/"reads.sorted.bam"
ref=$testdir/"humangenome.fa"
reads=$testdir/"reads.fastq"


for file in "${bamfile}" "${ref}" "${reads}"; do
    [[ -f "${file}" ]] || die "${file}: File does not exist"
done

t=64
/usr/bin/time -v "${f5cpath}" -b "${bamfile}" -g "${ref}" -r "${reads}" -t $t --secondary=yes --min-mapq=0 --print-scaling=yes -K4096 > /dev/null


for t in 8 16 24 32 48 64
do
     /usr/bin/time -v "${f5cpath}" -b "${bamfile}" -g "${ref}" -r "${reads}" -t $t --secondary=yes --min-mapq=0 --print-scaling=yes -K4096 > result.txt 2> f5c_$t.log
done

t=64
/usr/bin/time -v "${nanopath}" call-methylation -b "${bamfile}" -g "${ref}" -r "${reads}" -t $t  -K4096 > /dev/null
for t in 8 16 24 32 48 64
do
    /usr/bin/time -v "${nanopath}" call-methylation -b "${bamfile}" -g "${ref}" -r "${reads}" -t $t  -K4096 > result.txt 2> nano_$t.log
done
