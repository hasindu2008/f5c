#!/usr/bin/env bash

die() {
    local msg="${1}"
    echo "run.sh: ${msg}" >&2
    exit 1
}

set -x

#export LD_PRELOAD=/home/hasindu/installs/gperftools/.libs/libtcmalloc.so


nanopath_cpp_vec="/home/hasindu/nano_vector_test/nanopolish-cppvector/nanopolish"
nanopath_1d="/home/hasindu/nano_vector_test/nanopolish-1d/nanopolish"
nanopath_orig="/home/hasindu/nano_vector_test/nanopolish_orig/nanopolish"
testdir="test/chr22_meth_example/"

test -d $testdir || die "$testdir does not exist"

bamfile=$testdir/"reads.sorted.bam"
ref=$testdir/"humangenome.fa"
reads=$testdir/"reads.fastq"


for file in "${bamfile}" "${ref}" "${reads}"; do
    [[ -f "${file}" ]] || die "${file}: File does not exist"
done

t=64
/usr/bin/time -v "${nanopath_cpp_vec}" call-methylation -b "${bamfile}" -g "${ref}" -r "${reads}" -t $t  -K4096 > cpp.tsv
/usr/bin/time -v "${nanopath_1d}" call-methylation -b "${bamfile}" -g "${ref}" -r "${reads}" -t $t  -K4096 > oned.tsv
/usr/bin/time -v "${nanopath_orig}" call-methylation -b "${bamfile}" -g "${ref}" -r "${reads}" -t $t  -K4096 > orig.tsv
sort  cpp.tsv > cppsorted.tsv
sort oned.tsv > onedsorted.tsv
sort orig.tsv > origsorted.tsv

diff -q origsorted.tsv onedsorted.tsv || die "diff failed"
diff -q origsorted.tsv cppsorted.tsv || die "diff failed"

#--ploidy=2 -b test/chr22_meth_example//reads.sorted.bam -g test/chr22_meth_example//humangenome.fa -r test/chr22_meth_example//reads.fastq -w chr20

for t in 8 16 32 64
do
    /usr/bin/time -v "${nanopath_cpp_vec}" call-methylation -b "${bamfile}" -g "${ref}" -r "${reads}" -t $t  -K4096 > result.txt 2> nano_cpp_$t.log
done


for t in 8 16 32 64
do
    /usr/bin/time -v "${nanopath_1d}" call-methylation -b "${bamfile}" -g "${ref}" -r "${reads}" -t $t  -K4096 > result.txt 2> nano_1d_$t.log
done


for t in 8 16 32 64
do
    /usr/bin/time -v "${nanopath_orig}" call-methylation -b "${bamfile}" -g "${ref}" -r "${reads}" -t $t  -K4096 > result.txt 2> nano_orig_$t.log
done
