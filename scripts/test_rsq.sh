#!/bin/sh

# exit when command fails
set -e

. scripts/common.sh

# defaults
exepath=./f5c
testdir=test/ecoli_2kb_region

ref=${testdir}/draft.fa
reads=${testdir}/reads.fasta
slow5=${testdir}/reads.blow5
batchsize=256
max_bases=2M
if command -v nproc > /dev/null; then
	threads=$(nproc --all)
else
	threads=8
fi
# execution mode (valgrind/gdb/cpu/cuda/echo)
mode=
testset_url="https://f5c.page.link/f5c_ecoli_2kb_region_test"
fallback_url="https://f5c.page.link/f5c_ecoli_2kb_region_test_fallback"

# download test set given url
#
# $1 - URL
# $2 - Fallback URL (optional)
download_test_set() {
	# data set exists
	if [ -d ${testdir} ]; then
		return
	fi

	mkdir -p test
	tar_path=test/data.tgz
	wget -O $tar_path "$1" || wget -O $tar_path "$2" || rm -rf $tar_path ${testdir}
	echo "Extracting. Please wait."
	tar -xf $tar_path -C test || rm -rf $tar_path ${testdir}
	rm -f $tar_path
}



handle_tests() {
	numfailed=$(wc -l < ${testdir}/floatdiff.txt)
	numcases=$(wc -l < ${testdir}/meth_float.txt)
	numres=$(wc -l < ${testdir}/result_float.txt)
	echo "$numfailed of $numcases test cases deviated."
	missing=$(echo "$numcases-$numres" | bc)
	echo "$missing entries in the truthset are missing in the testset"
	failp=$(echo "$numfailed*100/$numcases" | bc)
	[ "$failp" -gt 5 ] && die "${1}: Validation failed"
	echo "Validation passed"
}

execute_test() {

	if [ $testdir = test/chr22_meth_example ]; then
		echo "A"
	elif [ $testdir = test/ecoli_2kb_region ]; then
		diff -q ${testdir}/result.txt ${testdir}_big_testresults/resquiggle.tsv || die "Validation of tsv failed"
		diff -q ${testdir}/result2.txt ${testdir}_big_testresults/resquiggle.paf || die "Validation of paf failed"
	elif [ $testdir = test/rna ]; then
		diff -q ${testdir}/result.txt ${testdir}/resquiggle.tsv || die "Validation of tsv failed"
		diff -q ${testdir}/result2.txt ${testdir}/resquiggle.paf || die "Validation of paf failed"
	fi

}

mode_test() {


	cmd="${exepath} resquiggle ${reads} ${slow5} -t ${threads} -K $batchsize -B $max_bases"

	case $1 in
		valgrind) valgrind --leak-check=full --error-exitcode=1 $cmd > /dev/null || die "valgrind failed" ;;
		gdb) gdb --args "$cmd";;
		cpu) $cmd --disable-cuda=yes > ${testdir}/result.txt; execute_test;;
		cuda) $cmd --disable-cuda=no > ${testdir}/result.txt; execute_test;;
		echo) echo "$cmd  > ${testdir}/result.txt";;
		nvprof) nvprof  -f --analysis-metrics -o profile.nvprof "$cmd" --disable-cuda=no --debug-break=5 > /dev/null;;
		custom) shift; $cmd "$@" > ${testdir}/result.txt; execute_test;;
		*) die "Unknown mode: $1";;
	esac
}

help_msg() {
	echo "Test script for f5c."
	echo "Usage: f5c_dir/script/test.sh [-c/-e] [-r fastq/fasta read] mode"
	echo
	echo "mode                 one of: valgrind/gdb/cpu/cuda/echo"
	echo
	echo "-c                   Uses chr22_meth_example test set."
	echo "-e                   Uses rna test set."
	echo "-K [n]               Same as f5c -K."
	echo "-B [n]               Same as f5c -B."
	echo "-r [read file]       Same as f5c -r."
	echo "-t [n]               Number of threads."
	echo "-d                   Download chr22_meth_example data set and exit."
	echo "-h                   Show this help message."
}

# parse options
while getopts b:g:r:t:K:B:cdhe opt
do
	case $opt in
		r) reads="$OPTARG";;
		t) threads="$OPTARG";;
		c) testdir=test/chr22_meth_example
		   reads=${testdir}/reads.fastq
		   slow5=${testdir}/reads.blow5
		   testset_url="https://f5c.page.link/f5c_na12878_test"
		   fallback_url="https://f5c.page.link/f5c_na12878_test_fallback";;
		e) testdir=test/rna
		   reads=${testdir}/reads.fastq
		   slow5=${testdir}/reads.blow5
		   testset_url="https://f5c.page.link/f5c_rna_test"
		   fallback_url="https://f5c.page.link/f5c_rna_test_fallback";;
		K) batchsize="$OPTARG";;
		B) max_bases="$OPTARG";;
		d) download_test_set "https://f5c.page.link/f5c_na12878_test" "https://f5c.page.link/f5c_na12878_test_fallback"
		   exit 0;;
		h) help_msg
		   exit 0;;
		?) printf "Usage: %s [-c] [-b bam file] [-g reference genome] [-r fastq/fasta read] args" "$0"
		   exit 2;;
	esac
done
shift $(($OPTIND - 1))
mode=$1

download_test_set $testset_url $fallback_url

# validate files
for file in ${reads} ${slow5}; do
	[ -f ${file} ] || die "${file}: File does not exist"
done

if [ -z "$mode" ]; then

	${exepath} index --slow5 ${slow5} ${reads} -t "$threads"
	${exepath} resquiggle ${reads} ${slow5} -t ${threads} -K $batchsize -B $max_bases > ${testdir}/result.txt
	${exepath} resquiggle ${reads} ${slow5} -t ${threads} -K $batchsize -B $max_bases -c > ${testdir}/result2.txt
	execute_test

else
	mode_test "$@"
fi
