#!/bin/sh

# exit when command fails
set -e

. scripts/common.sh

# defaults
exepath=./f5c
testdir=test/ecoli_2kb_region

bamfile=${testdir}/reads.sorted.bam
ref=${testdir}/draft.fa
reads=${testdir}/reads.fasta
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


#############################################################################################

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


download_ecoli_2kb_region_big_testresults() {
	mkdir -p test
	tar_path=test/data.tgz
	wget -O $tar_path "https://f5c.page.link/f5c_ecoli_2kb_region_big_testresults" || rm -rf $tar_path ${testdir}_big_testresults
	echo "Extracting. Please wait."
	tar -xf $tar_path -C test || rm -rf $tar_path ${testdir}_big_testresults
	rm -f $tar_path
}

#event alignsummaries
handle_tests() {
	numfailed=$(cat  ${testdir}/joined_diff.txt | awk '{print $NF}' | sort -u -k1,1 | wc -l)
	numcases=$(wc -l < ${testdir}/nanopolish.summary.txt)
	numres=$(wc -l < ${testdir}/f5c.summary.txt)
	echo "$numfailed of $numcases test cases deviated."
	missing=$(echo "$numcases-$numres" | bc)
	echo "$missing entries in the truthset are missing in the testset"
	failp=$(echo "$numfailed*100/$numcases" | bc)
	[ "$failp" -gt 5 ] && die "${1}: Validation failed"
	echo "Validation passed"
}

#full eventalign results
handle_tests2() {
	numfailed=$(cat  ${testdir}/joined_diff.txt |  wc -l)
	numcases=$(wc -l < ${testdir}/nanopolish.txt)
	numres=$(wc -l < ${testdir}/f5c.txt)
	echo "$numfailed of $numcases test cases deviated."
	missing=$(echo "$numcases-$numres" | bc)
	echo "$missing entries in the truthset are missing in the testset"
	failp=$(echo "$numfailed*100/$numcases" | bc)
	[ "$failp" -gt 5 ] && die "${1}: Validation failed"
	echo "Validation passed"
}


execute_test() {
	echo "----------------comparing summaries-----------------"
	tail -n +2 ${testdir}/eventalign.summary.exp | awk '{print $1"\t"$2"\t"$3"\tna\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}' > ${testdir}/nanopolish.summary.txt
	tail -n +2 ${testdir}/f5c_event_align.summary.txt | awk '{print $1"\t"$2"\t"$3"\tna\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14}' > ${testdir}/f5c.summary.txt
	#this test still passes if f5c.summary.txt is empty, but will be caught when comparing full summaries
	join ${testdir}/nanopolish.summary.txt ${testdir}/f5c.summary.txt > ${testdir}/joined_results.txt || echo "Join ran into an issue. Probably just a warning."
	awk -f  scripts/test_eventalign_summary.awk ${testdir}/joined_results.txt > ${testdir}/joined_diff.txt || handle_tests "${file}"

	if [ $testdir = test/chr22_meth_example ]; then
		echo "event by event test not implemented not yet implemented!"
	else
		echo "----------------comparing full results--------------"
		if [ $testdir = test/ecoli_2kb_region ]; then
			test -d ${testdir}_big_testresults || download_ecoli_2kb_region_big_testresults
			tail -n +2 ${testdir}_big_testresults/eventalign.exp   > ${testdir}/nanopolish.txt
		else
			tail -n +2 ${testdir}/eventalign.exp   > ${testdir}/nanopolish.txt
		fi
		tail -n +2 ${testdir}/result.txt  > ${testdir}/f5c.txt
		#todo : this must be fixed for join with two columns ideally
		paste ${testdir}/nanopolish.txt ${testdir}/f5c.txt > ${testdir}/joined_results.txt
		awk -f  scripts/test_eventalign.awk ${testdir}/joined_results.txt > ${testdir}/joined_diff.txt || handle_tests2 "${file}"
	fi


}

mode_test() {
	case $1 in
		valgrind) valgrind --leak-check=full --error-exitcode=1 $cmd > /dev/null || die "valgrind failed" ;;
		gdb) gdb --args "$cmd";;
		cpu) $cmd --disable-cuda=yes > ${testdir}/result.txt; execute_test;;
		cuda) $cmd --disable-cuda=no > ${testdir}/result.txt; execute_test;;
		echo) echo "$cmd > ${testdir}/result.txt";;
		nvprof) nvprof  -f --analysis-metrics -o profile.nvprof "$cmd" --disable-cuda=no --debug-break=5 > /dev/null;;
		custom) shift; $cmd "$@"  > ${testdir}/result.txt; execute_test;;
		*) die "Unknown mode: $1";;
	esac
}

help_msg() {
	echo "Test script for f5c."
	echo "Usage: f5c_dir/script/test_eventalign.sh [-c/-e] [-b bam file] [-g reference genome] [-r fastq/fasta read] mode"
	echo
	echo "mode                 one of: valgrind/gdb/cpu/cuda/echo"
	echo
	echo "-c                   Uses chr22_meth_example test set."
	echo "-e                   Uses rna test set."
	echo "-b [bam file]        Same as f5c -b."
	echo "-K [n]               Same as f5c -K."
	echo "-B [n]               Same as f5c -B."
	echo "-r [read file]       Same as f5c -r."
	echo "-g [ref genome]      Same as f5c -g."
	echo "-t [n]               Number of threads."
	echo "-d                   Download chr22_meth_example data set and exit."
	echo "-h                   Show this help message."
}

#############################################################################################
# parse options and change defaults if necessary
while getopts b:g:r:t:K:B:cdhe opt
do
	case $opt in
		b) bamfile="$OPTARG";;
		g) ref="$OPTARG";;
		r) reads="$OPTARG";;
		t) threads="$OPTARG";;
		c) testdir=test/chr22_meth_example
		   bamfile=${testdir}/reads.sorted.bam
		   ref=${testdir}/humangenome.fa
		   reads=${testdir}/reads.fastq
		   testset_url="https://f5c.page.link/f5c_na12878_test"
		   fallback_url="https://f5c.page.link/f5c_na12878_test_fallback";;
		e) testdir=test/rna
		   bamfile=${testdir}/reads.sorted.bam
		   ref=${testdir}/gencode.v35.transcripts.fa
		   reads=${testdir}/reads.fastq
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

#download the dataset if necessary
download_test_set $testset_url $fallback_url

#validate files
for file in ${bamfile} ${ref} ${reads}; do
	[ -f ${file} ] || die "${file}: File does not exist"
done

#generate the command
cmd="${exepath} eventalign -b ${bamfile} -g ${ref} -r ${reads} -t ${threads} -K $batchsize -B $max_bases --summary ${testdir}/f5c_event_align.summary.txt --secondary=yes --min-mapq=0"
if [ $testdir = test/rna ]; then
	${exepath} index -d ${testdir}/fast5_files ${reads}
	cmd="${cmd}"" --rna"
fi

#run test accordingly
if [ -z "$mode" ]; then
	if [ $testdir = test/chr22_meth_example ]; then
		${exepath} index -t12 --iop 12 -d ${testdir}/fast5_files/ ${reads}
		${cmd} > /dev/null
	else
		${exepath} index -d ${testdir}/fast5_files ${reads}
		${cmd} > ${testdir}/result.txt
	fi
	execute_test
else
	mode_test "$@"
fi
