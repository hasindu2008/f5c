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
testset_url="http://genome.cse.unsw.edu.au/tmp/f5c_ecoli_2kb_region_test.tgz"
fallback_url="https://ndownloader.figshare.com/files/13784075?private_link=b04e3976eaed2225b848"

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
	numfailed=$(cat  ${testdir}/joined_diff.txt | awk '{print $NF}' | sort -u -k1,1 | wc -l)
	numcases=$(wc -l < ${testdir}/nanopolish.summary.txt)
	numres=$(wc -l < ${testdir}/f5c.summary.txt)
	echo "$numfailed of $numcases test cases deviated."
	missing=$(echo "$numcases-$numres" | bc)
	echo "$missing entries in the truthset are missing in the testset"
	failp=$(echo "$numfailed/$numcases" | bc)
	[ "$failp" -gt 0 ] && die "${1}: Validation failed"
	echo "Validation passed"
}

handle_tests2() {
	numfailed=$(cat  ${testdir}/joined_diff.txt |  wc -l)
	numcases=$(wc -l < ${testdir}/nanopolish.txt)
	numres=$(wc -l < ${testdir}/f5c.txt)
	echo "$numfailed of $numcases test cases deviated."
	missing=$(echo "$numcases-$numres" | bc)
	echo "$missing entries in the truthset are missing in the testset"
	failp=$(echo "$numfailed/$numcases" | bc)
	[ "$failp" -gt 0 ] && die "${1}: Validation failed"
	echo "Validation passed"
}


execute_test() {
	echo "----------------comparing summaries---------------------------------------------"
	tail -n +2 ${testdir}/eventalign.summary.exp | awk '{print $2"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}' > ${testdir}/nanopolish.summary.txt
	tail -n +2 ${testdir}/f5c_event_align.summary.txt | awk '{print $2"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14}' > ${testdir}/f5c.summary.txt
	join ${testdir}/nanopolish.summary.txt ${testdir}/f5c.summary.txt > ${testdir}/joined_results.txt || echo "Join ran into an issue. Probably just a warning."
	#if [ $testdir = test/chr22_meth_example ]; then
	awk -f  scripts/test_eventalign_summary.awk ${testdir}/joined_results.txt > ${testdir}/joined_diff.txt || handle_tests "${file}"
	#else	
	#	awk -f  scripts/test_eventalign_summary.awk ${testdir}/joined_results.txt || die "${file}: Validation failed"
	#fi
	echo "----------------summaries are good---------------------------------------------"
	echo "----------------comparing full results-------------------------------------------"
	if [ $testdir = test/chr22_meth_example ]; then
		echo "event by event test not implemented not yet implemented!"
	else
		test -d ${testdir}_big_testresults || mkdir ${testdir}_big_testresults/
		test -e ${testdir}_big_testresults/eventalign.exp || wget "http://genome.cse.unsw.edu.au/tmp/f5c_ecoli_2kb_region_test_eventalign.exp.gz" -O ${testdir}_big_testresults/eventalign.exp.gz 
		test -e ${testdir}_big_testresults/eventalign.exp || gunzip ${testdir}_big_testresults/eventalign.exp.gz
		tail -n +2 ${testdir}_big_testresults/eventalign.exp | awk 		'{print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7"\t"$8"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14}'  > ${testdir}/nanopolish.txt
		tail -n +2 ${testdir}/result.txt | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7"\t"$8"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14}' > ${testdir}/f5c.txt
		paste ${testdir}/nanopolish.txt ${testdir}/f5c.txt > ${testdir}/joined_results.txt
		awk -f  scripts/test_eventalign.awk ${testdir}/joined_results.txt > ${testdir}/joined_diff.txt || handle_tests2 "${file}"
	fi


}

mode_test() {
	cmd="${exepath} eventalign -b ${bamfile} -g ${ref} -r ${reads} -t ${threads} -K $batchsize -B $max_bases"
	case $1 in
		valgrind) valgrind $cmd > /dev/null;;
		gdb) gdb --args "$cmd";;
		cpu) $cmd --disable-cuda=yes > ${testdir}/result.txt; execute_test;;
		cuda) $cmd --disable-cuda=no > ${testdir}/result.txt; execute_test;;
		echo) echo "$cmd -t $threads > ${testdir}/result.txt";;
		nvprof) nvprof  -f --analysis-metrics -o profile.nvprof "$cmd" --disable-cuda=no --debug-break=5 > /dev/null;;
		custom) shift; $cmd "$@" --secondary=yes --min-mapq=0 > ${testdir}/result.txt; execute_test;;
		*) die "Unknown mode: $1";;
	esac
}

help_msg() {
	echo "Test script for f5c."
	echo "Usage: f5c_dir/script/test_eventalign.sh [-c] [-b bam file] [-g reference genome] [-r fastq/fasta read] mode"
	echo
	echo "mode                 one of: valgrind/gdb/cpu/cuda/echo"
	echo
	echo "-c                   Uses chr22_meth_example test set."
	echo "-b [bam file]        Same as f5c -b."
	echo "-K [n]               Same as f5c -K."
	echo "-B [n]               Same as f5c -B."
	echo "-r [read file]       Same as f5c -r."
	echo "-g [ref genome]      Same as f5c -g."
	echo "-t [n]               Number of threads."
	echo "-d                   Download chr22_meth_example data set and exit."
	echo "-h                   Show this help message."
}

# parse options
while getopts b:g:r:t:K:B:cdh opt
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
		   testset_url="http://genome.cse.unsw.edu.au/tmp/f5c_na12878_test.tgz"
		   fallback_url="https://ndownloader.figshare.com/files/13784792?private_link=5dd2077f1041412a9518";;
		K) batchsize="$OPTARG";;
		B) max_bases="$OPTARG";;
		d) download_test_set "http://genome.cse.unsw.edu.au/tmp/f5c_na12878_test.tgz" "https://ndownloader.figshare.com/files/13784792?private_link=5dd2077f1041412a9518"
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
for file in ${bamfile} ${ref} ${reads}; do
	[ -f ${file} ] || die "${file}: File does not exist"
done

if [ -z "$mode" ]; then
	if [ $testdir = test/chr22_meth_example ]; then
		${exepath} eventalign -b ${bamfile} -g ${ref} -r ${reads} -t "$threads" -K "$batchsize" -B "$max_bases" --secondary=yes --min-mapq=0 --summary ${testdir}/f5c_event_align.summary.txt > /dev/null
	else
		#test -e ${testdir}/f5c_event_align.summary.txt && rm ${testdir}/f5c_event_align.summary.txt
		${exepath} index -d ${testdir}/fast5_files ${testdir}/reads.fasta
		${exepath} eventalign -b ${bamfile} -g ${ref} -r ${reads} --secondary=yes --min-mapq=0 -B "$max_bases" --summary ${testdir}/f5c_event_align.summary.txt> ${testdir}/result.txt
	fi
		execute_test
else
	mode_test "$@"
fi
