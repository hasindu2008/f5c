#!/bin/sh

# exit when command fails
set -e

# defaults
exepath=./f5c
testdir=test/ecoli_2kb_region

bamfile=${testdir}/reads.sorted.bam
ref=${testdir}/draft.fa
reads=${testdir}/reads.fasta
batchsize=256
if command -v nproc > /dev/null; then
	threads=$(nproc --all)
else
	threads=8
fi
# execution mode (valgrind/gdb/cpu/cuda/echo)
mode=
testset_url="http://genome.cse.unsw.edu.au/tmp/f5c_ecoli_2kb_region_test.tgz"
fallback_url="https://ndownloader.figshare.com/files/13784075?private_link=b04e3976eaed2225b848"

# terminate script
die() {
	echo "test.sh: $1" >&2
	exit 1
}

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
	wget -O $tar_path $1 || wget -O $tar_path $2 || rm -rf $tar_path ${testdir}
	echo "Extracting. Please wait."
	tar -xf $tar_path -C test || rm -rf $tar_path ${testdir}
	rm -f $tar_path
}

mode_test() {
	cmd="${exepath} call-methylation -b ${bamfile} -g ${ref} -r ${reads}"

	case $1 in
		valgrind) valgrind $cmd;;
		gdb) gdb --args $cmd;;
		cpu) $cmd --disable-cuda=yes > result.txt;;
		cuda) $cmd --disable-cuda=no > result.txt;;
		echo) echo "$cmd -t $threads > result.txt";;
		*) die "Unknown mode: $1";;
	esac
}

handle_tests() {
	numfailed=$(cat ${testdir}/floatdiff.txt | wc -l)
	numcases=$(cat ${testdir}/meth_float.txt | wc -l)
	echo "$numfailed of $numcases test cases failed."
	failp=$(echo "$numfailed/$numcases" | bc)
	[ $failp -gt 0 ] && die "${1}: Validation failed"
}

help_msg() {
	echo "Test script for f5c."
	echo "Usage: f5c_dir/script/test.sh [-c] [-b bam file] [-g reference genome] [-r fastq/fasta read] mode"
	echo
	echo "mode                 one of: valgrind/gdb/cpu/cuda/echo"
	echo
	echo "-c                   Uses chr22_meth_example test set."
	echo "-b [bam file]        Same as f5c -b."
	echo "-K [n]               Same as f5c -K."
	echo "-r [read file]       Same as f5c -r."
	echo "-g [ref genome]      Same as f5c -g."
	echo "-t [n]               Number of threads."
	echo "-d                   Download chr22_meth_example data set and exit."
	echo "-h                   Show this help message."
}

# parse options
while getopts b:g:r:t:K:cdh opt
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
		d) download_test_set "http://genome.cse.unsw.edu.au/tmp/f5c_na12878_test.tgz" "https://ndownloader.figshare.com/files/13784792?private_link=5dd2077f1041412a9518"
		   exit 0;;
		h) help_msg
		   exit 0;;
		?) printf "Usage: %s [-c] [-b bam file] [-g reference genome] [-r fastq/fasta read] args" $0
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

if [ -z $mode ]; then
	if [ $testdir = test/chr22_meth_example ]; then
		${exepath} call-methylation -b ${bamfile} -g ${ref} -r ${reads} -t ${threads} -K $batchsize > ${testdir}/result.txt
		grep -w "chr20" ${testdir}/result.txt | awk '{print $1$2$3$4$8$9$10"\t"$5"\t"$6"\t"$7}' > ${testdir}/result_float.txt
		grep -w "chr20" ${testdir}/meth.exp | awk '{print $1$2$3$4$8$9$10"\t"$5"\t"$6"\t"$7}'  > ${testdir}/meth_float.txt

		join -a 1 -a 2 ${testdir}/result_float.txt ${testdir}/meth_float.txt | awk -v thresh=0.1 '
		function abs(x){return (((x) < 0.0) ? -(x) : (x))}
		BEGIN{status=0}
		{
			if(abs($2-$5)>abs(thresh*$5)+0.02 || abs($3-$6)>abs(thresh*$6)+0.02 || abs($4-$7)>abs(thresh*$7)+0.02)
				{print $0,abs($2-$5)">"abs(thresh*$5)+0.02 ,abs($3-$6)">"abs(thresh*$6)+0.02,abs($4-$7)">"abs(thresh*$7)+0.02;status=1}
			}
		END{if(status>0){exit 1}}
		' > ${testdir}/floatdiff.txt || handle_tests "${file}: Validation failed"
	else
		${exepath} index -d ${testdir}/fast5_files ${testdir}/reads.fasta
		${exepath} call-methylation -b ${bamfile} -g ${ref} -r ${reads} --secondary=yes --min-mapq=0 > ${testdir}/result.txt
		awk '{print $1,$2,$3,$4,$8,$9,$10}' ${testdir}/result.txt > ${testdir}/result_exact.txt
		awk '{print $1,$2,$3,$4,$8,$9,$10}' ${testdir}/meth.exp > ${testdir}/meth_exact.txt
		diff -q ${testdir}/meth_exact.txt ${testdir}/result_exact.txt || die "diff ${testdir}/result_exact.txt ${testdir}/meth_exact.txt failed" 

		awk '{print $1$2$3$4$8$9$10"\t"$5"\t"$6"\t"$7}' ${testdir}/result.txt > ${testdir}/result_float.txt
		awk '{print $1$2$3$4$8$9$10"\t"$5"\t"$6"\t"$7}' ${testdir}/meth.exp > ${testdir}/meth_float.txt	

		join -a 1 -a 2 ${testdir}/result_float.txt ${testdir}/meth_float.txt | awk -v thresh=0.1 '
		function abs(x){return (((x) < 0.0) ? -(x) : (x))} 
		BEGIN{status=0} 
		{
			if(abs($2-$5)>abs(thresh*$5)+0.02 || abs($3-$6)>abs(thresh*$6)+0.02 || abs($4-$7)>abs(thresh*$7)+0.02)
				{print $0,abs($2-$5)">"abs(thresh*$5)+0.02 ,abs($3-$6)">"abs(thresh*$6)+0.02,abs($4-$7)">"abs(thresh*$7)+0.02;status=1} 
			} 
		END{if(status>0){exit 1}}
		' > ${testdir}/floatdiff.txt || die "${file}: Validation failed" 
	fi
else
	mode_test $mode
fi
