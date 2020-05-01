#!/bin/bash

#Run script from f5c directory. f5c, results, and data folders should all be located in the same directory
#Make sure simd version of f5c is built and running properly

set -e

#default arguments
exedir=../f5c
#where is the raw data for methylation calling located
testdir=''
#what is the name of the folder containing the data and where the results will be stored
dataset=''
#where is the directory containing the dataset directory
resultdir=../results
#number of times to repeat test
num_runs=3
#whether we are running simd or not
simd=1

#f5c parameters
bam_file=''
read_file=''
genome_file=''
K=256 #batchsize
B=2M #max_bases
if command -v nproc > /dev/null; then
	threads=$(nproc --all)
else
	threads=8
fi

#whether we are using a custom dataset or chr22 one
custom_dataset=0

# terminate script
die() {
	echo "$1" >&2
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
	wget -O $tar_path "$1" || wget -O $tar_path "$2" || rm -rf $tar_path ${testdir}
	echo "Extracting. Please wait."
	tar -xf $tar_path -C test || rm -rf $tar_path ${testdir}
	rm -f $tar_path
}

help_msg() {
	echo "Parameter tuning script for f5c."
	echo "Usage: f5c_dir/scripts/tune_parameters.sh args"
	echo
	echo "-T [test data dir]   Directory where test data is located. Must be specified if using custom dataset"
    echo "-s [resultdir]       Directory where test results are located. Default is ../results"
	echo "-d [dataset]   	   Name of dataset to be used. Must be specified if using custom dataset"
    echo "-b [bam file]        Path to the .bam file to be used. Default is {test data dir}/{dataset}/{dataset}.bam"
    echo "-r [read file]       Path to the .fastq file to be used. Default is {test data dir}/{dataset}/{dataset}.fastq"
    echo "-g [ref genome]      Path to the reference genome to be used. Default is {test data dir}/{dataset}/humangenome.fa"
    echo "-n [num_runs]        How many runs to do for each set of parameters. Default is 3"
	echo "-t [threads]         Number of threads to run with"
	echo "-K [batchsize]       Same as f5c -K."
	echo "-B [max_bases]       Same as f5c -B."
	echo "-S 		           Flag to specify that we are NOT using SIMD."
	echo "-h                   Show this help message."
}

# parse options
while getopts T:s:d:b:r:g:n:t:K:B:Sh opt
do
	case $opt in
		t) testdir="$OPTARG";;
        s) resultdir="$OPTARG";;
		d) dataset="$OPTARG";;
        b) bam_file="$OPTARG";;
        r) read_file="$OPTARG";;
        g) genome_file="$OPTARG";;
        n) num_runs=$((${OPTARG}));;
		t) threads="$OPTARG";;
		K) K="$OPTARG";;
		B) B="$OPTARG";;
		S) simd=0;;
		h) help_msg
		   exit 0;;
		?) printf "Usage: %s args\n" "$0"
		   exit 2;;
	esac
done

#set variables if using custom dataset
if [ ! -z ${dataset} ] && [ ! -z ${testdir} ]; then
	if [ -z ${bam_file} ]; then
		bam_file="${testdir}/${dataset}/${dataset}.bam"
	fi
	if [ -z ${read_file} ]; then
		read_file="${testdir}/${dataset}/${dataset}.fastq"
	fi
	if [ -z ${genome_file} ]; then
		genome_file="${testdir}/${dataset}/${dataset}.bam"
	fi
	custom_dataset=1
fi

#result directory configuration
if [ ! -d "${resultdir}" ]; then
	echo "Results directory does not exist. Creating one now..."
	mkdir ${resultdir}
fi

#set f5c run command
if [ ${custom_dataset} -eq 0 ]; then
	#set variables
	testdir=test/chr22_meth_example
	bam_file=${testdir}/reads.sorted.bam
	genome_file=${testdir}/humangenome.fa
	read_file=${testdir}/reads.fastq

	#download data set
	testset_url="http://genome.cse.unsw.edu.au/tmp/f5c_na12878_test.tgz"
	fallback_url="https://ndownloader.figshare.com/files/13784792?private_link=5dd2077f1041412a9518"
	download_test_set $testset_url $fallback_url

	# validate files
	for file in ${bamfile} ${ref} ${reads}; do
		[ -f ${file} ] || die "${file}: File does not exist"
	done
fi

command="./f5c call-methylation -b ${bam_file} -g ${genome_file} -r ${read_file} -t ${threads} -K ${K} -B ${B} --profile-cpu=yes"

if [ ${simd} -eq 1 ]; then
	#overwrite existing files
	-f "${resultdir}/simd_processing.txt" || rm "${resultdir}/simd_processing.txt"
	-f "${resultdir}/simd_align.txt" || rm "${resultdir}/simd_align.txt"

	#SIMD benchmarking
	for i in $( seq 1 ${num_runs} ); do
		echo
		echo "Running SIMD test ${i}"
		echo
		eval ${command} > "${resultdir}/result.txt" 2> "${resultdir}/simd_${i}.txt"

		#grep the output for useful metrics
		cat "${resultdir}/simd_${i}.txt" | grep 'Data processing time:' | grep -o -P '(?<=Data processing time: ).*(?= sec)' >> "${resultdir}/simd_processing.txt"
		cat "${resultdir}/simd_${i}.txt" | grep 'Alignment time:' | grep -o -P '(?<=Alignment time: ).*(?= sec)' >> "${resultdir}/simd_align.txt"
	done
	python3 scripts/average.py "${resultdir}/simd_processing.txt" 'simd' 'processing'
	python3 scripts/average.py "${resultdir}/simd_align.txt" 'simd' 'align'
else
	#overwrite existing files
	-f "${resultdir}/normal_processing.txt" || rm "${resultdir}/normal_processing.txt"
	-f "${resultdir}/normal_align.txt" || rm "${resultdir}/normal_align.txt"

	#Regular benchmarking
	for i in $( seq 1 ${num_runs} ); do
		echo
		echo "Running non-SIMD test ${i}"
		echo
		eval ${command} > "${resultdir}/result.txt" 2> "${resultdir}/normal_${i}.txt"

		#grep the output for useful metrics
		cat "${resultdir}/normal_${i}.txt" | grep 'Data processing time:' | grep -o -P '(?<=Data processing time: ).*(?= sec)' >> "${resultdir}/normal_processing.txt"
		cat "${resultdir}/normal_${i}.txt" | grep 'Alignment time:' | grep -o -P '(?<=Alignment time: ).*(?= sec)' >> "${resultdir}/normal_align.txt"
	done
	python3 scripts/average.py "${resultdir}/normal_processing.txt" 'normal' 'processing'
	python3 scripts/average.py "${resultdir}/normal_align.txt" 'normal' 'align'
fi
