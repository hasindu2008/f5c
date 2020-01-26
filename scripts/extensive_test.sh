#!/usr/bin/env bash

die() {
	local msg="${1}"
	echo "extensive_test: ${msg}" >&2
	exit 1
}

set_cpu_count() {
    if command -v nproc > /dev/null; then
        NCPU="$(nproc --all)"
    else
        NCPU=8
    fi
}

handle_tests(){
	numfailed=$(wc -l < "${testdir}"/floatdiff.txt)
	numcases=$(wc -l < "${testdir}"/meth_float.txt)
	echo "$numfailed of $numcases test cases failed."
	failp=$(echo "$numfailed*100/$numcases" | bc)
	[ "$failp" -gt 5 ] && die "${1}: Validation failed"
}


exepath="./f5c"
testdir="test/chr22_meth_example/"

test -d $testdir || die "$testdir does not exist"

bamfile=$testdir/"reads.sorted.bam"
ref=$testdir/"humangenome.fa"
reads=$testdir/"reads.fastq"

evaluate(){
	grep -w "chr20" ${testdir}/result.txt | awk '{print $1$2$3$4$8$9$10"\t"$5"\t"$6"\t"$7}' > ${testdir}/result_float.txt
	grep -w "chr20" ${testdir}/meth.exp | awk '{print $1$2$3$4$8$9$10"\t"$5"\t"$6"\t"$7}'  > ${testdir}/meth_float.txt

	join --nocheck-order ${testdir}/result_float.txt ${testdir}/meth_float.txt	 -a 1 -a 2 | awk -v thresh=0.1 '
	function abs(x){return (((x) < 0.0) ? -(x) : (x))}
	BEGIN{status=0}
	{
		if(abs($2-$5)>abs(thresh*$5)+0.02 || abs($3-$6)>abs(thresh*$6)+0.02 || abs($4-$7)>abs(thresh*$7)+0.02)
			{print $0,abs($2-$5)">"abs(thresh*$5)+0.02 ,abs($3-$6)">"abs(thresh*$6)+0.02,abs($4-$7)">"abs(thresh*$7)+0.02;status=1}
		}
	END{if(status>0){exit 1}}
	' > ${testdir}/floatdiff.txt || handle_tests "${file}: Validation failed"

}


set_cpu_count
for file in "${bamfile}" "${ref}" "${reads}"; do
	[[ -f "${file}" ]] || die "${file}: File does not exist"
done

test_suit1 () {

	make clean
	make -j8
	echo "-----------------Doing ecoli based CPU tests-------------------------"
	echo "Methylation calling"
	scripts/test.sh 2> ecoli_methcalling.log || die "failed"
	echo "____________________________________________________________________"
	echo "event alignment"
	scripts/test_eventalign.sh 2> ecoli_eventalign.log || die "failed"
	echo "____________________________________________________________________"	
	echo "methylation frequency"
	scripts/test_methfreq.sh 2> ecoli_methfreq.log || die "failed"
	echo "____________________________________________________________________"	
	echo "multi-fast5"
	scripts/test_multifast5.sh 2> ecoli_multifast5.log  || die "failed"
	echo "____________________________________________________________________"	

	echo ""
	echo "--------------------------------------------------------------------"
	echo ""

	echo "----------------Doing NA12878 based CPU tests-----------------------"
	echo "Methylation calling"
	scripts/test.sh -c 2> na12878_methcalling.log || die "failed"
	echo "____________________________________________________________________"
	echo "event alignment"
	scripts/test_eventalign.sh -c 2> na12878_eventalign.log || die "failed"
	echo "____________________________________________________________________"
	echo "methylation frequency"
	scripts/test_methfreq.sh -c 2> na12878_methfreq.log || die "failed"
	echo "____________________________________________________________________"
	echo "multi-fast5"
	scripts/test_multifast5.sh -c 2> na12878_multifast5.log || die "failed"
	echo "____________________________________________________________________"

}


test_suit1_cuda () {

	make clean
	make cuda=1 -j8

	echo "---------------Doing ecoli based CUDA tests-------------------------"
	echo "Methylation calling"
	scripts/test.sh 2> ecoli_methcalling_cuda.log || die "failed"
	echo "____________________________________________________________________"
	echo "event alignment"
	scripts/test_eventalign.sh 2> ecoli_eventalign_cuda.log || die "failed"
	echo "____________________________________________________________________"
	echo "methylation frequency"
	scripts/test_methfreq.sh 2> ecoli_methfreq_cuda.log || die "failed"
	echo "____________________________________________________________________"
	echo "multi-fast5"
	scripts/test_multifast5.sh 2> ecoli_multifast5_cuda.log  || die "failed"
	echo "____________________________________________________________________"

	echo ""
	echo "--------------------------------------------------------------------"
	echo ""

	echo "---------------Doing NA12878 based CUDA tests-----------------------"
	echo "Methylation calling"
	scripts/test.sh -c 2> na12878_methcalling_cuda.log || die "failed"
	echo "____________________________________________________________________"
	echo "event alignment"
	scripts/test_eventalign.sh -c 2> na12878_eventalign_cuda.log || die "failed"
	echo "____________________________________________________________________"
	echo "methylation frequency"
	scripts/test_methfreq.sh -c 2> na12878_methfreq_cuda.log || die "failed"
	echo "____________________________________________________________________"
	echo "multi-fast5"
	scripts/test_multifast5.sh -c 2> na12878_multifast5_cuda.log || die "failed"
	echo "____________________________________________________________________"

}


test_suit2 () {

	echo "---------------Doing NA12878 based CPU tests part 2-----------------------"
	echo "Default test"
	make clean && make
	"${exepath}" call-methylation -b "${bamfile}" -g "${ref}" -r "${reads}" -t "${NCPU}"  -K1024 -v5 > ${testdir}/result.txt 2> default.log
	evaluate
	echo ""
	echo "____________________________________________________________________"

	echo "sectional benchmark"
	"${exepath}" call-methylation -b "${bamfile}" -g "${ref}" -r "${reads}" -t "${NCPU}"  -K1024 -v5 --profile=yes > ${testdir}/result.txt 2> profile.log
	evaluate
	echo ""
	echo "____________________________________________________________________"

	echo "NO IO PROC INTERLEAVE test"
	make clean &&  CFLAGS+="-DIO_PROC_NO_INTERLEAVE=1" make
	"${exepath}" call-methylation -b "${bamfile}" -g "${ref}" -r "${reads}" -t "${NCPU}"  -K1024 -v5 > ${testdir}/result.txt 2> no_io_proc.log
	evaluate
	echo ""
	echo "____________________________________________________________________"


	echo "bad fast5 file"
	mv test/chr22_meth_example/fast5_files/DEAMERNANOPORE_20161117_FNFAB43577_MN16450_sequencing_run_MA_821_R9_4_NA12878_11_17_16_88738_ch1_read445_strand.fast5 test/chr22_meth_example/fast5_files/a.fast5
	make clean && make
	"${exepath}" call-methylation -b "${bamfile}" -g "${ref}" -r "${reads}" -t "${NCPU}"  -K1024 -v5 > ${testdir}/result.txt 2> badfast5.log
	mv test/chr22_meth_example/fast5_files/a.fast5 test/chr22_meth_example/fast5_files/DEAMERNANOPORE_20161117_FNFAB43577_MN16450_sequencing_run_MA_821_R9_4_NA12878_11_17_16_88738_ch1_read445_strand.fast5
	evaluate
	echo ""
	echo "____________________________________________________________________"

	echo "IOP test : I/O processes"
	make clean && make
	"${exepath}" call-methylation -b "${bamfile}" -g "${ref}" -r "${reads}" -t "${NCPU}"  -K256 -v5 --iop8 > ${testdir}/result.txt 2> cuda_malloc.log
	evaluate
	echo ""
	echo "____________________________________________________________________"

}

test_suit2_cuda () {

	echo "---------------Doing NA12878 based CUDA tests part 2 -----------------------"
	echo "CUDA test"
	make clean && make cuda=1
	"${exepath}" call-methylation -b "${bamfile}" -g "${ref}" -r "${reads}" -t "${NCPU}"  -K256 -v5 > ${testdir}/result.txt 2> default_cuda.log
	evaluate
	echo ""
	echo "____________________________________________________________________"


	echo "CUDA test : cuda disabled"
	make clean && make cuda=1
	"${exepath}" call-methylation -b "${bamfile}" -g "${ref}" -r "${reads}" -t "${NCPU}"  -K256 -v5 --disable-cuda=yes > ${testdir}/result.txt 2> cuda_disabled.log
	evaluate
	echo ""
	echo "____________________________________________________________________"

	echo "CUDA test : dynamic malloc"
	make clean && CFLAGS_CUDA+="-DCUDA_DYNAMIC_MALLOC=1" make cuda=1
	"${exepath}" call-methylation -b "${bamfile}" -g "${ref}" -r "${reads}" -t "${NCPU}"  -K256 -v5 > ${testdir}/result.txt 2> cuda_malloc.log
	evaluate
	echo ""
	echo "____________________________________________________________________"

}




help_msg() {
	echo "Extensive test script for f5c."
	echo "Usage: f5c_dir/script/extensive_test.sh mode"
	echo
	echo "mode                 one of: cpu/gpu/all"
	echo
	echo "-h                   Show this help message."
}

# parse options
while getopts h opt
do
	case $opt in
		h) help_msg
		   exit 0;;
		?) printf "Usage: %s args" "$0"
		   exit 2;;
	esac
done
shift $(($OPTIND - 1))
mode=$1

if [ "$mode" = "cpu" -o  "$mode" = "all" ]; then
	test_suit1
	test_suit2
fi
if [ "$mode" = "gpu" -o  "$mode" = "all" ]; then	
	test_suit1_cuda
	test_suit2_cuda
fi
if [ "$mode" = "cpu" -o "$mode" != "gpu" -o "$mode" != "all" ]; then	
	help_msg
	exit 2;
fi

exit 0
