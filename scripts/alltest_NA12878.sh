#!/usr/bin/env bash

die() {
	local msg="${1}"
	echo "run.sh: ${msg}" >&2
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
	failp=$(echo "$numfailed/$numcases" | bc)
	[ "$failp" -gt 0 ] && die "${1}: Validation failed"
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

echo "Default test"
make clean && make
"${exepath}" call-methylation -b "${bamfile}" -g "${ref}" -r "${reads}" -t "${NCPU}"  -K1024 -v5 > ${testdir}/result.txt 2> default.log
evaluate

echo "sectional benchmark"
"${exepath}" call-methylation -b "${bamfile}" -g "${ref}" -r "${reads}" -t "${NCPU}"  -K1024 -v5 --profile=yes > ${testdir}/result.txt 2> profile.log
evaluate

echo "NO IO PROC INTERLEAVE test"
make clean &&  CFLAGS+="-DIO_PROC_NO_INTERLEAVE=1" make
"${exepath}" call-methylation -b "${bamfile}" -g "${ref}" -r "${reads}" -t "${NCPU}"  -K1024 -v5 > ${testdir}/result.txt 2> no_io_proc.log
evaluate

echo "CUDA test"
make clean && make cuda=1
"${exepath}" call-methylation -b "${bamfile}" -g "${ref}" -r "${reads}" -t "${NCPU}"  -K256 -v5 > ${testdir}/result.txt 2> default_cuda.log
evaluate

echo "CUDA test : cuda disabled"
make clean && make cuda=1
"${exepath}" call-methylation -b "${bamfile}" -g "${ref}" -r "${reads}" -t "${NCPU}"  -K256 -v5 --disable-cuda=yes > ${testdir}/result.txt 2> cuda_disabled.log
evaluate


echo "CUDA test : dynamic malloc"
make clean && CFLAGS_CUDA+="-DCUDA_DYNAMIC_MALLOC=1" make cuda=1
"${exepath}" call-methylation -b "${bamfile}" -g "${ref}" -r "${reads}" -t "${NCPU}"  -K256 -v5 > ${testdir}/result.txt 2> cuda_malloc.log
evaluate

echo "bad fast5 file"
mv test/chr22_meth_example/fast5_files/DEAMERNANOPORE_20161117_FNFAB43577_MN16450_sequencing_run_MA_821_R9_4_NA12878_11_17_16_88738_ch1_read445_strand.fast5 test/chr22_meth_example/fast5_files/a.fast5
make clean && make
"${exepath}" call-methylation -b "${bamfile}" -g "${ref}" -r "${reads}" -t "${NCPU}"  -K1024 -v5 > ${testdir}/result.txt 2> badfast5.log
mv test/chr22_meth_example/fast5_files/a.fast5 test/chr22_meth_example/fast5_files/DEAMERNANOPORE_20161117_FNFAB43577_MN16450_sequencing_run_MA_821_R9_4_NA12878_11_17_16_88738_ch1_read445_strand.fast5
evaluate

exit 0
