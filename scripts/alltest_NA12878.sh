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
	numfailed=$(cat ${testdir}/floatdiff.txt | wc -l)
	numcases=$(cat ${testdir}/meth_float.txt | wc -l)
	echo "$numfailed of $numcases test cases failed."
	failp=$(echo "$numfailed/$numcases" | bc)
	[ $failp -gt 0 ] && die "${1}: Validation failed"
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
"${exepath}" call-methylation -b "${bamfile}" -g "${ref}" -r "${reads}" -t "${NCPU}"  -K256 -v5 > ${testdir}/result.txt
evaluate

echo "NO IO PROC INTERLEAVE test"
make clean &&  CFLAGS+="-DIO_PROC_NO_INTERLEAVE=1" make
"${exepath}" call-methylation -b "${bamfile}" -g "${ref}" -r "${reads}" -t "${NCPU}"  -K256 -v5 > ${testdir}/result.txt
evaluate



exit 0
