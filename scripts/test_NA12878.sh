#!/usr/bin/env bash

die() {
    local msg="${1}"
    echo "run.sh: ${msg}" >&2
    exit 1
}

exepath="./f5c"
testdir="test/chr22_meth_example/"

test -d $testdir || die "$testdir does not exist"

bamfile=$testdir/"reads.sorted.bam"
ref=$testdir/"humangenome.fa"
reads=$testdir/"reads.fastq"


for file in "${bamfile}" "${ref}" "${reads}"; do
    [[ -f "${file}" ]] || die "${file}: File does not exist"
done

if [[ "${#}" -eq 0 ]]; then
	echo "${exepath}" -b "${bamfile}" -g "${ref}" -r "${reads}" -t 8  -K256 ">" ${testdir}/result_exact.txt
    "${exepath}" -b "${bamfile}" -g "${ref}" -r "${reads}" -t 8  -K256 > ${testdir}/result_exact.txt
	awk '{print $1,$2,$3,$4,$8,$9,$10}' ${testdir}/result.txt > ${testdir}/result_exact.txt
	awk '{print $1,$2,$3,$4,$8,$9,$10}' ${testdir}/meth.exp > ${testdir}/meth_exact.txt
	diff -q ${testdir}/meth_exact.txt ${testdir}/result_exact.txt  || die "diff ${testdir}/result_exact.txt ${testdir}/meth_exact.txt failed" 
	
	grep -w "chr22" ${testdir}/result.txt | awk '{print $1$2$3$4$8$9$10"\t"$5"\t"$6"\t"$7}' > ${testdir}/result_float.txt
	grep "chr22" ${testdir}/meth.exp | awk '{print $1$2$3$4$8$9$10"\t"$5"\t"$6"\t"$7}'  > ${testdir}/meth_float.txt	
	
	join --nocheck-order ${testdir}/result_float.txt ${testdir}/meth_float.txt	 -a 1 -a 2 | awk -v thresh=0.02 '
		function abs(x){return ((x < 0.0) ? -x : x)} 
		BEGIN{status=0} 
		{
			if(abs($2-$5)>thresh || abs($3-$6)>thresh || abs($4-$7)>thresh)
				{print $0,abs($2-$5),abs($3-$6),abs($4-$7);status=1} 
		} 
		END{if(status>0){exit 1}}
		' > ${testdir}/floatdiff.txt || die "${file}: Validation failed" 	

elif [[ "${#}" -eq 1 ]]; then
    if [[ "${1}" == "valgrind" ]]; then

        valgrind "${exepath}" -b "${bamfile}" -g "${ref}" -r "${reads}" 
    elif [[ "${1}" == "gdb" ]]; then
        gdb --args "${exepath}" -b "${bamfile}" -g "${ref}" -r "${reads}" 
    elif [[ "${1}" == "cpu" ]]; then
        "${exepath}" -b "${bamfile}" -g "${ref}" -r "${reads}" -t 8  --disable-cuda=yes > result.txt	
    elif [[ "${1}" == "cuda" ]]; then
        "${exepath}" -b "${bamfile}" -g "${ref}" -r "${reads}" -t 8  --disable-cuda=no > result.txt		
    elif [[ "${1}" == "echo" ]]; then
		"${exepath}" -b "${bamfile}" -g "${ref}" -r "${reads}" -t 8  ">" result.txt
	else
        echo "wrong option"
		exit 1
    fi
else
    echo "wrong option"
	exit 1
fi

exit 0
