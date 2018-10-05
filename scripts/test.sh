#!/usr/bin/env bash

die() {
    local msg="${1}"
    echo "run.sh: ${msg}" >&2
    exit 1
}

exepath="./f5c"
testdir="test/ecoli_2kb_region/"

bamfile="${testdir}reads.sorted.bam"
ref="${testdir}draft.fa"
reads="${testdir}reads.fasta"


for file in "${bamfile}" "${ref}" "${reads}"; do
    [[ -f "${file}" ]] || die "${file}: File does not exist"
done

if [[ "${#}" -eq 0 ]]; then
    "${exepath}" -b "${bamfile}" -g "${ref}" -r "${reads}" --secondary=yes --min-mapq=0 --print-scaling=yes > ${testdir}/result.txt
	diff -q ${testdir}/testresult.exp ${testdir}/result.txt  || echo "${file}: Diff failed" 
	join --nocheck-order ${testdir}/testresult.exp ${testdir}/result.txt -a 1 -a 2 | awk 'function abs(x){return ((x < 0.0) ? -x : x)} BEGIN{status=0} {if(abs($2-$5)>0.01 || abs($3-$6)>0.01 || abs($4-$7)>0.01){print $0;status=1} } END{if(status>0){exit 1}}' || die "${file}: Validation failed" 
elif [[ "${#}" -eq 1 ]]; then
    if [[ "${1}" == "valgrind" ]]; then
        valgrind "${exepath}" -b "${bamfile}" -g "${ref}" -r "${reads}" --secondary=yes --min-mapq=0
    elif [[ "${1}" == "gdb" ]]; then
        gdb --args "${exepath}" -b "${bamfile}" -g "${ref}" -r "${reads}" --secondary=yes --min-mapq=0
    elif [[ "${1}" == "echo" ]]; then 
		echo "${exepath}" -b "${bamfile}" -g "${ref}" -r "${reads}" --secondary=yes --min-mapq=0 --print-scaling=yes ">" ${testdir}/result.txt
	else
        echo "wrong option"
		exit 1
    fi
else
    echo "wrong option"
	exit 1
fi

exit 0
