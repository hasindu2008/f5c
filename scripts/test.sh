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
    "${exepath}" -b "${bamfile}" -g "${ref}" -r "${reads}" --secondary=yes --min-mapq=0 --print-scaling=yes> ${testdir}/result.txt
	diff -q ${testdir}/testresult.exp ${testdir}/result.txt  || die "${file}: Diff failed" 
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
