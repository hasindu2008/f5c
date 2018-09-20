#!/usr/bin/env bash

die() {
    local msg="${1}"
    echo "run.sh: ${msg}" >&2
    exit 1
}

exepath="./f5c"
testdir="test/ecoli_2kb_region/"
#cd "${tstdir}"

bamdir="${testdir}single_read/read1.sorted.bam"
fadir="${testdir}draft.fa"
fastqdir="${testdir}single_read/read1.fasta"
for file in "${bamdir}" "${fadir}" "${fastqdir}"; do
    [[ -f "${file}" ]] || die "${file}: File does not exist"
done

if [[ "${#}" -eq 0 ]]; then
    "${exepath}" -b "${bamdir}" -g "${fadir}" -r "${fastqdir}" --secondary=yes --min-mapq=20  > ${testdir}/result.txt

elif [[ "${#}" -eq 1 ]]; then
    if [[ "${1}" == "valgrind" ]]; then
        valgrind "${exepath}" -b "${bamdir}" -g "${fadir}" -r "${fastqdir}"
    elif [[ "${1}" == "gdb" ]]; then
        gdb --args "${exepath}" -b "${bamdir}" -g "${fadir}" -r "${fastqdir}"
	elif [[ "${1}" == "echo" ]]; then	
		echo "${exepath}" -b "${bamdir}" -g "${fadir}" -r "${fastqdir}" --secondary=yes --min-mapq=20  ">" ${testdir}/result.txt
    else
        echo "wrong option"
		exit 1
    fi
else
    echo "wrong option"
	exit 1
fi

exit 0
