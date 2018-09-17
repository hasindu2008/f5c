#!/usr/bin/env bash

die() {
    local msg="${1}"
    echo "run.sh: ${msg}" >&2
    exit 1
}

exepath="../nanopolish/nanopolish"
testdir="test/ecoli_2kb_region/"
#cd "${tstdir}"

bamdir="${testdir}single_read/read1.sorted.bam"
fadir="${testdir}draft.fa"
fastqdir="${testdir}single_read/read1.fasta"
for file in "${bamdir}" "${fadir}" "${fastqdir}"; do
    [[ -f "${file}" ]] || die "${file}: File does not exist"
done

if [[ "${#}" -eq 0 ]]; then
    #"${exepath}" -b "${bamdir}" -g "${fadir}" -r "${fastqdir}" 
    "${exepath}" call-methylation -t 1 -r "$fastqdir" -b "${bamdir}" -g "${fadir} --secondary=yes --min-mapq=0 "

elif [[ "${#}" -eq 1 ]]; then
    if [[ "${1}" == "valgrind" ]]; then
        valgrind "${exepath}" -b "${bamdir}" -g "${fadir}" -r "${fastqdir}"
    elif [[ "${1}" == "gdb" ]]; then
        gdb --args "${exepath}" -b "${bamdir}" -g "${fadir}" -r "${fastqdir}"
    else
        echo "wrong option"
		exit 1
    fi
else
    echo "wrong option"
	exit 1
fi

exit 0
