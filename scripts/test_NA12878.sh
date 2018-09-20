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
	echo "${exepath}" -b "${bamfile}" -g "${ref}" -r "${reads}" -t 8 --secondary=yes --min-mapq=0 --print-scaling=yes ">" result.txt
    "${exepath}" -b "${bamfile}" -g "${ref}" -r "${reads}" -t 8 --secondary=yes --min-mapq=0 --print-scaling=yes > result.txt

elif [[ "${#}" -eq 1 ]]; then
    if [[ "${1}" == "valgrind" ]]; then
        valgrind "${exepath}" -b "${bamfile}" -g "${ref}" -r "${reads}" --secondary=yes --min-mapq=0
    elif [[ "${1}" == "gdb" ]]; then
        gdb --args "${exepath}" -b "${bamfile}" -g "${ref}" -r "${reads}" --secondary=yes --min-mapq=0
    elif [[ "${1}" == "cpu" ]]; then
        "${exepath}" -b "${bamfile}" -g "${ref}" -r "${reads}" -t 8 --secondary=yes --min-mapq=0 --print-scaling=yes --disable-cuda=yes > result.txt	
    elif [[ "${1}" == "cuda" ]]; then
        "${exepath}" -b "${bamfile}" -g "${ref}" -r "${reads}" -t 8 --secondary=yes --min-mapq=0 --print-scaling=yes --disable-cuda=no > result.txt		
    elif [[ "${1}" == "echo" ]]; then
		"${exepath}" -b "${bamfile}" -g "${ref}" -r "${reads}" -t 8 --secondary=yes --min-mapq=0 --print-scaling=yes ">" result.txt
	else
        echo "wrong option"
		exit 1
    fi
else
    echo "wrong option"
	exit 1
fi

exit 0
