#!/usr/bin/env bash
# Run f5c on the given input test file

exepath="../../f5c"
tstdir="test/ecoli_2kb_region/"
cd "${tstdir}"

bamdir="reads.sorted.bam"
fadir="draft.fa"
fastqdir="reads.fasta"
for file in "${bamdir}" "${fadir}" "${fastqdir}"; do
    [[ -f "${file}" ]] || die "${file}: File does not exist"
done

if [[ "${#}" -eq 0 ]]; then
    "${exepath}" -b "${bamdir}" -g "${fadir}" -r "${fastqdir}" -p

elif [[ "${#}" -eq 1 ]]; then
    if [[ "${1}" == "valgrind" ]]; then
        valgrind "${exepath}" -b "${bamdir}" -g "${fadir}" -r "${fastqdir}"
    elif [[ "${1}" == "gdb" ]]; then
        gdb --args "${exepath}" -b "${bamdir}" -g "${fadir}" -r "${fastqdir}"
    else
        echo "wrong option"
    fi
else
    echo ""
fi

exit 0
