#!/usr/bin/env bash
# Run f5c on the given input test file

exepath="f5c"

die() {
    local msg="${1}"
    echo "run.sh: ${msg}" >&2
    exit 1
}

usage() {
    echo "Usage: ${0} <test dir>" >&2
    exit 1
}

if [[ "${#}" -eq 0 ]]; then
    usage
elif [[ ! -f "${exepath}" call-methylation ]]; then
    #die "f5c must be compiled first"
    echo nothing > /dev/null
fi

tstdir="${1}"
[[ -d "${tstdir}" ]] || die "${tstdir}: Directory does not exist"

bamdir="${tstdir}/reads.sorted.bam"
fadir="${tstdir}/humangenome.fa"
fastqdir="${tstdir}/reads.fastq"
for file in "${bamdir}" "${fadir}" "${fastqdir}"; do
    [[ -f "${file}" ]] || die "${file}: File does not exist"
done

cd "${tstdir}"
if [[ "${#}" -eq 1 ]]; then
    "${exepath}" call-methylation -b "${bamdir}" -g "${fadir}" -r "${fastqdir}" -p

elif [[ "${#}" -eq 2 ]]; then
    if [[ "${2}" == "valgrind" ]]; then
        valgrind "${exepath}" call-methylation -b "${bamdir}" -g "${fadir}" -r "${fastqdir}"
    elif [[ "${2}" == "gdb" ]]; then
        gdb --args "${exepath}" call-methylation -b "${bamdir}" -g "${fadir}" -r "${fastqdir}"
    else
        usage
    fi
else
    usage
fi

exit 0
