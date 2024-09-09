#!/bin/sh

die() {
	echo "$1"
	exit 1
}

set -e

# download test set given url
#
# $1 - URL
# $2 - Fallback URL (optional)
download_test_set() {
	# data set exists
	if [ -d ${testdir} ]; then
		return
	fi

	mkdir -p test
	tar_path=test/data.tgz
	wget -O $tar_path "$1" || wget -O $tar_path "$2" || rm -rf $tar_path ${testdir}
	echo "Extracting. Please wait."
	tar -xf $tar_path -C test || rm -rf $tar_path ${testdir}
	rm -f $tar_path
}



handle_tests() {
	numfailed=$(wc -l < ${testdir}/diff.txt)
	numcases=$(wc -l < ${ORIG})
	numres=$(wc -l < ${RES})
	echo "$numfailed of $numcases test cases deviated."
	missing=$(echo "$numcases-$numres" | bc)
	echo "$missing entries in the truthset are missing in the testset"
	failp=$(echo "$numfailed*100/$numcases" | bc)
	[ "$failp" -gt ${THRESH} ] && die "${1}: Validation failed"
	echo "Validation passed"
}

execute_test() {

	ORIG=$2
	RES=$1
	THRESH=$3
	diff -y --suppress-common-lines ${ORIG} ${RES} > ${testdir}/diff.txt || handle_tests $testdir

}

testdir=test/hg2_lsk114_reads_1000
download_test_set "https://f5c.bioinf.science/hg2_lsk114_reads_1000"

./f5c eventalign -b ${testdir}/PGXX22394_reads_1000_6.4.2_sup.bam \
		-r ${testdir}/PGXX22394_reads_1000_6.4.2_sup.fastq -g test/chr22_meth_example/humangenome.fa \
		--slow5 ${testdir}/PGXX22394_reads_1000.blow5 -c > ${testdir}/result.txt || die "eventalign failed"
execute_test ${testdir}/result.txt ${testdir}/eventalign.paf 50

./f5c eventalign -b ${testdir}/PGXX22394_reads_1000_6.4.2_sup.bam \
		-r ${testdir}/PGXX22394_reads_1000_6.4.2_sup.fastq -g test/chr22_meth_example/humangenome.fa \
		--slow5 ${testdir}/PGXX22394_reads_1000.blow5 -a > ${testdir}/result.txt || die "eventalign failed"
execute_test ${testdir}/result.txt ${testdir}/eventalign.sam 50


testdir=test/rna
download_test_set "https://f5c.bioinf.science/f5c_rna_test"

./f5c eventalign -b ${testdir}/reads.sorted.bam -g ${testdir}/gencode.v35.transcripts.fa -r ${testdir}//reads.fastq  \
 --slow5 ${testdir}/reads.blow5 --rna -c > ${testdir}/result.txt || die "eventalign failed"
execute_test ${testdir}/result.txt ${testdir}/eventalign.paf 50

./f5c eventalign -b ${testdir}/reads.sorted.bam -g ${testdir}/gencode.v35.transcripts.fa -r ${testdir}//reads.fastq  \
 --slow5 ${testdir}/reads.blow5 --rna -a > ${testdir}/result.txt || die "eventalign failed"
execute_test ${testdir}/result.txt ${testdir}/eventalign.sam 5

testdir=test/uhr_rna004_1k
download_test_set "https://f5c.bioinf.science/uhr_rna004_1k"

./f5c eventalign -b ${testdir}/PNXRXX240011_reads_1k.bam -g ${testdir}/gencode.v40.transcripts.fa -r ${testdir}/PNXRXX240011_reads_1k.fastq  \
 --slow5 ${testdir}/PNXRXX240011_reads_1k.blow5 --rna -c > ${testdir}/result.txt || die "eventalign failed"
execute_test ${testdir}/result.txt ${testdir}/eventalign.paf 50

./f5c eventalign -b ${testdir}/PNXRXX240011_reads_1k.bam -g ${testdir}/gencode.v40.transcripts.fa -r ${testdir}/PNXRXX240011_reads_1k.fastq  \
 --slow5 ${testdir}/PNXRXX240011_reads_1k.blow5  --rna -a > ${testdir}/result.txt || die "eventalign failed"
execute_test ${testdir}/result.txt ${testdir}/eventalign.sam 5

./f5c eventalign -b ${testdir}/PNXRXX240011_reads_1k.bam -g ${testdir}/gencode.v40.transcripts.fa -r ${testdir}/PNXRXX240011_reads_1k.fastq  \
 --slow5 ${testdir}/PNXRXX240011_reads_1k.blow5  --rna --m6anet > ${testdir}/result.txt || die "eventalign failed"
execute_test ${testdir}/result.txt ${testdir}/m6anet.tsv 5