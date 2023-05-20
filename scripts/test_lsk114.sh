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

testdir=test/hg2_lsk114_reads_1000

download_test_set "https://f5c.page.link/hg2_lsk114_reads_1000"

./f5c resquiggle ${testdir}/PGXX22394_reads_1000_6.4.2_sup.fastq ${testdir}/PGXX22394_reads_1000.blow5 > ${testdir}/rsq.txt || die "resquiggle failed"

./f5c eventalign -b ${testdir}/PGXX22394_reads_1000_6.4.2_sup.bam \
		-r ${testdir}/PGXX22394_reads_1000_6.4.2_sup.fastq -g test/chr22_meth_example/humangenome.fa \
		--slow5 ${testdir}/PGXX22394_reads_1000.blow5 > ${testdir}/result.txt || die "eventalign failed"

./f5c call-methylation -b ${testdir}/PGXX22394_reads_1000_6.4.2_sup.bam \
		-r ${testdir}/PGXX22394_reads_1000_6.4.2_sup.fastq -g test/chr22_meth_example/humangenome.fa \
		--slow5 ${testdir}/PGXX22394_reads_1000.blow5 > ${testdir}/meth.tsv || die "eventalign failed"

./f5c eventalign -b ${testdir}/PGXX22394_reads_1000_6.4.2_sup.bam \
		-r ${testdir}/PGXX22394_reads_1000_6.4.2_sup.fastq -g test/chr22_meth_example/humangenome.fa \
		--slow5 ${testdir}/PGXX22394_reads_1000.blow5 -c > ${testdir}/result.txt || die "eventalign failed"

diff -q ${testdir}/result.txt ${testdir}/eventalign.paf || die "eventalign failed"