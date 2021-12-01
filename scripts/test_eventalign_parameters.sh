#!/bin/sh

# exit when command fails
set -e

. scripts/common.sh

# defaults
exepath=./f5c
testdir=test/ecoli_2kb_region

bamfile=${testdir}/reads.sorted.bam
ref=${testdir}/draft.fa
reads=${testdir}/reads.fasta
batchsize=256
max_bases=2M
if command -v nproc > /dev/null; then
	threads=$(nproc --all)
else
	threads=8
fi


download_ecoli_2kb_region_big_testresults() {
	mkdir -p test
	tar_path=test/data.tgz
	wget -O $tar_path "https://f5c.page.link/f5c_ecoli_2kb_region_big_testresults" || rm -rf $tar_path ${testdir}_big_testresults
	echo "Extracting. Please wait."
	tar -xf $tar_path -C test || rm -rf $tar_path ${testdir}_big_testresults
	rm -f $tar_path
}


handle_tests() {
	numfailed=$(cat  ${testdir}/joined_diff.txt | awk '{print $NF}' | sort -u -k1,1 | wc -l)
	numcases=$(wc -l < ${testdir}/nanopolish.summary.txt)
	numres=$(wc -l < ${testdir}/f5c.summary.txt)
	echo "$numfailed of $numcases test cases deviated."
	missing=$(echo "$numcases-$numres" | bc)
	echo "$missing entries in the truthset are missing in the testset"
	failp=$(echo "$numfailed*100/$numcases" | bc)
	[ "$failp" -gt 5 ] && die "${1}: Validation failed"
	echo "Validation passed"
}

#args
#1 - filename
#2 - threshold for allowing failures

handle_tests2() {
	numfailed=$(cat  ${testdir}/joined_diff.txt |  wc -l)
	numcases=$(wc -l < ${testdir}/nanopolish.txt)
	numres=$(wc -l < ${testdir}/f5c.txt)
	echo "$numfailed of $numcases test cases deviated."
	missing=$(echo "$numcases-$numres" | bc)
	echo "$missing entries in the truthset are missing in the testset"
	failp=$(echo "$numfailed*100/$numcases" | bc)
	[ "$failp" -gt $2 ] && die "${1}: Validation failed"
	echo "Validation passed"
}


#args
#1 - expected event alignment output
#2 - awk script
#3 - threshold for allowing failures in the  event alignment output
execute_test() {
	echo "----------------comparing summaries-----------------"
	tail -n +2 ${testdir}/eventalign.summary.exp | awk '{print $1"\t"$2"\t"$3"\tdna\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}' > ${testdir}/nanopolish.summary.txt
	tail -n +2 ${testdir}/f5c_event_align.summary.txt | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14}' > ${testdir}/f5c.summary.txt
	join ${testdir}/nanopolish.summary.txt ${testdir}/f5c.summary.txt > ${testdir}/joined_results.txt || echo "Join ran into an issue. Probably just a warning."
	awk -f  scripts/test_eventalign_summary.awk ${testdir}/joined_results.txt > ${testdir}/joined_diff.txt || handle_tests "${testdir}/joined_diff.txt"

	echo "----------------comparing full results--------------"
	tail -n +2 ${testdir}_big_testresults/$1  > ${testdir}/nanopolish.txt
	tail -n +2 ${testdir}/result.txt  > ${testdir}/f5c.txt
	paste ${testdir}/nanopolish.txt ${testdir}/f5c.txt > ${testdir}/joined_results.txt
	awk -f  scripts/$2 ${testdir}/joined_results.txt > ${testdir}/joined_diff.txt || handle_tests2 "${testdir}/joined_diff.txt" $3


}

test -d ${testdir}_big_testresults || download_ecoli_2kb_region_big_testresults

${exepath} index -d ${testdir}/fast5_files ${testdir}/reads.fasta
echo "************************** --print-read-names ****************************"
${exepath} eventalign -b ${bamfile} -g ${ref} -r ${reads} --secondary=yes --min-mapq=0 -B "$max_bases" --summary ${testdir}/f5c_event_align.summary.txt --print-read-names > ${testdir}/result.txt
execute_test "events_print_read_names.exp" "test_eventalign.awk" 5
echo ""
echo "****************************** --scale-events ****************************"
${exepath} eventalign -b ${bamfile} -g ${ref} -r ${reads} --secondary=yes --min-mapq=0 -B "$max_bases" --summary ${testdir}/f5c_event_align.summary.txt --scale-events > ${testdir}/result.txt
execute_test "events_scale_events.exp" "test_eventalign.awk" 5
echo ""
echo "****************************** --sam**************************************"
${exepath} eventalign -b ${bamfile} -g ${ref} -r ${reads} --secondary=yes --min-mapq=0 -B "$max_bases" --summary ${testdir}/f5c_event_align.summary.txt --sam > ${testdir}/result.txt
numfailed=$(diff -y --suppress-common-lines ${testdir}_big_testresults/events.sam.exp ${testdir}/result.txt  | wc -l)
numcases=$(wc -l < ${testdir}_big_testresults/events.sam.exp)
failp=$(echo "$numfailed*100/$numcases" | bc)
echo "$failp of $numcases test cases deviated."
[ "$failp" -gt 5 ] && die "${1}: Validation failed"
echo "Validation passed"
echo ""
echo "**************************** --signal-index********************************"
${exepath} eventalign -b ${bamfile} -g ${ref} -r ${reads} --secondary=yes --min-mapq=0 -B "$max_bases" --summary ${testdir}/f5c_event_align.summary.txt --signal-index > ${testdir}/result.txt
execute_test "event_signal-index.exp" "test_eventalign_signal_index.awk" 5
echo ""
echo "***************************** --samples************************************"
${exepath} eventalign -b ${bamfile} -g ${ref} -r ${reads} --secondary=yes --min-mapq=0 -B "$max_bases" --summary ${testdir}/f5c_event_align.summary.txt --samples > ${testdir}/result.txt
execute_test "event_samples.exp" "test_eventalign_samples.awk" 10
echo ""
echo "***************************** --signal-index --collapse-events************************************"
${exepath} eventalign -b ${bamfile} -g ${ref} -r ${reads} --secondary=yes --min-mapq=0 -B "$max_bases" --summary ${testdir}/f5c_event_align.summary.txt --signal-index --collapse-events > ${testdir}/result.txt
execute_test "collapsed_signal_index.exp" "test_eventalign_signal_index.awk" 10