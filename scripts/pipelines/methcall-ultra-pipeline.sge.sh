#!/bin/bash
#$ -cwd
#$ -V
#$ -N F5C-METHCALL
#$ -S /bin/bash
#$ -b y
#$ -pe smp 32
#$ -l mem_requested=4G
#$ -l h_vmem=4G
#$ -l tmp_requested=150G
#$ -l nvgpu=1

###################################################################

## example SGE script - resource efficient methylation calling pipeline for a dataset with many ultra-long reads
## tested on the hpc cluster at Garvan Institute of Medical Research (on a GPU node with Tesla V100 - 16GB GPUs)
## 32 CPU cores have been requested using -pe smp 32
## 1 GPU has been requested using -l nvgpu=1
## note that mem_requested=4G, h_vmem=4G, tmp_requested=150G values are for a single core
## thus the total memory requested is 4GBx32=128GB and total temporary local scratch space requested is 150GBx32=4800GB
## it is assumed that the compute node has adequate local scratch storage space and is accessible via $TMPDIR
## the local scratch space should be adequate to hold the reference genome, FASTQ file, all the fast5 files and the intermediate files: generated f5c index files (~0.6 of FASTQ size), SAM file, BAM file, methylation calls (tsv files), etc.

#change this to where the fast5 files are located on the shared cluster storage (relative or absolute path)
#this directory will be later copied to the local scratch storage
FAST5DIR=fast5/

#change this to where the FASTQ file is located on the shared cluster storage (relative or absolute path)
#this file will be later copied to the local scratch storage
FASTQFILE=na12878_rel6_all.fastq

#change this to the directory where you want the output to be copied on the shared cluster storage (relative or absolute path)
#the directory will be created if not already existing
#the SAM file (reads.sam), BAM file (reads.bam, tmp.bam), per-read methylation calls (meth.tsv, meth-ultra.tsv), methylation frequencies (meth-freq.tsv, meth-ultra-freq.tsv, meth-freq-combined.tsv) will be first written to the local scratch storage and then will be copied to here at the end of each pipeline step.
#content may be overwritten if they already exists
#the location should be writable and should have adequate free space
OUTPUT_DIR=./output

#the reference genome
#this file will be later copied to the local scratch storage
REFERENCE=hg38noAlt.fa

#command names or paths. note that you may have to load appropriate modules prior to these.
MINIMAP=minimap2
SAMTOOLS=samtools
F5C=f5c

#number of threads to launch - should match with -pe smp
num_threads=32

###################################################################
###################################################################

# terminate script
die() {
	echo "$1" >&2
	echo
	exit 1
}

#check if the files/directories exist
test -d $FAST5DIR || die "directory $FAST5DIR not found"
test -e $FASTQFILE || die "input file $FASTQFILE not found"
test -e $REFERENCE || die "reference genome $REFERENCE not found"
test -d $OUTPUT_DIR || mkdir $OUTPUT_DIR || die "Could not create directory $OUTPUT_DIR"

#temporary location of the inputs and outputs on local scratch storage. 
REF_LOCAL=$TMPDIR/ref.fa
FASTQ_LOCAL=$TMPDIR/reads.fastq
SAM_LOCAL=$TMPDIR/reads.sam
BAM_LOCAL=$TMPDIR/reads.bam
METH_LOCAL=$TMPDIR/meth.tsv
METH_ULTRA_LOCAL=$TMPDIR/meth-ultra.tsv
METH_FREQ_LOCAL=$TMPDIR/meth-freq.tsv
METH_ULTRA_FREQ_LOCAL=$TMPDIR/meth-ultra-freq.tsv
TMP_BAM_LOCAL=$TMPDIR/tmp.bam
METH_FREQ_COMBINED_LOCAL=$TMPDIR/meth-freq-combined.tsv

mkdir $TMPDIR/fast5  || die "Could not create temporary directories at $TMPDIR" 

#copy the reference and the FASTQ file to the local scratch storage
echo "Copying reference $REFERENCE to $REF_LOCAL"
cp $REFERENCE  $REF_LOCAL || die "Copying $REFERENCE to $REF_LOCAL failed"
echo "Copying fastq file"
cp $FASTQFILE  $FASTQ_LOCAL || die "Copying $FASTQFILE to $FASTQ_LOCAL failed"	

#minimap
echo "Minimap2 running"
/usr/bin/time -v $MINIMAP -x map-ont -a -t$num_threads --secondary=no $REF_LOCAL $FASTQ_LOCAL > $SAM_LOCAL  || die "minimap2 failed"
echo "copying minimap2 results data back"
cp $SAM_LOCAL $OUTPUT_DIR/ || die "copying minimap2 results data back failed"

#sorting and indexing the BAM
echo "samtools sorting"
/usr/bin/time -v $SAMTOOLS sort -@$num_threads $SAM_LOCAL > $BAM_LOCAL || die "samtools failed"
/usr/bin/time -v $SAMTOOLS index $BAM_LOCAL || die "samtools index failed"
echo "copying samtools results back"
cp $BAM_LOCAL $OUTPUT_DIR/ || die "copying samtools results back"

#copying the fast5 files to local scratch space
echo "Copying fast5. This will take ages."
/usr/bin/time -v cp -r $FAST5DIR $TMPDIR/fast5/ || die "Copying $FAST5DIR to $TMPDIR/fast5/ failed"

#f5c index
echo "Indexing the fast5"
/usr/bin/time -v $F5C index -d $TMPDIR/fast5/ -t $num_threads --iop $num_threads $FASTQ_LOCAL || die "f5c index failed"

#methylation calling and counting methylation frequencies (except ultra-long reads) 
echo "Methylation calling"
/usr/bin/time -v $F5C call-methylation -x hpc-low -t $num_threads --iop $num_threads -r  $FASTQ_LOCAL -g $REF_LOCAL -b $BAM_LOCAL -K 1024 -B 10M --skip-ultra $TMP_BAM_LOCAL > $METH_LOCAL || die "f5c methylation calling failed"
echo "copying methylation calling results"
cp $METH_LOCAL $TMP_BAM_LOCAL $OUTPUT_DIR/ || die "copying methylation calling results failed"
echo "Methylation frequencies"
/usr/bin/time -v $F5C meth-freq -i $METH_LOCAL -s > $METH_FREQ_LOCAL || die "f5c methylation frequency failed"
echo "copying methylation frequency results"
cp $METH_FREQ_LOCAL $OUTPUT_DIR/ || die "copying methylation frequency results failed"

#methylation calling and counting methylation frequencies (ultra-long reads)
echo "Handling ultra-long reads"
samtools index $TMP_BAM_LOCAL
/usr/bin/time -v $F5C call-methylation -t $num_threads --iop $num_threads -r  $FASTQ_LOCAL -g $REF_LOCAL -b $TMP_BAM_LOCAL -K 512 -B 50M  --disable-cuda=yes > $METH_ULTRA_LOCAL || die "f5c methylation calling on ultra-long reads failed"
cp $METH_ULTRA_LOCAL $OUTPUT_DIR/ 
/usr/bin/time -v $F5C meth-freq -i $METH_ULTRA_LOCAL -s > $METH_ULTRA_FREQ_LOCAL || die "f5c methylation frequency on ultra-long reads failed"
cp $METH_ULTRA_FREQ_LOCAL $OUTPUT_DIR/ || die "copying methylation frequency results failed"

#merging methylation frequencies
echo "Merging methylation frequencies"
/usr/bin/time -v $F5C freq-merge $METH_FREQ_LOCAL $METH_ULTRA_FREQ_LOCAL > $METH_FREQ_COMBINED_LOCAL
cp $METH_FREQ_COMBINED_LOCAL $OUTPUT_DIR/

echo "all done"
