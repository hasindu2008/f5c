#!/bin/bash
#PBS -N F5C-METHCALL
#PBS -q gpuvolta
#PBS -l ncpus=12
#PBS -l ngpus=1
#PBS -l mem=96GB
#PBS -l jobfs=350GB
#PBS -l walltime=48:00:00
#PBS -l wd

###################################################################

## example PBS script - resource efficient methylation calling pipeline for a dataset with many ultra-long reads
## tested on the nci-gadi cluster (on a GPU node equipped with Tesla V100 - 32 GB GPUs) [https://nci.org.au/our-systems/hpc-systems]
## 12 CPU cores have been requested using -l ncpus=12
## 1 GPU has been requested using -l ngpus=1
## 96 GB of RAM has been requested using -l mem=96GB 
## 350 GB of local scratch storage space (fast SSD storage on the node) has been requested using -l jobfs=350GB
## note that mem=96GB and jobfs=350GB values are the totals (i.e. NOT per single core)
## it is assumed that the compute node has adequate local scratch storage space and is accessible via $PBS_JOBFS
## the local scratch space should be adequate to hold the reference genome, FASTQ file and the generated f5c index files (~0.6 of FASTQ size)

#change this to where the fast5 files are located on the shared cluster storage (relative or absolute path)
FAST5DIR=fast5/

#change this to where the FASTQ file is located on the shared cluster storage (relative or absolute path)
#this file will be later copied to the local scratch storage
FASTQFILE=na12878_rel6_all.fastq

#change this to the directory where you want the output to be written on the shared cluster storage (relative or absolute path)
#the directory will be created if not already existing
#the SAM file (reads.sam), BAM file (reads.bam, tmp.bam), per-read methylation calls (meth.tsv, meth-ultra.tsv), methylation frequencies (meth-freq.tsv, meth-ultra-freq.tsv, meth-freq-combined.tsv) will be written to here.
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

###################################################################
###################################################################

num_threads=$PBS_NCPUS

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

#temporary location of the reference and the FASTQ file on local scratch storage where they are to be copied
REF_LOCAL=$PBS_JOBFS/ref.fa
FASTQ_LOCAL=$PBS_JOBFS/reads.fastq

#output files
SAM=$OUTPUT_DIR/reads.sam
BAM=$OUTPUT_DIR/reads.bam
TMP_BAM=$OUTPUT_DIR/tmp.bam
METH=$OUTPUT_DIR/meth.tsv
METH_ULTRA=$OUTPUT_DIR/meth-ultra.tsv
METH_FREQ=$OUTPUT_DIR/meth-freq.tsv
METH_ULTRA_FREQ=$OUTPUT_DIR/meth-ultra-freq.tsv
METH_FREQ_COMBINED=$OUTPUT_DIR/meth-freq-combined.tsv

#copy the reference and the FASTQ file to the local scratch storage
echo "Copying reference $REFERENCE to $REF_LOCAL"
cp $REFERENCE  $REF_LOCAL || die "Copying $REFERENCE to $REF_LOCAL failed"
echo "Copying $FASTQFILE to $FASTQ_LOCAL"
cp $FASTQFILE  $FASTQ_LOCAL || die "Copying $FASTQFILE to $FASTQ_LOCAL failed"	

#minimap
echo "Minimap2 running"
/usr/bin/time -v $MINIMAP -x map-ont -a -t$num_threads --secondary=no $REF_LOCAL $FASTQ_LOCAL > $SAM  || die "minimap2 failed"

#sorting and indexing the BAM
echo "samtools sorting"
/usr/bin/time -v $SAMTOOLS sort -@$num_threads $SAM > $BAM || die "samtools failed"
/usr/bin/time -v $SAMTOOLS index $BAM || die "samtools index failed"

#f5c index
echo "Indexing the fast5"
/usr/bin/time -v $F5C index -d $FAST5DIR -t $num_threads --iop 64 $FASTQ_LOCAL || die "f5c index failed"

#methylation calling and counting methylation frequencies (except ultra-long reads) 
echo "Methylation calling"
/usr/bin/time -v $F5C call-methylation -x nci-gadi -t $num_threads --iop 64 -r $FASTQ_LOCAL -g $REF_LOCAL -b $BAM -K 2048 -B 20M --skip-ultra $TMP_BAM > $METH || die "f5c methylation calling failed"
echo "Methylation frequencies"
/usr/bin/time -v $F5C meth-freq -i $METH -s > $METH_FREQ || die "f5c methylation frequency failed"

#methylation calling and counting methylation frequencies (ultra-long reads)
echo "Handling ultra-long reads"
/usr/bin/time -v $SAMTOOLS index $TMP_BAM
/usr/bin/time -v $F5C call-methylation -t $num_threads --iop 64 -r  $FASTQ_LOCAL -g $REF_LOCAL -b $TMP_BAM -K 512 -B 20M --disable-cuda=yes > $METH_ULTRA || die "f5c methylation calling on ultra-long reads failed"
/usr/bin/time -v $F5C meth-freq -i $METH_ULTRA -s > $METH_ULTRA_FREQ || die "f5c methylation frequency on ultra-long reads failed"

#merging methylation frequencies
echo "Merging methylation frequencies"
/usr/bin/time -v $F5C freq-merge $METH_FREQ $METH_ULTRA_FREQ > $METH_FREQ_COMBINED

echo "all done"
