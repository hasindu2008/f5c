if [ $# -eq 0 ]
then
	cd /mnt/f/share/778Nanopore/fastq/ && /mnt/f/hasindug.git/f5c/f5c -b /mnt/f/share/778Nanopore/fastq/740475-67.bam -g /mnt/f/share/reference/hg38noAlt.fa -r /mnt/f/share/778Nanopore/fastq/740475-67.fastq
fi

if [ $# -eq 1 ]
then
	if [ $1 == "valgrind" ]
	then
		cd /mnt/f/share/778Nanopore/fastq/ && valgrind /mnt/f/hasindug.git/f5c/f5c -b /mnt/f/share/778Nanopore/fastq/740475-67.bam -g /mnt/f/share/reference/hg38noAlt.fa -r /mnt/f/share/778Nanopore/fastq/740475-67.fastq
	elif [ $1 == "gdb" ]
	then
		cd /mnt/f/share/778Nanopore/fastq/ && gdb --args /mnt/f/hasindug.git/f5c/f5c -b /mnt/f/share/778Nanopore/fastq/740475-67.bam -g /mnt/f/share/reference/hg38noAlt.fa -r /mnt/f/share/778Nanopore/fastq/740475-67.fastq
	else
		echo "Wrong usage"
	fi	
fi




