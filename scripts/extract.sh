

awk 'NR%4==1{a=$0} NR%4==2{b=$0} NR%4==3{c=$0} NR%4==0&&length(b)<30000{print a"\n"b"\n"c"\n"$0;}' file.fq > result.fq

minimap2 -a -x map-ont -t12 humangenome.fa reads30k.fq | samtools sort -T tmp -o reads30k.bam
samtools index reads30k.bam

nanopolish index -d test/chr22_meth_example/fast5_files/ test/chr22_meth_example/reads30k.fq
