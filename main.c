#include <stdio.h>
#include <stdlib.h>

#include "f5c.h"
#include "f5cmisc.h"

double realtime0=realtime();

int main(){

	const char *bamfilename="/mnt/f/share/778Nanopore/fastq/740475-67.bam";
	const char *fastafile="/mnt/f/share/reference/hg38noAlt.fa";
	const char *fastqfile="/mnt/f/share/778Nanopore/fastq/740475-67.fastq";

	data_global* dg =init_files(bamfilename,fastafile,fastqfile);
	data_batch* db=init_databatch();

	int status;
	while(1){
		status=load_databatch(db,dg);
		if(status<0){
			break;
		}
		fprintf(stderr,"[%s::%.3f*%.2f] %d Entries loaded\n",__func__,realtime()-realtime0,cputime()/(realtime()-realtime0),status);	
		process_databatch(db,dg);
		fprintf(stderr,"[%s::%.3f*%.2f] %d Entries processed\n",__func__,realtime()-realtime0,cputime()/(realtime()-realtime0),status);	
	}

	fprintf(stderr, "\n[%s] Real time: %.3f sec; CPU time: %.3f sec\n", __func__, realtime() - realtime0, cputime());


	return 0;
}
