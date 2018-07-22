#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <getopt.h>

#include "f5c.h"
#include "f5cmisc.h"

double realtime0=realtime();

static struct option long_options[] = {
	{ "reads",    required_argument, 0, 'r' },
	{ "bam",      required_argument, 0, 'b' },
	{ "genome",   required_argument, 0, 'g' },
	{ "aaaaaa",   no_argument,       0,  0  },
	{ "help",     no_argument,       0, 'h' },
	{ "version",  no_argument,       0, 'V' },
	{ 0, 0, 0, 0}
};

int main(int argc, char *argv[]){

	const char *optstring = "r:b:g:h:v:";
	int longindex=0;
	int c=-1;
	
	char *bamfilename=NULL;
	char *fastafile=NULL;
	char *fastqfile=NULL;
	
	while ((c = getopt_long(argc, argv, optstring, long_options, &longindex)) >= 0){
		if (c == 'r') {
			fastqfile = optarg;
			//if(fastqfile==NULL) fprintf(stderr,ERR "Fastq cannot be null" CEND);
		}
		else if (c == 'b') {
			bamfilename = optarg;
		}
		else if (c == 'g') {
			fastafile = optarg;
		}
		else{
			//puts("dsfgsdgsd");
			//fprintf(stderr,ERR "usage : " CEND);	
			//exit(EXIT_FAILURE);
		}
	}
	
	if(fastqfile==NULL || bamfilename==NULL || fastafile==NULL){
		fprintf(stderr,"Usage: %s [OPTIONS] -r reads.fa -b alignments.bam -g genome.fa\n",argv[0]);
		exit(EXIT_FAILURE);
	}
	
	
	// const char *bamfilename="/mnt/f/share/778Nanopore/fastq/740475-67.bam";
	// const char *fastafile="/mnt/f/share/reference/hg38noAlt.fa";
	// const char *fastqfile="/mnt/f/share/778Nanopore/fastq/740475-67.fastq";

	core_t* core =init_core(bamfilename,fastafile,fastqfile);
	db_t* db=init_db();

	int32_t status;
	while(1){
		status=load_db(core,db);
		if(status<0){
			break;
		}
		fprintf(stderr,"[%s::%.3f*%.2f] %d Entries loaded\n",__func__,realtime()-realtime0,cputime()/(realtime()-realtime0),status);	
		process_db(core,db);
		fprintf(stderr,"[%s::%.3f*%.2f] %d Entries processed\n",__func__,realtime()-realtime0,cputime()/(realtime()-realtime0),status);	
	}

	fprintf(stderr, "\n[%s] Real time: %.3f sec; CPU time: %.3f sec\n", __func__, realtime() - realtime0, cputime());


	return 0;
}
