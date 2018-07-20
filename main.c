
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include <errno.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/faidx.h>
#include "fast5lite.h"
#include "common.h"
#include "nanopolish_read_db.h"
#include "main.h"

#define m_min_mapping_quality 30

data_batch* init_databatch(){

	data_batch* db=(data_batch *)(malloc(sizeof(data_batch)));
	errorCheckNULL(db);

	db->capacity_bam_rec=512;
	db->n_bam_rec=0;

	db->bam_rec=(bam1_t**)(malloc(sizeof(bam1_t*) * db->capacity_bam_rec));
	errorCheckNULL(db->bam_rec);

	int	 i=0;
	for(i = 0; i < db->capacity_bam_rec; ++i) {
		db->bam_rec[i] = bam_init1();
		errorCheckNULL(db->bam_rec[i]);
	}

	db->fasta_cache=(char**)(malloc(sizeof(char*) * db->capacity_bam_rec));
	errorCheckNULL(db->fasta_cache);

	db->f5=(fast5**)malloc(sizeof(fast5*) * db->capacity_bam_rec);

	return db;
}

int load_databatch(data_batch* db,data_global* dg){


	//get bams
	bam1_t* record;
	int result=0;
	db->n_bam_rec=0;
	int i=0;
	while(db->n_bam_rec < db->capacity_bam_rec){
		record=db->bam_rec[db->n_bam_rec];
		result = sam_itr_next(dg->m_bam_fh, dg->itr, record);

		if(result<0){
			break;
		}
		else{
			if ((record->core.flag & BAM_FUNMAP) == 0 && record->core.qual >= m_min_mapping_quality){
				//printf("%s\t%d\n",bam_get_qname(db->bam_rec[db->n_bam_rec]),result);
				db->n_bam_rec++;
			}

		}



	}
	fprintf(stderr,"%d queries read\n",db->n_bam_rec);


	//get refsequences (can emake efficient by taking the the start and end of the sorted bam)
	for(i=0;i<db->n_bam_rec;i++){
		bam1_t* record=db->bam_rec[i];
		char *ref_name=dg->m_hdr->target_name[record->core.tid];
		//printf("refname : %s\n",ref_name);
		int ref_start_pos = record->core.pos;
		int ref_end_pos =  bam_endpos(record);
		assert(ref_end_pos >= ref_start_pos);

		// Extract the reference sequence for this region
		int fetched_len = 0;
		char *refseq=faidx_fetch_seq(dg->fai, ref_name, ref_start_pos, ref_end_pos, &fetched_len);
		db->fasta_cache[i]=refseq;
	   // printf("seq : %s\n",db->fasta_cache[i]);

		//get the fast5

		// Get the read type from the fast5 file
		std::string qname=bam_get_qname(db->bam_rec[i]);
		char *fast5_path=(char *)malloc(1024); //this is really horrible, Need to get size
		strcpy(fast5_path,dg->readbb->get_signal_path(qname).c_str());

		hid_t hdf5_file = fast5_open(fast5_path);
		if(hdf5_file>=0){
			db->f5[i]=(fast5*)calloc(1,sizeof(fast5));
			fast5_read(hdf5_file, db->f5[i]);
			fast5_close(hdf5_file);
		}
		else{
			 fprintf(stderr, "[warning] fast5 file is unreadable and will be skipped: %s\n", fast5_path);
		}

		// printf("%s : %s : %ld\n",qname.c_str(),fast5_path,db->f5[i]->nsample);
	 //	   int j=0;
	 //	   for(j=0;j<db->f5[i]->nsample;j++){
	 //		   printf("%f ",db->f5[i]->rawptr[j]);
	 //	   }
	 //	   printf("\n");
	}
	fprintf(stderr,"%d fast5 read\n",db->n_bam_rec);


	return result;
}



// destroy_databatch();



data_global* init_files(const char *bamfilename, const char *fastafile,const char *fastqfile){

	data_global* dg=(data_global *)malloc(sizeof(data_global));
	errorCheckNULL(dg);

	// load bam file
	dg->m_bam_fh = sam_open(bamfilename, "r");
	errorCheckNULL(dg->m_bam_fh);

	// load bam index file
	dg->m_bam_idx = sam_index_load(dg->m_bam_fh, bamfilename);
	errorCheckNULL(dg->m_bam_idx);

	// read the bam header
	dg->m_hdr = sam_hdr_read(dg->m_bam_fh);
	errorCheckNULL(dg->m_hdr);


	dg->itr = sam_itr_queryi(dg->m_bam_idx, HTS_IDX_START, 0, 0);
	errorCheckNULL(dg->itr);


	//reference file
	dg->fai = fai_load(fastafile); //error check done inside inside

	//readbb
	dg->readbb=new ReadDB;
	dg->readbb->load(fastqfile);



	return dg;
}

int process_databatch(data_batch* db,data_global* dg);

int main(){

	const char *bamfilename="/mnt/f/share/778Nanopore/fastq/740475-67.bam";
	const char *fastafile="/mnt/f/share/reference/hg38noAlt.fa";
	const char *fastqfile="/mnt/f/share/778Nanopore/fastq/740475-67.fastq";

	data_global* dg =init_files(bamfilename,fastafile,fastqfile);
	data_batch* db=init_databatch();

	while(load_databatch(db,dg)>=0){
			process_databatch(db,dg);
			fprintf(stderr,"Processed\n") ;
	}


	// // Initialize iteration
	// std::vector<bam1_t*> records(m_batch_size, NULL);
	// for(size_t i = 0; i < records.size(); ++i) {
	//	   records[i] = bam_init1();
	// }

	// int result;
	// size_t num_reads_realigned = 0;
	// size_t num_records_buffered = 0;




	// do {
	//	   assert(num_records_buffered < records.size());

	//	   // read a record into the next slot in the buffer
	//	   result = sam_itr_next(m_bam_fh, itr, records[num_records_buffered]);
	//	   num_records_buffered += result >= 0;

	//	   // realign if we've hit the max buffer size or reached the end of file
	//	   if(num_records_buffered == records.size() || result < 0 || (num_records_buffered + num_reads_realigned == m_max_reads)) {
	//		   #pragma omp parallel for schedule(dynamic)
	//		   for(size_t i = 0; i < num_records_buffered; ++i) {
	//			   bam1_t* record = records[i];
	//			   size_t read_idx = num_reads_realigned + i;
	//			   if( (record->core.flag & BAM_FUNMAP) == 0 && record->core.qual >= m_min_mapping_quality) {
	//				   func(m_hdr, record, read_idx, clip_start, clip_end);
	//			   }
	//		   }

	//		   num_reads_realigned += num_records_buffered;
	//		   num_records_buffered = 0;
	//	   }
	// } while(result >= 0 && num_reads_realigned < m_max_reads);


	return 0;
}
