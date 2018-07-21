#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "f5c.h"
#include "f5cmisc.h"

#define m_min_mapping_quality 30

data_batch* init_databatch(){

	data_batch* db=(data_batch *)(malloc(sizeof(data_batch))); errorCheckNULL(db);

	db->capacity_bam_rec=512;
	db->n_bam_rec=0;

	db->bam_rec=(bam1_t**)(malloc(sizeof(bam1_t*) * db->capacity_bam_rec));	MALLOC_CHK(db->bam_rec);

	int	 i=0;
	for(i = 0; i < db->capacity_bam_rec; ++i) {
		db->bam_rec[i] = bam_init1(); // does bam_init1 already perform an error check
		errorCheckNULL(db->bam_rec[i]);
	}

	db->fasta_cache=(char**)(malloc(sizeof(char*) * db->capacity_bam_rec)); MALLOC_CHK(db->fasta_cache);
	db->f5=(fast5**)malloc(sizeof(fast5*) * db->capacity_bam_rec); MALLOC_CHK(db->f5);

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
			if ((record->core.flag & BAM_FUNMAP) == 0 && record->core.qual >= m_min_mapping_quality){ //remove secondraies?
				//printf("%s\t%d\n",bam_get_qname(db->bam_rec[db->n_bam_rec]),result);
				db->n_bam_rec++;
			}

		}
	}
	//fprintf(stderr,"%s:: %d queries read\n",__func__,db->n_bam_rec);

	//get ref sequences (can make efficient by taking the the start and end of the sorted bam)
	for(i=0;i<db->n_bam_rec;i++){
		bam1_t* record=db->bam_rec[i];
		char *ref_name=dg->m_hdr->target_name[record->core.tid];
		//printf("refname : %s\n",ref_name);
		int ref_start_pos = record->core.pos;
		int ref_end_pos =  bam_endpos(record);
		assert(ref_end_pos >= ref_start_pos);

		// Extract the reference sequence for this region
		int fetched_len = 0;
		char *refseq=faidx_fetch_seq(dg->fai, ref_name, ref_start_pos, ref_end_pos, &fetched_len); //error handle?
		db->fasta_cache[i]=refseq;
		// printf("seq : %s\n",db->fasta_cache[i]);

		//get the fast5

		// Get the read type from the fast5 file
		std::string qname=bam_get_qname(db->bam_rec[i]);
		char *fast5_path=(char *)malloc(dg->readbb->get_signal_path(qname).size()+10); //is +10 needed? do errorcheck
		strcpy(fast5_path,dg->readbb->get_signal_path(qname).c_str());

		hid_t hdf5_file = fast5_open(fast5_path);
		if(hdf5_file>=0){
			db->f5[i]=(fast5*)calloc(1,sizeof(fast5)); //todo : errorcheck
			fast5_read(hdf5_file, db->f5[i]); //todo : errorhandle
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
	//fprintf(stderr,"%s:: %d fast5 read\n",__func__,db->n_bam_rec);


	if(result<0){
		return result;
	}
	return db->n_bam_rec;
}


//todo : 
// destroy_databatch();

data_global* init_files(const char *bamfilename, const char *fastafile,const char *fastqfile){

	data_global* dg=(data_global *)malloc(sizeof(data_global)); MALLOC_CHK(dg);

	// load bam file
	dg->m_bam_fh = sam_open(bamfilename, "r"); errorCheckNULL(dg->m_bam_fh);

	// load bam index file
	dg->m_bam_idx = sam_index_load(dg->m_bam_fh, bamfilename); errorCheckNULL(dg->m_bam_idx);

	// read the bam header
	dg->m_hdr = sam_hdr_read(dg->m_bam_fh); errorCheckNULL(dg->m_hdr); 

	dg->itr = sam_itr_queryi(dg->m_bam_idx, HTS_IDX_START, 0, 0); errorCheckNULL(dg->itr);

	//reference file
	dg->fai = fai_load(fastafile); //error check done inside?

	//readbb
	dg->readbb=new ReadDB;
	dg->readbb->load(fastqfile);

	return dg;
}



void* process_databatch(data_batch* db,data_global* dg){

	event_table* et=(event_table*)malloc(sizeof(event_table)*db->n_bam_rec); MALLOC_CHK(et);

	int i;
	for(i=0 ; i<db->n_bam_rec ; i++){

		float *rawptr=db->f5[i]->rawptr;
		float range=db->f5[i]->range;
		float digitisation=db->f5[i]->digitisation;
		float offset = db->f5[i]->offset;
		int nsample=db->f5[i]->nsample;

		// convert to pA
		float raw_unit = range / digitisation;
		for (int j = 0; j < nsample; j++) {
			rawptr[j] = (rawptr[j] + offset) * raw_unit;
		}
		et[i]=getevents(db->f5[i]->nsample,rawptr);
	}

	return (void*)et;
}





