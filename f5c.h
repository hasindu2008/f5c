#ifndef F5C_H
#define F5C_H

#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/faidx.h>

#include "fast5lite.h"
#include "nanopolish_read_db.h"
 
//from scrappie
typedef struct {
	uint64_t start;
	float length;
	float mean;
	float stdv;
	int pos;
	int state;
} event_t;

//from scrappie
typedef struct {
	size_t n;
	size_t start;
	size_t end;
	event_t *event;
} event_table;


struct data_batch_t{

	//region string
	char *region;

	//bam records
	bam1_t** bam_rec;
	int capacity_bam_rec;
	int n_bam_rec;

	//fasta cache //can optimise later by caching a common string for all records in th ebatch
	char** fasta_cache;

	//fast5 file
	fast5** f5;
};
typedef struct data_batch_t data_batch;

struct data_global_t{
	
	//bam file related
	htsFile* m_bam_fh;
    hts_idx_t* m_bam_idx;
    bam_hdr_t* m_hdr;
    hts_itr_t* itr;

    //fa related
    faidx_t *fai;

    //readbb
    ReadDB *readbb;
	
};
typedef struct data_global_t data_global;


data_batch* init_databatch();
int load_databatch(data_batch* db,data_global* dg);
data_global* init_files(const char *bamfilename, const char *fastafile,const char *fastqfile);
void* process_databatch(data_batch* db,data_global* dg);

#endif