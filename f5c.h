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
	int32_t pos;
	int32_t state;
} event_t;

//from scrappie
typedef struct {
	size_t n;
	size_t start;
	size_t end;
	event_t *event;
} event_table;


typedef struct{

	//region string
	char *region;

	//bam records
	bam1_t** bam_rec;
	int32_t capacity_bam_rec; //will these overflow?
	int32_t n_bam_rec;

	//fasta cache //can optimise later by caching a common string for all records in th ebatch
	char** fasta_cache;

	//fast5 file
	fast5_t** f5;
	
} db_t;


typedef struct{
	
	//bam file related
	htsFile* m_bam_fh;
    hts_idx_t* m_bam_idx;
    bam_hdr_t* m_hdr;
    hts_itr_t* itr;

    //fa related
    faidx_t *fai;

    //readbb
    ReadDB *readbb;
	
	//options
	int8_t print;
	int8_t print_events;
	
} core_t;



db_t* init_db();
int32_t load_db(core_t* dg,db_t* db);
core_t* init_core(const char *bamfilename, const char *fastafile,const char *fastqfile);
void* process_db(core_t* dg,db_t* db);

#endif