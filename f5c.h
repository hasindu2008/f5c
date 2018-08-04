#ifndef F5C_H
#define F5C_H

#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/sam.h>

#include "fast5lite.h"
#include "nanopolish_read_db.h"

#define m_min_mapping_quality 30
#define KMER_SIZE 6 //hard coded for now; todo : change to dynamic
#define NUM_KMER 4096

//flags
#define F5C_PRINT_RAW 0x001     //print the raw signal to stdio
#define F5C_SECONDARY_YES 0x002 //consider secondary reads
#define F5C_SKIP_UNREADABLE                                                    \
    0x004 //Skip unreadable fast5 and continue rather than exiting
#define F5C_PRINT_EVENTS 0x008

typedef struct {
    int32_t min_mapq;       //minimum mapq
    const char* model_file; //name of the model file
    uint32_t flag;
    int32_t batch_size;
} opt_t;

// from scrappie
typedef struct {
    uint64_t start;
    float length; //cant be made int
    float mean;
    float stdv;
    int32_t pos;   //always -1 can be removed
    int32_t state; //always -1 can be removed
} event_t;

// from scrappie
typedef struct {
    size_t n;
    size_t start; //always 0
    size_t end;   //always eqial to n
    event_t* event;
} event_table;

//model
typedef struct {
    //char kmer[KMER_SIZE + 1]; //KMER_SIZE+null character //can get rid of this
    float level_mean;
    float level_stdv;
    float sd_mean;
    float sd_stdv;
    //float weight;
} model_t;

typedef struct {
    // region string
    char* region;

    // bam records
    bam1_t** bam_rec;
    int32_t capacity_bam_rec; // will these overflow?
    int32_t n_bam_rec;

    // fasta cache //can optimise later by caching a common string for all
    // records in the batch
    char** fasta_cache;

    // fast5 file //should flatten this to reduce mallocs
    fast5_t** f5;

    //event table
    event_table* et;

} db_t;

typedef struct {
    // bam file related
    htsFile* m_bam_fh;
    hts_idx_t* m_bam_idx;
    bam_hdr_t* m_hdr;
    hts_itr_t* itr;

    // fa related
    faidx_t* fai;

    // readbb
    ReadDB* readbb;

    // models
    model_t* model;

    // options
    opt_t opt;

} core_t;

db_t* init_db(core_t* core);
int32_t load_db(core_t* dg, db_t* db);
core_t* init_core(const char* bamfilename, const char* fastafile,
                  const char* fastqfile, opt_t opt);
void process_db(core_t* dg, db_t* db);
void output_db(core_t* core, db_t* db);
void free_core(core_t* core);
void free_db_tmp(db_t* db);
void free_db(db_t* db);
void init_opt(opt_t* opt);

#endif
