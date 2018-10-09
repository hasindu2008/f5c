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
//#define HAVE_CUDA 1 //if compile for CUDA or not
#define ALN_BANDWIDTH 100 // the band size in banded_alignment

//option flags
#define F5C_PRINT_RAW 0x001     //print the raw signal to stdio
#define F5C_SECONDARY_YES 0x002 //consider secondary reads
#define F5C_SKIP_UNREADABLE                                                    \
    0x004 //Skip unreadable fast5 and continue rather than exiting
#define F5C_PRINT_EVENTS 0x008
#define F5C_PRINT_BANDED_ALN 0x010
#define F5C_PRINT_SCALING 0x020
#define F5C_DISABLE_CUDA 0x040
#define F5C_DEBUG_BRK 0x080


//flags for a read
#define FAILED_CALIBRATION 0x001 //if the calibration failed
#define FAILED_ALIGNMENT 0x002
#define FAILED_QUALITY_CHK  0x004


#define WORK_STEAL 1
#define STEAL_THRESH 5

#define IO_PROC_INTERLEAVE 1
#define SECTIONAL_BENCHMARK 1   

//#define ALIGN_2D_ARRAY 1 //for CPU whether to use a 1D array or a 2D array



typedef struct {
    int32_t min_mapq;       //minimum mapq
    const char* model_file; //name of the model file
    uint32_t flag;
    int32_t batch_size;
    int32_t num_thread;
    int32_t cuda_block_size;
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

//taken from nanopolish
typedef struct {
    // direct parameters that must be set
    float scale;
    float shift;
    //float drift; = 0 always?
    float var; // set later
    //float scale_sd;
    //float var_sd;

    // derived parameters that are cached for efficiency
    //float log_var;
    //float scaled_var;
    //float log_scaled_var;

} scalings_t;


typedef struct {
        int event_idx;
        int kmer_idx;
} EventKmerPair;

//from nanopolish
typedef struct {
    int ref_pos;
    int read_pos;
} AlignedPair;

//from nanopolish
typedef struct {
    int32_t start;
    int32_t stop; // inclusive
} index_pair_t;


//from nanopolish
typedef struct {
    // ref data
    //char* ref_name;
    char ref_kmer[KMER_SIZE + 1];
    int32_t ref_position;

    // event data
    int32_t read_idx;
    //int32_t strand_idx;
    int32_t event_idx;
    bool rc;

    // hmm data
    char model_kmer[KMER_SIZE + 1];
    char hmm_state;
} event_alignment_t;

typedef struct {
    // region string
    //char* region;

    // bam records
    bam1_t** bam_rec;
    int32_t capacity_bam_rec; // will these overflow?
    int32_t n_bam_rec;

    // fasta cache //can optimise later by caching a common string for all
    // records in the batch
    char** fasta_cache;

    //read sequence //todo : optimise by grabbing it from bam seq. is it possible due to clipping?
    char** read;
    int32_t* read_len;

    // fast5 file //should flatten this to reduce mallocs
    fast5_t** f5;

    //event table
    event_table* et;

    //scalings
    scalings_t* scalings;

    //aligned pairs
    AlignedPair** event_align_pairs;
    int32_t* n_event_align_pairs;

    //event alignments
    event_alignment_t** event_alignment;
    int32_t* n_event_alignment;
    double* events_per_base; //todo : do we need double?

    index_pair_t** base_to_event_map;

    int32_t* read_stat_flag;
    

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
    model_t* model; //dna model
    model_t* cpgmodel;

    // options
    opt_t opt;

    //realtime0
    double realtime0;
    double align_time;

} core_t;

typedef struct {
    core_t* core;
    db_t* db;
    int32_t starti;
    int32_t endi;
    void (*func)(core_t*,db_t*,int);
#ifdef WORK_STEAL
    void *all_pthread_args;
#endif
} pthread_arg_t;

typedef struct {
    core_t* core;
    db_t* db;
    //conditional variable for notifying the processing to the output threads
    pthread_cond_t cond;
    pthread_mutex_t mutex;
} pthread_arg2_t;

db_t* init_db(core_t* core);
int32_t load_db(core_t* dg, db_t* db);
core_t* init_core(const char* bamfilename, const char* fastafile,
                  const char* fastqfile, opt_t opt,double realtime0);
void process_db(core_t* dg, db_t* db);
void pthread_db(core_t* core, db_t* db, void (*func)(core_t*,db_t*,int));
void align_db(core_t* core, db_t* db);
void output_db(core_t* core, db_t* db);
void free_core(core_t* core);
void free_db_tmp(db_t* db);
void free_db(db_t* db);
void init_opt(opt_t* opt);

#endif
