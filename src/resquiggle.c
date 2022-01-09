/* @file resquiggle.c
**
** align raw signal to basecalled read
** @author: Hiruna Samarakoon (hiruna@unsw.edu.au)
** @@author: Hasindu Gamaarachchi (hasindu@unsw.edu.au)
******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <unistd.h>
#include <assert.h>
#include <getopt.h>

#include "error.h"
#include "f5c.h"
#include "f5cmisc.h"
#include "kseq.h"

#define TO_PICOAMPS(RAW_VAL,DIGITISATION,OFFSET,RANGE) (((RAW_VAL)+(OFFSET))*((RANGE)/(DIGITISATION)))

KSEQ_INIT(gzFile, gzread)

static const char *RESQUIGGLE_USAGE_MESSAGE =
    "Usage: f5c signal-read-align [OPTIONS] [SLOW5_FILE/DIR] reads.fastq ...\n"
    "Align raw signal to the basecalled read.\n\n"
    "   -o FILE         output file. Write to stdout if not specified\n"
    "   -h              help\n"
    "   --version       print version\n"
    "\nSee the manual page for details (`man ./docs/f5c.1' or https://f5c.page.link/man)."
    "\n\n"
    ;


/* initialise the core data structure */
core_t* init_core2(opt_t opt, int8_t mode) {

    core_t* core = (core_t*)malloc(sizeof(core_t));
    MALLOC_CHK(core);

    core->opt = opt;

    //model
    core->model = (model_t*)malloc(sizeof(model_t) * MAX_NUM_KMER); //4096 is 4^6 which is hardcoded now
    MALLOC_CHK(core->model);
    core->cpgmodel = (model_t*)malloc(sizeof(model_t) * MAX_NUM_KMER_METH); //15625 is 4^6 which os hardcoded now
    MALLOC_CHK(core->cpgmodel);

    //load the model from files
    uint32_t kmer_size=0;
    uint32_t kmer_size_meth=0;
    if (opt.model_file) {
        kmer_size=read_model(core->model, opt.model_file, MODEL_TYPE_NUCLEOTIDE);
    } else {
        if(opt.flag & F5C_RNA){
            INFO("%s","builtin RNA nucleotide model loaded");
            kmer_size=set_model(core->model, MODEL_ID_RNA_NUCLEOTIDE);
        }
        else{
            kmer_size=set_model(core->model, MODEL_ID_DNA_NUCLEOTIDE);
        }
    }
    if (opt.meth_model_file) {
        kmer_size_meth=read_model(core->cpgmodel, opt.meth_model_file, MODEL_TYPE_METH);
    } else {
        kmer_size_meth=set_model(core->cpgmodel, MODEL_ID_DNA_CPG);
    }
    if( mode==0 && kmer_size != kmer_size_meth){
        ERROR("The k-mer size of the nucleotide model (%d) and the methylation model (%d) should be the same.",kmer_size,kmer_size_meth);
        exit(EXIT_FAILURE);
    }
    core->kmer_size = kmer_size;
    return core;

}

void free_core2(core_t* core) {
    free(core->model);
    free(core->cpgmodel);
//    delete core->readbb;

    free(core);
}

/* initialise a data batch */
db_t* init_db2(core_t* core) {
    int32_t i = 0;
    db_t* db = (db_t*)(malloc(sizeof(db_t)));
    MALLOC_CHK(db);

    db->capacity_bam_rec = core->opt.batch_size;
    db->n_bam_rec = 0;

    db->read = (char**)(malloc(sizeof(char*) * db->capacity_bam_rec));
    MALLOC_CHK(db->read);
    db->read_len = (int32_t*)(malloc(sizeof(int32_t) * db->capacity_bam_rec));
    MALLOC_CHK(db->read_len);
//    db->read_idx = (int64_t*)(malloc(sizeof(int64_t) * db->capacity_bam_rec));
//    MALLOC_CHK(db->read_idx);

    db->f5 = (fast5_t**)malloc(sizeof(fast5_t*) * db->capacity_bam_rec);
    MALLOC_CHK(db->f5);

    db->et = (event_table*)malloc(sizeof(event_table) * db->capacity_bam_rec);
    MALLOC_CHK(db->et);

    db->scalings =
        (scalings_t*)malloc(sizeof(scalings_t) * db->capacity_bam_rec);
    MALLOC_CHK(db->scalings);

    db->event_align_pairs =
        (AlignedPair**)malloc(sizeof(AlignedPair*) * db->capacity_bam_rec);
    MALLOC_CHK(db->event_align_pairs);
    db->n_event_align_pairs =
        (int32_t*)malloc(sizeof(int32_t) * db->capacity_bam_rec);
    MALLOC_CHK(db->n_event_align_pairs);

    db->event_alignment = (event_alignment_t**)malloc(
        sizeof(event_alignment_t*) * db->capacity_bam_rec);
    MALLOC_CHK(db->event_alignment);
    db->n_event_alignment =
        (int32_t*)malloc(sizeof(int32_t) * db->capacity_bam_rec);
    MALLOC_CHK(db->n_event_alignment);

    db->events_per_base =
        (double*)malloc(sizeof(double) * db->capacity_bam_rec);
    MALLOC_CHK(db->events_per_base);

    db->base_to_event_map =
        (index_pair_t**)malloc(sizeof(index_pair_t*) * db->capacity_bam_rec);
    MALLOC_CHK(db->base_to_event_map);

    db->read_stat_flag = (int32_t *)malloc(sizeof(int32_t) * db->capacity_bam_rec);
    MALLOC_CHK(db->read_stat_flag);

    db->site_score_map = (std::map<int, ScoredSite> **)malloc(sizeof(std::map<int, ScoredSite> *) * db->capacity_bam_rec);
    MALLOC_CHK(db->site_score_map);

    db->total_reads=0;
    db->bad_fast5_file=0;
    db->ultra_long_skipped=0;

    //eventalign related
//    if(core->mode==1){
    if(0){
        db->eventalign_summary = (EventalignSummary *)malloc(sizeof(EventalignSummary) * db->capacity_bam_rec);
        MALLOC_CHK(db->eventalign_summary);

        db->event_alignment_result = (std::vector<event_alignment_t> **)malloc(sizeof(std::vector<event_alignment_t> *) * db->capacity_bam_rec);
        MALLOC_CHK(db->event_alignment_result);

        db->event_alignment_result_str = (char **)malloc(sizeof(char *) * db->capacity_bam_rec);
        MALLOC_CHK(db->event_alignment_result_str);

        for (i = 0; i < db->capacity_bam_rec; ++i) {
            db->event_alignment_result[i] = new std::vector<event_alignment_t> ;
            NULL_CHK(db->event_alignment_result[i]);
            (db->eventalign_summary[i]).num_events=0; //done here in the same loop for efficiency
            db->event_alignment_result_str[i] = NULL;
        }

    }
    else{
        db->eventalign_summary = NULL;
        db->event_alignment_result = NULL;
        db->event_alignment_result_str = NULL;
    }

    return db;
}

void free_db2(db_t* db) {
    int32_t i = 0;
    free(db->read);
    free(db->read_len);
    free(db->et);
    free(db->f5);
    free(db->scalings);
    free(db->event_align_pairs);
    free(db->n_event_align_pairs);
    free(db->event_alignment);
    free(db->n_event_alignment);
    free(db->events_per_base);
    free(db->base_to_event_map);
    free(db->read_stat_flag);
//    for (i = 0; i < db->capacity_bam_rec; ++i) {
//        delete db->site_score_map[i];
//    }
    free(db->site_score_map);
    //eventalign related
    if(db->eventalign_summary){
        free(db->eventalign_summary);
    }
    if(db->event_alignment_result){
        for (i = 0; i < db->capacity_bam_rec; ++i) {
            delete db->event_alignment_result[i];
        }
        free(db->event_alignment_result);
    }
    if(db->event_alignment_result_str){
        free(db->event_alignment_result_str);
    }
    free(db);
}

/* partially free a data batch - only the read dependent allocations are freed */
void free_db_tmp2(db_t* db) {
    int32_t i = 0;
    for (i = 0; i < db->capacity_bam_rec; ++i) {
        free(db->read[i]);
        free(db->f5[i]->rawptr);
        free(db->f5[i]);
        free(db->et[i].event);
        free(db->event_align_pairs[i]);
//        free(db->base_to_event_map[i]);
//        delete db->site_score_map[i];
//        db->site_score_map[i] = new std::map<int, ScoredSite>;

        if(db->event_alignment_result){ //eventalign related
            delete db->event_alignment_result[i];
            db->event_alignment_result[i] = new std::vector<event_alignment_t>;
        }
        if(db->event_alignment_result_str){ //eventalign related
            free(db->event_alignment_result_str[i]);
            db->event_alignment_result_str[i]=NULL;
        }
    }
}

void process_single2(core_t* core, db_t* db,int32_t i) {
    event_single(core,db,i);
    align_single(core, db,i);
//    scaling_single(core,db,i);
//    meth_single(core,db,i);
}

void process_db2(core_t* core, db_t* db) {
    if (core->opt.num_thread == 1) {
        int32_t i=0;
        for (i = 0; i < db->capacity_bam_rec; i++) {
            process_single2(core,db,i);
        }
    }
    else {
        pthread_db(core,db,process_single2);
    }
}

int resquiggle_main(int argc, char **argv) {
    fprintf(stderr,"resquiggle_main\n");

    // No arguments given
    if (argc <= 1) {
        fprintf(stderr, RESQUIGGLE_USAGE_MESSAGE, argv[0]);
        exit(EXIT_FAILURE);
    }

    opt_t opt;
    init_opt(&opt); //initialise options to defaults
//    opt.num_thread = 2;

    // open fastq
    gzFile fp;
    kseq_t *seq;
    int ret_kseq_read;
    if (argc == 1) {
        fprintf(stderr, "Usage: %s <in.fasta>\n", argv[1]);
        return 1;
    }
    fp = gzopen(argv[1], "r");
    seq = kseq_init(fp);

    //open slow5
    slow5_file_t *sp = slow5_open(argv[2],"r");
    if(sp==NULL){
        fprintf(stderr,"Error in opening file\n");
        exit(EXIT_FAILURE);
    }
    slow5_rec_t *rec = NULL;
    int ret=0;

    ret = slow5_idx_load(sp);
    if(ret<0){
        fprintf(stderr,"Error in loading index\n");
        exit(EXIT_FAILURE);
    }


    int64_t batch_size = opt.batch_size;

    //initialise the core data structure
    core_t* core = init_core2(opt, 1);

    int flag_EOF = 0;
    while (1){
        db_t* db_t = init_db2(core);
        int64_t record_count = 0;
        while (record_count < batch_size) {
            if ((ret_kseq_read = kseq_read(seq)) < 0) {
                if (ret_kseq_read < -1) {
                    return EXIT_FAILURE;
                } else { //EOF file reached
                    flag_EOF = 1;
                    db_t->capacity_bam_rec = record_count;
                    db_t->n_bam_rec = record_count;
                    break;
                }
            } else {
                db_t->read[record_count] = strdup(seq->seq.s);
                db_t->read_len[record_count] = strlen(db_t->read[record_count]);
                ret = slow5_get(seq->name.s, &rec, sp);
                if(ret < 0){
                    fprintf(stderr,"Error in when fetching the read\n");
                }
                else{

                    //        printf("name: %s\n", seq->name.s);
                    //        if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
                    //        printf("seq: %s\n", seq->seq.s);
                    //        if (seq->qual.l) printf("qual: %s\n", seq->qual.s);

                    printf("%s\t",rec->read_id);
                    uint64_t len_raw_signal = rec->len_raw_signal;
                    fprintf(stdout,"%" PRIu64 "\n",len_raw_signal);

                    fast5_t *f5 = (fast5_t*)calloc(1, sizeof(fast5_t));
                    MALLOC_CHK(f5);

                    f5->nsample = len_raw_signal;
                    f5->rawptr = (float*)calloc(len_raw_signal, sizeof(float));
                    f5->digitisation = rec->digitisation;
                    f5->offset = rec->offset;
                    f5->range = rec->range;
                    f5->sample_rate = rec->sampling_rate;
                    db_t->f5[record_count] = f5;

                    for(uint64_t i=0;i<len_raw_signal;i++){
                        db_t->f5[record_count]->rawptr[i] = rec->raw_signal[i];
                    }
                }
                record_count++;
            }
        }

        process_db2(core,db_t);

        free_db_tmp2(db_t);
        free_db2(db_t);

        if(flag_EOF){
            break;
        }
    }

//    printf("return value: %d\n", ret_kseq_read);
    kseq_destroy(seq);
    gzclose(fp);
//
    slow5_rec_free(rec);
    slow5_idx_unload(sp);
    slow5_close(sp);

    free_core2(core);

    return 0;

}




























