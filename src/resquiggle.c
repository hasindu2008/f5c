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

void pthread_db(core_t* core, db_t* db, void (*func)(core_t*,db_t*,int));
void event_single(core_t* core,db_t* db, int32_t i);
void scaling_single(core_t* core, db_t* db, int32_t i);

//todo: CUDA support


static const char *RESQUIGGLE_USAGE_MESSAGE =
    "Usage: f5c resquiggle [OPTIONS] reads.fastq signals.blow5\n"
    "Align raw signals to basecalled reads.\n\n"
    "Options:\n"
    "   -t INT              number of processing threads\n"
    "   -K INT              batch size (max number of reads loaded at once)\n"
    "   -o FILE             output file. Write to stdout if not specified\n"
    "   --rna               the dataset is direct RNA\n"
    "   -h                  help\n"
    "   --version           print version\n"
    "   --verbose INT       verbosity level\n"
    "   --kmer-model FILE   custom nucleotide k-mer model file (format similar to test/r9-models/r9.4_450bps.nucleotide.6mer.template.model)\n"
    "\n\n"
    ;


static struct option long_options[] = {
    {"threads", required_argument, 0, 't'},        //0 number of threads [8]
    {"batchsize", required_argument, 0, 'K'},      //1 batchsize - number of reads loaded at once [512]
    {"verbose", required_argument, 0, 'v'},        //2 verbosity level [1]
    {"help", no_argument, 0, 'h'},                 //3
    {"version", no_argument, 0, 'V'},              //4
    {"kmer-model", required_argument, 0, 0},       //5 custom nucleotide k-mer model file
    {"output",required_argument, 0, 'o'},          //6 output to a file [stdout]
    {"rna",no_argument,0,0},                       //7 if RNA
    {0, 0, 0, 0}
};

/* initialise the core data structure */
core_t* init_core_rsq(opt_t opt, const char *slow5file) {

    core_t* core = (core_t*)malloc(sizeof(core_t));
    MALLOC_CHK(core);

    core->opt = opt;

    //model
    core->model = (model_t*)malloc(sizeof(model_t) * MAX_NUM_KMER); //4096 is 4^6 which is hardcoded now
    MALLOC_CHK(core->model);

    //load the model from files
    uint32_t kmer_size=0;
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
    core->kmer_size = kmer_size;

    core->sf = slow5_open(slow5file,"r");
    if (core->sf == NULL) {
        STDERR("Error opening SLOW5 file %s\n",slow5file);
        exit(EXIT_FAILURE);
    }
    int ret=slow5_idx_load(core->sf);
    if(ret<0){
        STDERR("Error in loading SLOW5 index for %s\n",slow5file);
        exit(EXIT_FAILURE);
    }

    return core;
}

void free_core_rsq(core_t* core) {
    free(core->model);
    slow5_idx_unload(core->sf);
    slow5_close(core->sf);
    free(core);

}

/* initialise a data batch */
db_t* init_db_rsq(core_t* core) {
    db_t* db = (db_t*)(malloc(sizeof(db_t)));
    MALLOC_CHK(db);

    db->capacity_bam_rec = core->opt.batch_size;
    db->n_bam_rec = 0;

    db->read = (char**)(malloc(sizeof(char*) * db->capacity_bam_rec));
    MALLOC_CHK(db->read);
    db->read_id = (char**)(malloc(sizeof(char*) * db->capacity_bam_rec));
    MALLOC_CHK(db->read_id);
    db->read_len = (int32_t*)(malloc(sizeof(int32_t) * db->capacity_bam_rec));
    MALLOC_CHK(db->read_len);

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

    db->read_stat_flag = (int32_t *)calloc(db->capacity_bam_rec,sizeof(int32_t));
    MALLOC_CHK(db->read_stat_flag);

    db->total_reads=0;
    db->bad_fast5_file=0;
    //db->ultra_long_skipped=0;

    return db;
}

void free_db_rsq(db_t* db) {

    free(db->read);
    free(db->read_id);
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

    free(db);
}

/* partially free a data batch - only the read dependent allocations are freed */
void free_db_tmp_rsq(db_t* db) {
    int32_t i = 0;
    for (i = 0; i < db->n_bam_rec; ++i) {
        free(db->read[i]);
        free(db->read_id[i]);
        free(db->f5[i]->rawptr);
        free(db->f5[i]);
        free(db->et[i].event);
        free(db->event_align_pairs[i]);
        free(db->base_to_event_map[i]);
    }
}

void process_single_rsq(core_t* core, db_t* db,int32_t i) {
    event_single(core,db,i);
    align_single(core, db,i);
    scaling_single(core,db,i);
}

void process_db_rsq(core_t* core, db_t* db) {
    if (core->opt.num_thread == 1) {
        int32_t i=0;
        for (i = 0; i < db->capacity_bam_rec; i++) {
            process_single_rsq(core,db,i);
        }
    }
    else {
        pthread_db(core,db,process_single_rsq);
    }
}

void output_db_rsq(core_t* core, db_t* db) {
    for (int i = 0; i < db->n_bam_rec; i++) {
        if((db->read_stat_flag[i]) & FAILED_ALIGNMENT){
            continue;
        }
        event_table et = db->et[i];

        index_pair_t *indexPair = db->base_to_event_map[i];
        int32_t n_kmers = db->read_len[i] - core->kmer_size + 1;
        int32_t prev_st_event_idx = -1; //prev_valid_start_event_idx
        //nt merged_events = 0;    //multiple_bases_per_event_count

        for (int32_t j=0; j<n_kmers; j++){
            int32_t start_event_idx = indexPair[j].start;
            int end_event_idx = indexPair[j].stop;
            if(start_event_idx == -1 && prev_st_event_idx != -1){
                //merged_events++;
                start_event_idx = prev_st_event_idx;
                if(end_event_idx != -1){
                    fprintf(stderr,"assertion failed\n");
                    exit(EXIT_FAILURE);
                }
                end_event_idx = prev_st_event_idx;
            }else{
                prev_st_event_idx = start_event_idx;
            }
            int signal_start_point = et.event[start_event_idx].start;
            int signal_end_point = et.event[end_event_idx].start + (int)et.event[end_event_idx].length;
            printf("%s\t%d\t%d\t%d\t%d\t%d\n", db->read_id[i], j, start_event_idx, end_event_idx, signal_start_point, signal_end_point);
        }
        //fprintf(stderr,"merged_events=%d\n", merged_events);
    }
}

ret_status_t load_db_rsq(core_t* core, db_t* db, gzFile fp) {

    ret_status_t status={0,0};

    kseq_t *seq;
    int ret_kseq_read;
    seq = kseq_init(fp);
    if(seq == NULL){
        ERROR("%s","Failed to init kseq");
        exit(EXIT_FAILURE);
    }

    int64_t i = 0;
    int ret = 0;
    slow5_rec_t *rec = NULL;

    while (i < core->opt.batch_size && status.num_bases<core->opt.batch_size_bases) {
        if ((ret_kseq_read = kseq_read(seq)) < 0) {
            if (ret_kseq_read < -1) {
                ERROR("Failed to read kseq. Error code is %d", ret_kseq_read);
                exit(EXIT_FAILURE);
            } else { //EOF file reached
                break;
            }
        } else {
            db->read[i] = strdup(seq->seq.s);
            NULL_CHK(db->read[i]);
            db->read_id[i] = strdup(seq->name.s);
            NULL_CHK(db->read_id[i]);
            db->read_len[i] = strlen(db->read[i]);
            status.num_bases += db->read_len[i];
            if(core->opt.flag & F5C_RNA){
                replace_char(db->read[i], 'U', 'T');
            }
            ret = slow5_get(seq->name.s, &rec, core->sf);
            if(ret < 0){
                fprintf(stderr,"Error in when fetching the read\n");
            }
            else{
                uint64_t len_raw_signal = rec->len_raw_signal;
                fast5_t *f5 = (fast5_t*)calloc(1, sizeof(fast5_t));
                MALLOC_CHK(f5);

                f5->nsample = len_raw_signal;
                f5->rawptr = (float*)calloc(len_raw_signal, sizeof(float));
                MALLOC_CHK(f5->rawptr);
                f5->digitisation = rec->digitisation;
                f5->offset = rec->offset;
                f5->range = rec->range;
                f5->sample_rate = rec->sampling_rate;
                db->f5[i] = f5;

                for(uint64_t j=0;j<len_raw_signal;j++){
                    db->f5[i]->rawptr[j] = rec->raw_signal[j];
                }
            }
            i++;
        }
    }
    slow5_rec_free(rec);
    kseq_destroy(seq);
    db->n_bam_rec = i;
    status.num_reads = i;

    return status;
}

int resquiggle_main(int argc, char **argv) {

    double realtime0 = realtime();


    FILE *fp_help = stderr;
    const char* optstring = "t:K:v:o:hV";
    int longindex = 0;
    int32_t c = -1;
    char *slow5file = NULL;

    opt_t opt;
    init_opt(&opt); //initialise options to defaults

    //parse the user args
    while ((c = getopt_long(argc, argv, optstring, long_options, &longindex)) >= 0) {
        if (c == 'K') {
            opt.batch_size = atoi(optarg);
            if (opt.batch_size < 1) {
                ERROR("Batch size should larger than 0. You entered %d",opt.batch_size);
                exit(EXIT_FAILURE);
            }
        } else if (c == 't') {
            opt.num_thread = atoi(optarg);
            if (opt.num_thread < 1) {
                ERROR("Number of threads should larger than 0. You entered %d", opt.num_thread);
                exit(EXIT_FAILURE);
            }
        }
        else if (c=='v'){
            opt.verbosity = atoi(optarg);
        }
        else if (c=='V'){
            fprintf(stdout,"F5C %s\n",F5C_VERSION);
            exit(EXIT_SUCCESS);
        }
        else if (c=='h'){
            fp_help = stdout;
        }
        else if(c=='o'){
			if (strcmp(optarg, "-") != 0) {
				if (freopen(optarg, "wb", stdout) == NULL) {
					ERROR("failed to write the output to file %s : %s",optarg, strerror(errno));
					exit(EXIT_FAILURE);
				}
			}
        } else if (c == 0 && longindex == 5) { //custom nucleotide model file
            opt.model_file = optarg;
        }
        else if (c == 0 && longindex == 7) {
            opt.flag |= F5C_RNA;
        }
    }

    // No arguments given
    if (argc - optind != 2 || fp_help == stdout) {
        fprintf(fp_help, RESQUIGGLE_USAGE_MESSAGE, argv[0]);
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    // open fastq
    gzFile fp;

    fp = gzopen(argv[optind], "r");
    if(fp == NULL){
        ERROR("Failed to open reads file %s : %s",argv[optind], strerror(errno));
        exit(EXIT_FAILURE);
    }

    WARNING("%s","f5c resquiggle is experimental. Use with caution. Report any bugs under GitHub issues.");

    //open slow5
    slow5file = argv[optind+1];

    ret_status_t status = {opt.batch_size,opt.batch_size_bases};

    char resquiggle_header0 [] = "read_id\tkmer_idx\tstart_event_idx\tend_event_idx\tstart_raw_idx\tend_signal_idx";
    fprintf(stdout, "%s\n", resquiggle_header0);

    //initialise the core data structure
    core_t* core = init_core_rsq(opt, slow5file);

    while (status.num_reads >= core->opt.batch_size || status.num_bases>=core->opt.batch_size_bases){
        db_t* db = init_db_rsq(core);
        status = load_db_rsq(core,db,fp);

        fprintf(stderr, "[%s::%.3f*%.2f] %d Entries (%.1fM bases) loaded\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                status.num_reads,status.num_bases/(1000.0*1000.0));

        process_db_rsq(core,db);

        fprintf(stderr, "[%s::%.3f*%.2f] %d Entries (%.1fM bases) processed\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                status.num_reads,status.num_bases/(1000.0*1000.0));

        output_db_rsq(core, db);
        free_db_tmp_rsq(db);
        free_db_rsq(db);
    }

    //these must be updated and fixed
    fprintf(stderr, "[%s] total entries: %ld, qc fail: %ld, could not calibrate: %ld, no alignment: %ld, bad read: %ld",
             __func__,(long)core->total_reads, (long)core->qc_fail_reads, (long)core->failed_calibration_reads, (long)core->failed_alignment_reads, (long)core->bad_fast5_file);
    fprintf(stderr,"\n[%s] total bases: %.1f Mbases",__func__,core->sum_bases/(float)(1000*1000));

    fprintf(stderr, "\n[%s] Data loading time: %.3f sec", __func__,core->load_db_time);
    fprintf(stderr, "\n[%s]     - fasta load time: %.3f sec", __func__, core->db_fasta_time);
    fprintf(stderr, "\n[%s]     - slow5 load time: %.3f sec", __func__, core->db_fast5_time);
    fprintf(stderr, "\n[%s] Data processing time: %.3f sec", __func__,core->process_db_time);
    fprintf(stderr, "\n[%s] Data output time: %.3f sec", __func__,core->output_time);

    fprintf(stderr,"\n");

    gzclose(fp);
    free_core_rsq(core);

    return 0;

}