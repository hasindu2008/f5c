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
    "Usage: f5c resquiggle [OPTIONS] reads.fastq signals.blow5\n"
    "Align raw signal to the basecalled read.\n\n"
    "   -o FILE         output file. Write to stdout if not specified\n"
    "   -h              help\n"
    "   --version       print version\n"
    "\nSee the manual page for details (`man ./docs/f5c.1' or https://f5c.page.link/man)."
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
    db->read_id = (char**)(malloc(sizeof(char*) * db->capacity_bam_rec));
    MALLOC_CHK(db->read_id);
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
    for (i = 0; i < db->capacity_bam_rec; ++i) {
        db->read_stat_flag[i] = RESET_READ_STATUS;
    }

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
    for (i = 0; i < db->n_bam_rec; ++i) {
        free(db->read[i]);
        free(db->read_id[i]);
        free(db->f5[i]->rawptr);
        free(db->f5[i]);
        free(db->et[i].event);
        free(db->event_align_pairs[i]);
        free(db->base_to_event_map[i]);
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

void calculate_dwell_time_single(core_t* core,db_t* db, int32_t i){
    index_pair_t * indexPair = db->base_to_event_map[i];
    int32_t n_kmers = db->read_len[i] - core->kmer_size + 1;
    int32_t prev_valid_start_event_index = -1;
    int32_t prev_valid_j = -1;
    event_table et = db->et[i];

    int multiple_bases_per_event_count = 0;
    for (int32_t j=0; j<n_kmers; j++){
        int32_t start_event_index = indexPair[j].start;
        int end_event_index = indexPair[j].stop;
        int divider = 1;
        if(start_event_index == -1 && prev_valid_start_event_index != -1){
            multiple_bases_per_event_count++;
//            printf("start_event_index == -1 at basecalled_read_index %d\n", j);
            start_event_index = prev_valid_start_event_index;
            if(end_event_index != -1){
                //                assert(end_event_index == -1);
                fprintf(stderr,"assertion failed\n");
                exit(EXIT_FAILURE);
            }
            end_event_index = prev_valid_start_event_index;
            divider = j-prev_valid_j+1;
        }else{
            prev_valid_start_event_index = start_event_index;
            prev_valid_j = j;
        }
        int signal_start_point = et.event[start_event_index].start;
        int signal_end_point = et.event[end_event_index].start + (int)et.event[end_event_index].length;
        int dwell_time = signal_end_point - signal_start_point;
        for(int32_t k=prev_valid_j; k<=j; k++){
            indexPair[k].dwell_time = dwell_time/divider;
        }
    }
}

void process_single2(core_t* core, db_t* db,int32_t i) {
    event_single(core,db,i);
    align_single(core, db,i);
    scaling_single(core,db,i);
    //calculate_dwell_time_single(core,db,i);
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

void output_db2(core_t* core, db_t* db) {
    for (int i = 0; i < db->n_bam_rec; i++) {
        if((db->read_stat_flag[i]) & FAILED_ALIGNMENT){
            continue;
        }
        event_table et = db->et[i];
//        AlignedPair* event_align_pairs = db->event_align_pairs[i];
//        for (int32_t j = 0; j < db->n_event_align_pairs[i]; j++) {
//            printf("%d\t%d\t%d\n", event_align_pairs[j].ref_pos, event_align_pairs[j].read_pos, (int)et.event[event_align_pairs[j].read_pos].length);
//        }
        index_pair_t * indexPair = db->base_to_event_map[i];
        int32_t n_kmers = db->read_len[i] - core->kmer_size + 1;
        int32_t prev_valid_start_event_index = -1;
        int multiple_bases_per_event_count = 0;
        for (int32_t j=0; j<n_kmers; j++){
            int32_t start_event_index = indexPair[j].start;
            int end_event_index = indexPair[j].stop;
            if(start_event_index == -1 && prev_valid_start_event_index != -1){
                multiple_bases_per_event_count++;
//                printf("start_event_index == -1 at basecalled_read_index %d\n", j);
                start_event_index = prev_valid_start_event_index;
                if(end_event_index != -1){
//                assert(end_event_index == -1);
                    fprintf(stderr,"assertion failed\n");
                    exit(EXIT_FAILURE);
                }
                end_event_index = prev_valid_start_event_index;
            }else{
                prev_valid_start_event_index = start_event_index;
            }
            int signal_start_point = et.event[start_event_index].start;
            int signal_end_point = et.event[end_event_index].start + (int)et.event[end_event_index].length;
            int dwell_time = signal_end_point - signal_start_point;
            float corrected_dwell_time = indexPair[j].dwell_time;
            //printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%f\n", db->read_id[i], j, start_event_index, end_event_index, signal_start_point, signal_end_point, dwell_time, corrected_dwell_time);
            printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%f\n", db->read_id[i], j, start_event_index, end_event_index, signal_start_point, signal_end_point, 0, 0);
        }
        fprintf(stderr,"multiple_bases_per_event_count=%d\n", multiple_bases_per_event_count);
    }
}

int resquiggle_main(int argc, char **argv) {

    double realtime0 = realtime();

    FILE *fp_help = stderr;
    const char* optstring = "t:K:v:o:hV";
    int longindex = 0;
    int32_t c = -1;

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
        }
        else if (c == 0 && longindex == 7) {
            fprintf(stderr,"RNA set\n");
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
    kseq_t *seq;
    int ret_kseq_read;
    if (argc == 1) {
        fprintf(stderr, "Usage: %s <in.fasta>\n", argv[0]);
        return 1;
    }
    fp = gzopen(argv[optind], "r");
    seq = kseq_init(fp);

    //open slow5
    slow5_file_t *sp = slow5_open(argv[optind+1],"r");
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
    char resquiggle_header0 [] = "read_id\tbasecalled_read_index\tstart_event_index\tend_event_index\tstart_signal_index\tend_signal_index\tdwell_time\tcorrected_dwell_time";
    fprintf(stdout, "%s\n", resquiggle_header0);

    int64_t batch_size = opt.batch_size;

    //initialise the core data structure
    core_t* core = init_core2(opt, 1);

    int flag_EOF = 0;
    while (1){
        db_t* db = init_db2(core);
        int64_t record_count = 0;

        while (record_count < batch_size) {
            if ((ret_kseq_read = kseq_read(seq)) < 0) {
                if (ret_kseq_read < -1) {
                    return EXIT_FAILURE;
                } else { //EOF file reached
                    flag_EOF = 1;
                    break;
                }
            } else {
                db->read[record_count] = strdup(seq->seq.s);
                db->read_id[record_count] = strdup(seq->name.s);
                db->read_len[record_count] = strlen(db->read[record_count]);
                if(core->opt.flag & F5C_RNA){
                    replace_char(db->read[record_count], 'U', 'T');
                }
                ret = slow5_get(seq->name.s, &rec, sp);
                fprintf(stderr,"%s loaded\n",seq->name.s);
                if(ret < 0){
                    fprintf(stderr,"Error in when fetching the read\n");
                }
                else{
                     //       printf("name: %s\n", seq->name.s);
                    //        if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
                    //        printf("seq: %s\n", seq->seq.s);
                    //        if (seq->qual.l) printf("qual: %s\n", seq->qual.s);

                    uint64_t len_raw_signal = rec->len_raw_signal;
                    //fprintf(stdout,"%" PRIu64 "\n",len_raw_signal);

                    fast5_t *f5 = (fast5_t*)calloc(1, sizeof(fast5_t));
                    MALLOC_CHK(f5);

                    f5->nsample = len_raw_signal;
                    f5->rawptr = (float*)calloc(len_raw_signal, sizeof(float));
                    f5->digitisation = rec->digitisation;
                    f5->offset = rec->offset;
                    f5->range = rec->range;
                    f5->sample_rate = rec->sampling_rate;
                    db->f5[record_count] = f5;

                    for(uint64_t i=0;i<len_raw_signal;i++){
                        db->f5[record_count]->rawptr[i] = rec->raw_signal[i];
                    }
                }
                record_count++;
            }
        }
        db->n_bam_rec = record_count;

        process_db2(core,db);
        output_db2(core, db);
        free_db_tmp2(db);
        free_db2(db);

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




























