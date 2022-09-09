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
#include "str.h"

#define TO_PICOAMPS(RAW_VAL,DIGITISATION,OFFSET,RANGE) (((RAW_VAL)+(OFFSET))*((RANGE)/(DIGITISATION)))

KSEQ_INIT(gzFile, gzread)

//TODO add to f5c.h and document
void pthread_db(core_t* core, db_t* db, void (*func)(core_t*,db_t*,int));
void event_single(core_t* core,db_t* db, int32_t i);
void scaling_single(core_t* core, db_t* db, int32_t i);

//todo: CUDA support
//TODO there are lot of copy pasted sections from f5c.c - can be modularised

static inline int8_t drna_detect(slow5_file_t *sp){

    const slow5_hdr_t* hdr = sp->header;
    int8_t rna = 0;
    char *exp =slow5_hdr_get("experiment_type", 0, hdr);
    if(exp==NULL){
        WARNING("%s","experiment_type not found in SLOW5 header. Assuming genomic_dna");
        return 0;
    }
    if (strcmp(exp,"genomic_dna")==0){
        rna = 0;
    }else if (strcmp(exp,"rna")==0){
        rna = 1;
    } else {
        WARNING("Unknown experiment type: %s. Assuming genomic_dna", exp);
    }

    for(uint32_t  i=1; i < hdr->num_read_groups; i++){
        char *curr =slow5_hdr_get("experiment_type", i, hdr);
        if (strcmp(curr, exp)){
            WARNING("Experiment type mismatch: %s != %s in read group %d. Defaulted to %s", curr, exp, i, exp);
        }
    }
    return rna;
}


static void print_help_msg(FILE *fp_help, opt_t opt){
    fprintf(fp_help,"Usage: f5c [OPTIONS] reads.fastq signals.blow5\n");
    fprintf(fp_help,"\noptions:\n");
    fprintf(fp_help,"   -t INT                     number of processing threads [%d]\n",opt.num_thread);
    fprintf(fp_help,"   -K INT                     batch size (max number of reads loaded at once) [%d]\n",opt.batch_size);
    fprintf(fp_help,"   -B FLOAT[K/M/G]            max number of bases loaded at once [%.1fM]\n",opt.batch_size_bases/(float)(1000*1000));
    fprintf(fp_help,"   -h                         help\n");
    fprintf(fp_help,"   -o FILE                    output to file [stdout]\n");
    fprintf(fp_help,"   -x STR                     parameter profile to be used for better performance (always applied before other options)\n"); //Added option in help
    fprintf(fp_help,"                              e.g., laptop, desktop, hpc; see https://f5c.page.link/profiles for the full list\n");
    fprintf(fp_help,"   -c                         print in paf format\n");
    fprintf(fp_help,"   --verbose INT              verbosity level [%d]\n",opt.verbosity);
    fprintf(fp_help,"   --version                  print version\n");
    fprintf(fp_help,"   --kmer-model FILE          custom nucleotide k-mer model file (format similar to test/r9-models/r9.4_450bps.nucleotide.6mer.template.model)\n");
    fprintf(fp_help,"   --rna                      the dataset is direct RNA\n");
#ifdef HAVE_CUDA
        fprintf(fp_help,"   --disable-cuda=yes|no      disable running on CUDA [%s]\n",(opt.flag&F5C_DISABLE_CUDA?"yes":"no"));
        fprintf(fp_help,"   --cuda-dev-id INT          CUDA device ID to run kernels on [%d]\n",opt.cuda_dev_id);
        fprintf(fp_help,"   --cuda-mem-frac FLOAT      Fraction of free GPU memory to allocate [0.9 (0.7 for tegra)]\n");
#endif

}

static struct option long_options[] = {
    {"threads", required_argument, 0, 't'},        //0 number of threads [8]
    {"batchsize", required_argument, 0, 'K'},      //1 batchsize - number of reads loaded at once [512]
    {"verbose", required_argument, 0, 'v'},        //2 verbosity level [1]
    {"help", no_argument, 0, 'h'},                 //3
    {"version", no_argument, 0, 'V'},              //4
    {"kmer-model", required_argument, 0, 0},       //5 custom nucleotide k-mer model file
    {"output",required_argument, 0, 'o'},          //6 output to a file [stdout]
    {"rna",no_argument,0,0},                       //7 if RNA
    {"max-bases", required_argument, 0, 'B'},      //8 batchsize - number of bases loaded at once
    {"profile",required_argument, 0,'x'},          //9 profile used to tune parameters for GPU
    {"disable-cuda", required_argument, 0, 0},     //10 disable running on CUDA [no] (only if compiled for CUDA)
    {"cuda-dev-id",required_argument, 0, 0},       //11 cuda device ID to run on (only if compiled for CUDA)
    {"cuda-mem-frac",required_argument, 0, 0},     //12 fraction of the free GPU memory to use (only if compiled for CUDA)
    {"paf",no_argument,0,'c'},                       //13 if print in paf format
    {0, 0, 0, 0}
};

/* initialise the core data structure */
core_t* init_core_rsq(opt_t opt, const char *slow5file, double realtime0) {

    core_t* core = (core_t*)malloc(sizeof(core_t));
    MALLOC_CHK(core);

    core->opt = opt;

    //slow5 FILE
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
    if(drna_detect(core->sf)) {
        core->opt.flag |= F5C_RNA;
        opt.flag |= F5C_RNA;
    }

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

    //realtime0
    core->realtime0=realtime0;

    core->load_db_time=0;
    core->process_db_time=0;
    core->db_fasta_time=0;
    core->db_fast5_time=0;
    core->event_time=0;
    core->align_time=0;
    core->est_scale_time=0;
    core->output_time=0;

    //cuda stuff
#ifdef HAVE_CUDA
    if (!(core->opt.flag & F5C_DISABLE_CUDA)) {
        init_cuda(core);
    }
#endif

    core->sum_bases=0;
    core->total_reads=0; //total number mapped entries in the bam file (after filtering based on flags, mapq etc)
    core->bad_fast5_file=0; //empty fast5 path returned by readdb, could not open fast5
    core->ultra_long_skipped=0;
    core->qc_fail_reads=0;
    core->failed_calibration_reads=0;
    core->failed_alignment_reads=0;

    return core;
}

void free_core_rsq(core_t* core) {
    free(core->model);
    slow5_idx_unload(core->sf);
    slow5_close(core->sf);

#ifdef HAVE_CUDA
    if (!(core->opt.flag & F5C_DISABLE_CUDA)) {
        free_cuda(core);
    }
#endif

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

    db->sig = (signal_t**)malloc(sizeof(signal_t*) * db->capacity_bam_rec);
    MALLOC_CHK(db->sig);

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
    db->ultra_long_skipped=0;

    return db;
}

void free_db_rsq(db_t* db) {

    free(db->read);
    free(db->read_id);
    free(db->read_len);
    free(db->et);
    free(db->sig);
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
        free(db->sig[i]->rawptr);
        free(db->sig[i]);
        free(db->et[i].event);
        free(db->event_align_pairs[i]);
        free(db->base_to_event_map[i]);
    }
}


//high similarity with f5c.c - can be modularised
void process_db_rsq(core_t* core, db_t* db) {

    double process_start = realtime();
    double realtime0=core->realtime0;

    double event_start = realtime();
    pthread_db(core,db,event_single);
    double event_end = realtime();
    core->event_time += (event_end-event_start);

    fprintf(stderr, "[%s::%.3f*%.2f] Events computed\n", __func__,
            realtime() - realtime0, cputime() / (realtime() - realtime0));

    double align_start = realtime();
    align_db(core, db);
    double align_end = realtime();
    core->align_time += (align_end-align_start);

    fprintf(stderr, "[%s::%.3f*%.2f] Banded alignment done\n", __func__,
            realtime() - realtime0, cputime() / (realtime() - realtime0));

    double est_scale_start = realtime();
    pthread_db(core,db,scaling_single);
    double est_scale_end = realtime();
    core->est_scale_time += (est_scale_end-est_scale_start);

    fprintf(stderr, "[%s::%.3f*%.2f] Scaling calibration done\n", __func__,
            realtime() - realtime0, cputime() / (realtime() - realtime0));

    double process_end= realtime();
    core->process_db_time += (process_end-process_start);

}

void output_db_rsq(core_t* core, db_t* db, int8_t fmt) {

    double output_start = realtime();

    for (int i = 0; i < db->n_bam_rec; i++) {
        if(!db->read_stat_flag[i]){
            event_table et = db->et[i];

            index_pair_t *base_to_event_map = db->base_to_event_map[i];
            int32_t n_kmers = db->read_len[i] - core->kmer_size + 1;
            int64_t signal_start_point = -1;
            int64_t signal_start_point2 = -1; //very first, to be printed later
            int64_t signal_end_point = -1;
            int64_t signal_end_point2 = -1;
            int64_t read_start = -1;
            int64_t read_end = -1;

            kstring_t str;
            kstring_t *sp = &str;
            str_init(sp, db->read_len[i]*3);
            int64_t ci = 0; //current index
            int64_t mi = 0;
            int64_t d = 0; //deletion count
            int8_t ff = 1; //first start flag
            int matches = 0;
            int8_t rna = core->opt.flag & F5C_RNA;

            if (rna){
                for (int j = 0; j < n_kmers/2; ++j) {
                    index_pair_t tmp= base_to_event_map[j];
                    base_to_event_map[j]=base_to_event_map[n_kmers-1-j];
                    base_to_event_map[n_kmers-1-j]=tmp;
                }
                for (int j = 0; j < n_kmers; ++j) {
                    int32_t tmp = base_to_event_map[j].start;
                    base_to_event_map[j].start = base_to_event_map[j].stop;
                    base_to_event_map[j].stop = tmp;
                }
            }

            for (int j=0;j<n_kmers; j++){
                int32_t start_event_idx = base_to_event_map[j].start;
                int32_t end_event_idx = base_to_event_map[j].stop;
                if(start_event_idx == -1){ //deletion from ref
                    assert(end_event_idx == -1);
                    signal_start_point = signal_end_point = -1;
                    if(!ff){
                        assert(j!=0);
                        d++;
                    }

                } else {
                    assert(end_event_idx != -1);
                    //assert(start_event_idx <= end_event_idx);

                    signal_start_point = et.event[start_event_idx].start; //inclusive
                    if(ff) {
                        signal_start_point2 = signal_start_point;
                        read_start = j;
                        ff=0;
                    }
                    signal_end_point2=signal_end_point = et.event[end_event_idx].start + (int)et.event[end_event_idx].length; //non-inclusive
                    read_end = j;
                    if(fmt){
                        if(d>0){
                            sprintf_append(sp,"%dD",d);
                            d=0;
                        }
                        if(j==0) ci = signal_start_point;
                        ci += (mi =  signal_start_point - ci);
                        if(mi) sprintf_append(sp,"%dI",(int)mi);
                        ci += (mi = signal_end_point-signal_start_point);

                        if(mi) {
                            matches++;
                            sprintf_append(sp,"%d,",(int)mi);
                        }
                    }

                }

                if(fmt==0){
                    printf("%s\t%d\t",db->read_id[i], rna? n_kmers-j-1 : j);
                    if(signal_start_point<0){
                        printf(".\t");
                    } else {
                        printf("%ld\t", (long)signal_start_point);
                    }
                    if(signal_end_point<0){
                        printf(".\t");
                    } else {
                        printf("%ld", (long)signal_end_point);
                    }
                    printf("\n");
                    if(signal_start_point>=0 && signal_end_point>=0){
                        if(signal_end_point<=signal_start_point){
                            fprintf(stderr,"ERROR: signal_end_point(%ld)<=signal_start_point(%ld)\n",signal_end_point,signal_start_point);
                            fprintf(stderr,"read_id=%s\tj=%d\tstartevent=%d\tendevent=%d\n",db->read_id[i],j,start_event_idx,end_event_idx);
                            exit(EXIT_FAILURE);
                        }
                    }

                }

            }

            if(fmt==1){

                assert(n_kmers>0);
                str.s[str.l] = '\0';

                assert(signal_start_point2!=-1);
                assert(signal_end_point2!=-1);
                //query: name, start, end, strand
                printf("%s\t%ld\t%ld\t%ld\t+\t", db->read_id[i], (long)db->sig[i]->nsample,(long)signal_start_point2, (long)signal_end_point2);
                //target: name, start, end
                printf("%s\t%d\t%ld\t%ld\t", db->read_id[i], n_kmers,(long)(rna?n_kmers-read_start:read_start), (long)(rna?n_kmers-1-read_end:read_end+1));
                //residue matches, block len, mapq
                printf("%d\t%d\t%d\t",matches,n_kmers,255);
                printf("sc:f:%f\t",db->scalings->scale);
                printf("sh:f:%f\t",db->scalings->shift);
                printf("ss:Z:%s\n",str.s);

            }
            str_free(sp);

            //fprintf(stderr,"merged_events=%d\n", merged_events);
        } else {  //common section with f5c.c - can be modularised
            if((db->read_stat_flag[i])&FAILED_CALIBRATION){
                core->failed_calibration_reads++;
            }
            else if ((db->read_stat_flag[i])&FAILED_ALIGNMENT){
                core->failed_alignment_reads++;
            }
            else if ((db->read_stat_flag[i])&FAILED_QUALITY_CHK){
                core->qc_fail_reads++;
            }
            else{
                assert(0);
            }
        }
    }

    //common section with f5c.c - can be modularised
    core->sum_bases += db->sum_bases;
    core->total_reads += db->total_reads;
    core->bad_fast5_file += db->bad_fast5_file;
    core->ultra_long_skipped += db->ultra_long_skipped;

    double output_end = realtime();
    core->output_time += (output_end-output_start);

}

static void read_slow5_single(core_t* core, db_t* db, int i){

    slow5_rec_t *record=NULL;
    int len=0;
    len = slow5_get(db->read_id[i], &record, core->sf);

    db->sig[i] = (signal_t*)calloc(1, sizeof(signal_t));
    MALLOC_CHK(db->sig[i]);

    if(record==NULL || len <0){ //todo : should we free if len<0
        db->bad_fast5_file++;
        if (core->opt.flag & F5C_SKIP_UNREADABLE) {
            WARNING("Slow5 record for read [%s] is unavailable/unreadable and will be skipped", db->read_id[i]);
            db->sig[i]->nsample = 0;
            db->sig[i]->rawptr = NULL;
        } else {
            ERROR("Slow5 record for read [%s] is unavailable/unreadable", db->read_id[i]);
            exit(EXIT_FAILURE);
        }
    }
    else{
        assert(strcmp(db->read_id[i],record->read_id)==0);
        db->sig[i]->nsample = record->len_raw_signal;  //n_samples
        assert(db->sig[i]->nsample>0);
        db->sig[i]->rawptr = (float*)calloc(db->sig[i]->nsample, sizeof(float));
        MALLOC_CHK( db->sig[i]->rawptr);

        db->sig[i]->digitisation = (float)record->digitisation;
        db->sig[i]->offset = (float)record->offset;
        db->sig[i]->range = (float)record->range;
        db->sig[i]->sample_rate = (float)record->sampling_rate;

        for (int j = 0; j < (int)db->sig[i]->nsample; j++) { //check for int overflow
            db->sig[i]->rawptr[j] = (float)record->raw_signal[j];
        }
        // ret=1;
        slow5_rec_free(record);
    }

}

ret_status_t load_db_rsq(core_t* core, db_t* db, gzFile fp, kseq_t *seq) {

    double load_start = realtime();

    ret_status_t status={0,0};
    db->n_bam_rec = 0;
    db->sum_bases = 0;
    db->total_reads = 0;
    db->bad_fast5_file = 0;
    db->ultra_long_skipped =0;

    double t = realtime();
    int64_t i = 0;
    int ret_kseq_read;

    while (i < core->opt.batch_size && status.num_bases<core->opt.batch_size_bases) {
        if ((ret_kseq_read = kseq_read(seq)) < 0) {
            if (ret_kseq_read < -1) {
                ERROR("Failed to read kseq. Error code is %d", ret_kseq_read);
                exit(EXIT_FAILURE);
            } else { //EOF file reached
                break;
            }
        } else {
            //fprintf(stderr,"Loaded: %s\n", seq->name.s);
            db->read[i] = strdup(seq->seq.s);
            NULL_CHK(db->read[i]);
            db->read_id[i] = strdup(seq->name.s);
            NULL_CHK(db->read_id[i]);
            db->read_len[i] = strlen(db->read[i]);
            status.num_bases += db->read_len[i];
            if(core->opt.flag & F5C_RNA){
                replace_char(db->read[i], 'U', 'T');
            }
            db->total_reads++;
            db->sum_bases += db->read_len[i];
            db->read_stat_flag[i] = 0; //reset the flag

            i++;
        }
    }

    core->db_fasta_time += realtime() - t;

    db->n_bam_rec = i;
    status.num_reads = i;

    t=realtime();
    pthread_db(core, db, read_slow5_single);
    core->db_fast5_time += realtime() - t;

    double load_end = realtime();
    core->load_db_time += (load_end-load_start);

    return status;
}

//parse yes or no arguments : taken from minimap2
static inline void yes_or_no(opt_t* opt, uint64_t flag, int long_idx,
                             const char* arg,
                             int yes_to_set)
{
    if (yes_to_set) {
        if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) {
            opt->flag |= flag;
        } else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) {
            opt->flag &= ~flag;
        } else {
            WARNING("option '--%s' only accepts 'yes' or 'no'.",
                    long_options[long_idx].name);
        }
    } else {
        if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) {
            opt->flag &= ~flag;
        } else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) {
            opt->flag |= flag;
        } else {
            WARNING("option '--%s' only accepts 'yes' or 'no'.",
                    long_options[long_idx].name);
        }
    }
}

int resquiggle_main(int argc, char **argv) {

    double realtime0 = realtime();


    FILE *fp_help = stderr;
    const char* optstring = "t:K:B:v:o:hVx:c";
    int longindex = 0;
    int32_t c = -1;
    char *slow5file = NULL;
    char* profilename = NULL; //Create variable to store profile arg
    int8_t fmt = 0;

    opt_t opt;
    init_opt(&opt); //initialise options to defaults

    //parse the user args
    while ((c = getopt_long(argc, argv, optstring, long_options, &longindex)) >= 0) {
        if (c == 'x') {  //Set profile values. Any user-specified arguments will override the profile values.
            profilename = optarg;
            int profile_return = set_profile(profilename,&opt);
            if(profile_return == 1){
                ERROR("%s","Error occurred setting the profile.");
                exit(EXIT_FAILURE);
            }
        } else if (c == 'K') {
            opt.batch_size = atoi(optarg);
            if (opt.batch_size < 1) {
                ERROR("Batch size should larger than 0. You entered %d",opt.batch_size);
                exit(EXIT_FAILURE);
            }
        } else if (c == 'B') {
            opt.batch_size_bases = mm_parse_num(optarg);
            if(opt.batch_size_bases<=0){
                ERROR("%s","Maximum number of bases should be larger than 0.");
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
        } else if (c == 'c'){
            fmt=1;
        } else if (c == 0 && longindex == 5) { //custom nucleotide model file
            opt.model_file = optarg;
        }
        else if (c == 0 && longindex == 7) { //if RNA
            opt.flag |= F5C_RNA;
        } else if (c == 0 && longindex == 10) {
#ifdef HAVE_CUDA
            yes_or_no(&opt, F5C_DISABLE_CUDA, longindex, optarg, 1);
#else
            WARNING("%s", "disable-cuda has no effect when compiled for the CPU");
#endif
        } else if(c == 0 && longindex == 11){ //todo : warning for CPU mode
            opt.cuda_dev_id = atoi(optarg);
        } else if(c == 0 && longindex == 12){ //todo : warning for CPU mode, warning for dynamic malloc mode
            opt.cuda_mem_frac = atof(optarg);
        }
    }

    // No arguments given
    if (argc - optind != 2 || fp_help == stdout) {
        print_help_msg(fp_help, opt);
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
    kseq_t *seq;
    seq = kseq_init(fp);
    if(seq == NULL){
        ERROR("%s","Failed to init kseq");
        exit(EXIT_FAILURE);
    }

    WARNING("%s","f5c resquiggle is experimental. Use with caution. Report any bugs under GitHub issues.");

    //open slow5
    slow5file = argv[optind+1];

    ret_status_t status = {opt.batch_size,opt.batch_size_bases};

    if(fmt==0){
        const char *resquiggle_header0 = "read_id\tkmer_idx\tstart_raw_idx\tend_raw_idx";
        fprintf(stdout, "%s\n", resquiggle_header0);
    }

    //initialise the core data structure
    core_t* core = init_core_rsq(opt, slow5file, realtime0); //use opt inside core afte this call

    while (status.num_reads >= core->opt.batch_size || status.num_bases>=core->opt.batch_size_bases){
        db_t* db = init_db_rsq(core);
        status = load_db_rsq(core,db,fp,seq);

        fprintf(stderr, "[%s::%.3f*%.2f] %d Entries (%.1fM bases) loaded\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                status.num_reads,status.num_bases/(1000.0*1000.0));

        process_db_rsq(core,db);

        fprintf(stderr, "[%s::%.3f*%.2f] %d Entries (%.1fM bases) processed\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                status.num_reads,status.num_bases/(1000.0*1000.0));

        output_db_rsq(core, db, fmt);
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
    fprintf(stderr, "\n[%s]     - Events time: %.3f sec", __func__, core->event_time);
    fprintf(stderr, "\n[%s]     - Alignment time: %.3f sec", __func__, core->align_time);
    #ifdef HAVE_CUDA
        if (!(core->opt.flag & F5C_DISABLE_CUDA)) {
            fprintf(stderr, "\n[%s]           -cpu preprocess time: %.3f sec",
                __func__, core->align_cuda_preprocess);
            fprintf(stderr, "\n[%s]           -cuda malloc time: %.3f sec",
                __func__, core->align_cuda_malloc);
            fprintf(stderr, "\n[%s]           -cuda data transfer time: %.3f sec",
                __func__, core->align_cuda_memcpy);
            fprintf(stderr, "\n[%s]           -cuda kernel time: %.3f sec",
                __func__, core->align_kernel_time);
            fprintf(stderr, "\n[%s]                -align-pre kernel only time: %.3f sec",
                __func__, core->align_pre_kernel_time);
            fprintf(stderr, "\n[%s]                -align-core kernel only time: %.3f sec",
                __func__, core->align_core_kernel_time);
            fprintf(stderr, "\n[%s]                -align-post kernel only time: %.3f sec",
                __func__, core->align_post_kernel_time);

            fprintf(stderr, "\n[%s]           -cpu postprocess time: %.3f sec",
                __func__, core->align_cuda_postprocess);
            fprintf(stderr, "\n[%s]           -additional cpu processing time (load imbalance): %.3f sec",
                __func__, core->extra_load_cpu);
        }
    #endif
    fprintf(stderr, "\n[%s]     - Estimate scaling time: %.3f sec", __func__, core->est_scale_time);
    fprintf(stderr, "\n[%s] Data output time: %.3f sec", __func__,core->output_time);

    fprintf(stderr,"\n");
    kseq_destroy(seq);
    gzclose(fp);
    free_core_rsq(core);

    return 0;

}