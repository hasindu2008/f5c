/* @file meth_main.c
**
** f5c call-methylation/eventalign entry point implementation
** @author: Hasindu Gamaarachchi (hasindu@unsw.edu.au)
** @@
******************************************************************************/

#include "f5c.h"
#include "f5cmisc.h"
#include <assert.h>
#include <getopt.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "logsum.h"

/* Input/processing/output interleave framework :
unless IO_PROC_NO_INTERLEAVE is set input, processing and output are interleaved
main thread
1. allocates and loads a databatch
2. create `pthread_processor` thread which will perform the processing
(note that `pthread_processor` is the process-controller that will spawn user specified number of processing threads - see below)
3. create the `pthread_post_processor` thread that will print the output and free the databatch one the `pthread_processor` is done
4. allocates and load another databatch
5. wait till the previous `pthread_processor` is done and perform 2
6. wait till the previous `pthread_post_processor` is done and perform 3
7. goto 4 untill all the input is processed
*/

/* Multithreaded framework for processing


*/


/* CUDA acceleration


*/

/*
TODO :
--print-* options should take a filename
a step controller if one requires to perform only upto a certain step such as alignment, ignoring the reset
--the names of functions and variablrs starting with meth_ are confusing as it stands for both meth and event now. should be fixed later
*/

/*
LIMITATIONS :
Does not support multi strand reads (2D and 1D^2 reads) at the moment
*/

//fast logsum data structure
float flogsum_lookup[p7_LOGSUM_TBL]; //todo : get rid of global vars

static struct option long_options[] = {
    {"reads", required_argument, 0, 'r'},          //0 fastq/fasta read file
    {"bam", required_argument, 0, 'b'},            //1 sorted bam file
    {"genome", required_argument, 0, 'g'},         //2 reference genome
    {"threads", required_argument, 0, 't'},        //3 number of threads [8]
    {"batchsize", required_argument, 0, 'K'},      //4 batchsize - number of reads loaded at once [512]
    {"max-bases", required_argument, 0, 'B'},      //5 batchsize - number of bases loaded at once
    {"verbose", required_argument, 0, 'v'},        //6 verbosity level [1]
    {"help", no_argument, 0, 'h'},                 //7
    {"version", no_argument, 0, 'V'},              //8
    {"min-mapq", required_argument, 0, 0},         //9 consider only reads with MAPQ>=min-mapq [30]
    {"secondary", required_argument, 0, 0},        //10 consider secondary alignments or not [yes]
    {"kmer-model", required_argument, 0, 0},       //11 custom nucleotide k-mer model file
    {"skip-unreadable", required_argument, 0, 0},  //12 skip any unreadable fast5 or terminate program [yes]
    {"print-events", required_argument, 0, 0},     //13 prints the event table (used for debugging)
    {"print-banded-aln", required_argument, 0, 0}, //14 prints the event alignment (used for debugging)
    {"print-scaling", required_argument, 0, 0},    //15 prints the estimated scalings (used for debugging)
    {"print-raw", required_argument, 0, 0},        //16 prints the raw signal (used for debugging)
    {"disable-cuda", required_argument, 0, 0},     //17 disable running on CUDA [no] (only if compiled for CUDA)
    {"cuda-block-size",required_argument, 0, 0},   //18
    {"debug-break",required_argument, 0, 0},       //19 break after processing the first batch (used for debugging)
    {"profile-cpu",required_argument, 0, 0},       //20 perform section by section (used for profiling - for CPU only)
    {"cuda-max-lf",required_argument, 0, 0},       //21 reads <= cuda-max-lf*avg_readlen on GPU, rest on CPU (only if compiled for CUDA)
    {"cuda-avg-epk",required_argument, 0, 0},      //22 average number of events per kmer - for allocating GPU arrays (only if compiled for CUDA)
    {"cuda-max-epk",required_argument, 0, 0},      //23 reads <= cuda_max_epk on GPU, rest on CPU (only if compiled for CUDA)
    {"cuda-dev-id",required_argument, 0, 0},       //24 cuda device ID to run on (only if compiled for CUDA)
    {"cuda-mem-frac",required_argument, 0, 0},     //25 fraction of the free GPU memory to use (only if compiled for CUDA)
    {"skip-ultra",required_argument, 0, 0},        //26 skip the ultra long reads for better load balancing
    {"ultra-thresh",required_argument, 0, 0},      //27 the threadshold for skipping ultra long reads
    {"write-dump",required_argument, 0, 0},        //28 write the raw data as a dump
    {"read-dump",required_argument, 0, 0},         //29 read the raw data as a dump
    {"output",required_argument, 0, 'o'},          //30 output to a file [stdout]
    {"iop",required_argument, 0, 0},               //31 number of I/O processes
    {"window",required_argument, 0, 'w'},          //32 the genomic window (region)
    {"summary",required_argument,0,0},             //33 summarize the alignment of each read/strand in FILE (eventalign only)
    {"sam",no_argument,0,'a'},                     //34 write output in SAM format (eventalign only)
    {"scale-events",no_argument,0,0},              //35 scale events to the model, rather than vice-versa (eventalign only)
    {"print-read-names",no_argument,0,0},          //36 print read names instead of indexes (eventalign only)
    {"samples",no_argument,0,0},                   //37 write the raw samples for the event to the tsv output (eventalign only)
    {"meth-out-version",required_argument,0,0},    //38 specify the version of the tsv output for methylation (call-methylation only)
    {"profile",required_argument, 0,'x'},          //39 profile used to tune parameters for GPU
    {"meth-model",required_argument,0,0},          //40 custom methylation k-mer model file
    {"signal-index", no_argument,0,0},             //41 write the raw signal start and end index values for the event to the tsv output (eventalign only)
    {"rna",no_argument,0,0},                       //42 if RNA (eventalign only)
    {"slow5",required_argument,0,0},               //43 read from a slow5 file
    {"min-recalib-events",required_argument,0,0},  //44 minimum number of events to recalibrate
    {"collapse-events",no_argument,0,0},           //45 collapse events that stays on the same reference k-mer
    {"pore",required_argument,0,0},                //46 pore
    {"paf",no_argument,0,'c'},                     //47 if print in paf format (only for eventalign)
    {"sam-out-version",required_argument,0,0},     //48 specify the version of the sam output for eventalign (eventalign only)
    {"m6anet",no_argument,0,0},                    //49 m6anet output (eventalign only)
    {0, 0, 0, 0}};


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


//function that processes a databatch - for pthreads when I/O and processing are interleaved
void* pthread_processor(void* voidargs) {
    pthread_arg2_t* args = (pthread_arg2_t*)voidargs;
    db_t* db = args->db;
    core_t* core = args->core;
    double realtime0=core->realtime0;

    //process
    process_db(core, db);

    fprintf(stderr, "[%s::%.3f*%.2f] %d Entries (%.1fM bases) processed\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                db->n_bam_rec,db->sum_bases/(1000.0*1000.0));

    //need to inform the output thread that we completed the processing
    pthread_mutex_lock(&args->mutex);
    pthread_cond_signal(&args->cond);
    args->finished=1;
    pthread_mutex_unlock(&args->mutex);

    if(core->opt.verbosity > 1){
        fprintf(stderr, "[%s::%.3f*%.2f] Signal sent by processor thread!\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0));
    }

    pthread_exit(0);
}


//function that prints the output and free - for pthreads when I/O and processing are interleaved
void* pthread_post_processor(void* voidargs){
    pthread_arg2_t* args = (pthread_arg2_t*)voidargs;
    db_t* db = args->db;
    core_t* core = args->core;
    double realtime0=core->realtime0;

    //wait until the processing thread has informed us
    pthread_mutex_lock(&args->mutex);
    while(args->finished==0){
        pthread_cond_wait(&args->cond, &args->mutex);
    }
    pthread_mutex_unlock(&args->mutex);

    if(core->opt.verbosity > 1){
        fprintf(stderr, "[%s::%.3f*%.2f] Signal got by post-processor thread!\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0));
    }

    //output and free
    output_db(core, db);
    free_db_tmp(db);
    free_db(db);
    free(args);
    pthread_exit(0);
}

void slow_fast5_warn(core_t *core){

    if(core->db_fast5_time  > core->process_db_time * 1.5){
        if(!(core->opt.flag & F5C_RD_SLOW5)){
            INFO("%s","Fast5 reading took more time than processing. Try increasing --iop. See http://bit.ly/f5cperf");
        }
    }
}

//todo : need to print error message and arg check with respect to eventalign
int meth_main(int argc, char* argv[], int8_t mode) {

    double realtime0 = realtime();

    //signal(SIGSEGV, sig_handler);

    const char* optstring = "r:b:g:t:B:K:v:o:x:w:hVca";

    int longindex = 0;
    int32_t c = -1;

    char* bamfilename = NULL;
    char* fastafile = NULL;
    char* fastqfile = NULL;
    char *tmpfile = NULL;
    char* profilename = NULL; //Create variable to store profile arg
    char *eventalignsummary = NULL;
    char *slow5file = NULL;

    int8_t is_ultra_thresh_set = 0;
    int8_t is_meth_out_version_set = 0;
    int8_t is_sam_out_version_set = 0;
    int8_t is_iop_set = 0;

    FILE *fp_help = stderr;

    opt_t opt;
    init_opt(&opt); //initialise options to defaults

    //parse the user args
    while ((c = getopt_long(argc, argv, optstring, long_options, &longindex)) >= 0) {
        if (c == 'r') {
            fastqfile = optarg;
        } else if (c == 'b') {
            bamfilename = optarg;
        } else if (c == 'g') {
            fastafile = optarg;

        } else if (c == 'x') {  //Set profile values. Any user-specified arguments will override the profile values.
            profilename = optarg;
            int profile_return = set_profile(profilename,&opt);
            if(profile_return == 1){
                ERROR("%s","Error occurred setting the profile.");
                exit(EXIT_FAILURE);
            }
        } else if (c == 'w') {
            opt.region_str = optarg;
        } else if (c == 'B') {
            opt.batch_size_bases = mm_parse_num(optarg);
            if(opt.batch_size_bases<=0){
                ERROR("%s","Maximum number of bases should be larger than 0.");
                exit(EXIT_FAILURE);
            }
        } else if (c == 'K') {
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
        else if (c=='c'){
            if(mode!=1){
                ERROR("%s","Option -c is available only in eventalign");
                exit(EXIT_FAILURE);
            }
            yes_or_no(&opt, F5C_PAF, longindex, "yes", 1);
        } else if (c == 'a'){ //sam output
            if(mode!=1){
                ERROR("%s","Option --sam is available only in eventalign");
                exit(EXIT_FAILURE);
            }
            yes_or_no(&opt, F5C_SAM, longindex, "yes", 1);
        } else if (c == 0 && longindex == 9) {
            opt.min_mapq = atoi(optarg); //todo : check whether this is between 0 and 60
        } else if (c == 0 && longindex == 10) { //consider secondary mappings or not
            yes_or_no(&opt, F5C_SECONDARY_YES, longindex, optarg, 1);
        } else if (c == 0 && longindex == 11) { //custom nucleotide model file
            opt.model_file = optarg;
        } else if (c == 0 && longindex == 12) {
            yes_or_no(&opt, F5C_SKIP_UNREADABLE, longindex, optarg, 1);
        } else if (c == 0 && longindex == 13) {
            yes_or_no(&opt, F5C_PRINT_EVENTS, longindex, optarg, 1);
        } else if (c == 0 && longindex == 14) {
            yes_or_no(&opt, F5C_PRINT_BANDED_ALN, longindex, optarg, 1);
        } else if (c == 0 && longindex == 15) {
            yes_or_no(&opt, F5C_PRINT_SCALING, longindex, optarg, 1);
        } else if (c == 0 && longindex == 16) {
            yes_or_no(&opt, F5C_PRINT_RAW, longindex, optarg, 1);
        } else if (c == 0 && longindex == 17) {
#ifdef HAVE_CUDA
            yes_or_no(&opt, F5C_DISABLE_CUDA, longindex, optarg, 1);
#else
            WARNING("%s", "disable-cuda has no effect when compiled for the CPU");
#endif
        } else if(c == 0 && longindex == 18){
            opt.cuda_block_size = atoi(optarg); //todo : warnining for cpu only mode, check limits
        }else if(c == 0 && longindex == 19){ //debug break
            //yes_or_no(&opt, F5C_DEBUG_BRK, longindex, optarg, 1);
            opt.debug_break = atoi(optarg);
        }else if(c == 0 && longindex == 20){ //sectional benchmark todo : warning for gpu mode
            yes_or_no(&opt, F5C_SEC_PROF, longindex, optarg, 1);
        }else if(c == 0 && longindex == 21){ //cuda todo : warning for cpu mode, error check
            opt.cuda_max_readlen = atof(optarg);
        }else if(c == 0 && longindex == 22){ //cuda todo : warning for cpu mode, warning for dynamic cuda malloc mode, error check
            opt.cuda_avg_events_per_kmer = atof(optarg);
        }else if(c == 0 && longindex == 23){ //cuda todo : warning for cpu mode, error check
            opt.cuda_max_avg_events_per_kmer = atof(optarg);
        }else if(c == 0 && longindex == 24){
            opt.cuda_dev_id = atoi(optarg);
        } else if(c == 0 && longindex == 25){ //todo : warning for CPU mode, warning for dynamic malloc mode
            opt.cuda_mem_frac = atof(optarg);
        } else if(c == 0 && longindex == 26){ //check for empty strings
            tmpfile = optarg;
        } else if(c == 0 && longindex == 27){
            opt.ultra_thresh = atoi(optarg);
            is_ultra_thresh_set = 1;
        } else if(c == 0 && longindex == 28){ //write the raw dump of the fast5 files
            yes_or_no(&opt, F5C_WR_RAW_DUMP, longindex, optarg, 1);
        } else if(c == 0 && longindex == 29){ //read the raw dump of the fast5 files
            yes_or_no(&opt, F5C_RD_RAW_DUMP, longindex, optarg, 1);
        } else if(c=='o'){
			if (strcmp(optarg, "-") != 0) {
				if (freopen(optarg, "wb", stdout) == NULL) {
					ERROR("failed to write the output to file %s : %s",optarg, strerror(errno));
					exit(EXIT_FAILURE);
				}
			}
        } else if (c == 0 && longindex == 31) {  //I/O procs
            opt.num_iop = atoi(optarg);
            if (opt.num_iop < 1) {
                ERROR("Number of I/O processes should be larger than 0. You entered %d", opt.num_iop);
                exit(EXIT_FAILURE);
            }
            is_iop_set = 1;
        } else if (c == 0 && longindex == 33){ //eventalign summary
            if(mode!=1){
                ERROR("%s","Option --summary is available only in eventalign");
                exit(EXIT_FAILURE);
            }
            eventalignsummary=optarg;
        } else if (c == 0 && longindex == 35){ //scale events
            if(mode!=1){
                ERROR("%s","Option --scale-events is available only in eventalign");
                exit(EXIT_FAILURE);
            }
            yes_or_no(&opt, F5C_SCALE_EVENTS, longindex, "yes", 1);
        } else if (c == 0 && longindex == 36){ //print read names
            if(mode!=1){
                ERROR("%s","Option --print-read-names is available only in eventalign");
                exit(EXIT_FAILURE);
            }
            yes_or_no(&opt, F5C_PRINT_RNAME, longindex, "yes", 1);
        } else if (c == 0 && longindex == 37){ //print samples
            if(mode!=1){
                ERROR("%s","Option --samples is available only in eventalign");
                exit(EXIT_FAILURE);
            }
            yes_or_no(&opt, F5C_PRINT_SAMPLES, longindex, "yes", 1);
            //ERROR ("%s","--samples not yet implemented. Submit a github request if you need this feature.");
            //exit(EXIT_FAILURE);
        } else if (c == 0 && longindex == 38){ //specify the version of the tsv output for methylation (call-methylation only)
            opt.meth_out_version=atoi(optarg);
            if(mode!=0){
                ERROR("%s","Option --meth-out-version is available only in call-methylation");
                exit(EXIT_FAILURE);
            }
            if(opt.meth_out_version<1 || opt.meth_out_version>2){
                ERROR("--meth-out-version accepts only 1 or 2. You entered %d",opt.meth_out_version);
                exit(EXIT_FAILURE);
            }
            is_meth_out_version_set = 1;
        } else if (c == 0 && longindex == 40){ //custom methylation k-mer model
            if(mode!=0){
                ERROR("%s","Option --meth-model is available only in call-methylation");
                exit(EXIT_FAILURE);
            }
            opt.meth_model_file = optarg;
        } else if (c == 0 && longindex == 41){ //write the raw signal start and end index values for the event to the tsv output
            if(mode!=1){
                ERROR("%s","Option --signal-index is available only in eventalign");
                exit(EXIT_FAILURE);
            }
            yes_or_no(&opt, F5C_PRINT_SIGNAL_INDEX, longindex, "yes", 1);
        } else if (c == 0 && longindex == 42){ //if RNA
            if(mode!=1){
                ERROR("%s","Option --rna is available only in eventalign");
                exit(EXIT_FAILURE);
            }
            yes_or_no(&opt, F5C_RNA, longindex, "yes", 1);
        } else if(c == 0 && longindex == 43){ //read from a slow5
            slow5file=optarg;
            assert(slow5file!=NULL);
            yes_or_no(&opt, F5C_RD_SLOW5, longindex, "yes", 1);
        } else if(c == 0 && longindex == 44){ //minimum number of events to rescale
            opt.min_num_events_to_rescale = atoi(optarg);
            if (opt.min_num_events_to_rescale < 1) {
                ERROR("Minimum number of events to rescale should be larger than 0. You entered %d", opt.min_num_events_to_rescale);
                exit(EXIT_FAILURE);
            }
        } else if (c == 0 && longindex == 45){ //collapse events
            if(mode!=1){
                ERROR("%s","Option --collapse-events is available only in eventalign");
                exit(EXIT_FAILURE);
            }
            yes_or_no(&opt, F5C_COLLAPSE_EVENTS, longindex, "yes", 1);
        } else if (c==0 && longindex == 46){ //pore
            opt.pore = optarg;
            if(!(strcmp(opt.pore,"r9")==0 || strcmp(opt.pore,"r10")==0 || strcmp(opt.pore,"rna004")==0)){
                ERROR("%s","Pore model should be r9, r10 or rna004");
                exit(EXIT_FAILURE);
            }
            if(strcmp(opt.pore,"r10")==0){
                opt.flag |= F5C_R10;
            } else if (strcmp(opt.pore,"rna004")==0){
                opt.flag |= F5C_RNA;
                opt.flag |= F5C_R10;
            }
        } else if (c == 0 && longindex == 48){ //specify the version of the sam output for eventalign (eventalign only)
            opt.sam_out_version=atoi(optarg);
            if(mode!=1){
                ERROR("%s","Option --sam-out-version is available only in eventalign");
                exit(EXIT_FAILURE);
            }
            if(opt.sam_out_version<1 || opt.sam_out_version>2){
                ERROR("--sam-out-version accepts only 1 or 2. You entered %d",opt.sam_out_version);
                exit(EXIT_FAILURE);
            }
            is_sam_out_version_set = 1;
        } else if (c == 0 && longindex == 49){ //m6anet output
            if(mode!=1){
                ERROR("%s","Option --m6anet is available only in eventalign");
                exit(EXIT_FAILURE);
            }
            yes_or_no(&opt, F5C_M6ANET, longindex, "yes", 1);
        }
    }

    if(is_ultra_thresh_set ==1 && tmpfile==NULL){
        WARNING("%s", "--ultra-thresh has no effect without --skip-ultra");
    }

    if(mode==0 && is_meth_out_version_set==0){
        INFO("%s","Default methylation tsv output format is changed from f5c v0.7 onwards to match latest nanopolish output. Set --meth-out-version=1 to fall back to the old format.");
    }

    if(mode==1 && is_sam_out_version_set==0 && (opt.flag&F5C_SAM)){
        INFO("%s","Default sam output format is changed from f5c v1.3 onwards. Set --sam-out-version=1 to fall back to the old format.");
    }

    if(slow5file!=NULL){
        if (is_iop_set){
            WARNING("%s","--iop option is only meant to be used for fast5 files, not slow5 files. Ignoring --iop option.");
        }
        opt.num_iop = 1;
    }

    if (fastqfile == NULL || bamfilename == NULL || fastafile == NULL || fp_help == stdout) {
        fprintf(fp_help,"Usage: f5c %s [OPTIONS] -r reads.fa -b alignments.bam -g genome.fa\n",mode==1 ? "eventalign" : "call-methylation");
        fprintf(fp_help,"       f5c %s [OPTIONS] -r reads.fa -b alignments.bam -g genome.fa --slow5 signals.blow5\n",mode==1 ? "eventalign" : "call-methylation");
        fprintf(fp_help,"\nbasic options:\n");
        fprintf(fp_help,"   -r FILE                    fastq/fasta read file\n");
        fprintf(fp_help,"   -b FILE                    sorted bam file\n");
        fprintf(fp_help,"   -g FILE                    reference genome\n");
        fprintf(fp_help,"   -w STR                     limit processing to a genomic region specified as chr:start-end or a list of regions in a .bed file\n");
        fprintf(fp_help,"   -t INT                     number of processing threads [%d]\n",opt.num_thread);
        fprintf(fp_help,"   -K INT                     batch size (max number of reads loaded at once) [%d]\n",opt.batch_size);
        fprintf(fp_help,"   -B FLOAT[K/M/G]            max number of bases loaded at once [%.1fM]\n",opt.batch_size_bases/(float)(1000*1000));
        fprintf(fp_help,"   -h                         help\n");
        fprintf(fp_help,"   -o FILE                    output to file [stdout]\n");
        fprintf(fp_help,"   -x STR                     parameter profile to be used for better performance (always applied before other options)\n"); //Added option in help
        fprintf(fp_help,"                              e.g., laptop, desktop, hpc; see https://f5c.bioinf.science/profiles for the full list\n");
        fprintf(fp_help,"   --iop INT                  number of I/O processes to read fast5 files [%d]\n",opt.num_iop);
        fprintf(fp_help,"   --pore STR                 set the pore chemistry (r9, r10 or rna004) [r9]\n");
        fprintf(fp_help,"   --slow5 FILE               read from a slow5 file\n");
        fprintf(fp_help,"   --min-mapq INT             minimum mapping quality [%d]\n",opt.min_mapq);
        fprintf(fp_help,"   --secondary=yes|no         consider secondary mappings or not [%s]\n",(opt.flag&F5C_SECONDARY_YES)?"yes":"no");
        fprintf(fp_help,"   --verbose INT              verbosity level [%d]\n",opt.verbosity);
        fprintf(fp_help,"   --version                  print version\n");

#ifdef HAVE_CUDA
        fprintf(fp_help,"   --disable-cuda=yes|no      disable running on CUDA [%s]\n",(opt.flag&F5C_DISABLE_CUDA?"yes":"no"));
        fprintf(fp_help,"   --cuda-dev-id INT          CUDA device ID to run kernels on [%d]\n",opt.cuda_dev_id);
        fprintf(fp_help,"   --cuda-max-lf FLOAT        reads with length <= cuda-max-lf*avg_readlen on GPU, rest on CPU [%.1f]\n",opt.cuda_max_readlen);
        fprintf(fp_help,"   --cuda-avg-epk FLOAT       average number of events per kmer - for allocating GPU arrays [%.1f]\n",opt.cuda_avg_events_per_kmer);
        fprintf(fp_help,"   --cuda-max-epk FLOAT       reads with events per kmer <= cuda_max_epk on GPU, rest on CPU [%.1f]\n",opt.cuda_max_avg_events_per_kmer);
#endif
        fprintf(fp_help,"\nadvanced options:\n");
        fprintf(fp_help,"   --skip-ultra FILE          skip ultra long reads and write those entries to the bam file provided as the argument\n");
        fprintf(fp_help,"   --ultra-thresh INT         threshold to skip ultra long reads [%ld]\n",(long)opt.ultra_thresh);
        fprintf(fp_help,"   --skip-unreadable=yes|no   skip any unreadable fast5 or terminate program [%s]\n",(opt.flag&F5C_SKIP_UNREADABLE?"yes":"no"));
        fprintf(fp_help,"   --kmer-model FILE          custom nucleotide k-mer model file (format similar to test/r9-models/r9.4_450bps.nucleotide.6mer.template.model)\n");
    if(mode==0){
        fprintf(fp_help,"   --meth-model FILE          custom methylation k-mer model file (format similar to test/r9-models/r9.4_450bps.cpg.6mer.template.model)\n");
        fprintf(fp_help,"   --meth-out-version INT     methylation tsv output version (set 2 to print the strand column) [%d]\n",opt.meth_out_version);
    }
    if(mode==1){
        fprintf(fp_help,"   --summary FILE             summarise the alignment of each read/strand in FILE\n");
        fprintf(fp_help,"   --paf                      write output in PAF format\n");
        fprintf(fp_help,"   --sam                      write output in SAM format\n");
        fprintf(fp_help,"   --sam-out-version INT      sam output version (set 1 to revert to old nanopolish style format) [%d]\n",opt.meth_out_version);
        fprintf(fp_help,"   --m6anet                   write output in m6anet format\n");

        fprintf(fp_help,"   --print-read-names         print read names instead of indexes\n");
        fprintf(fp_help,"   --scale-events             scale events to the model, rather than vice-versa\n");
        fprintf(fp_help,"   --samples                  write the raw samples for the event to the tsv output\n");
        fprintf(fp_help,"   --signal-index             write the raw signal start and end index values for the event to the tsv output\n");
        fprintf(fp_help,"   --rna                      the dataset is direct RNA\n");
        fprintf(fp_help,"   --collapse-events          collapse events that stays on the same reference k-mer\n");
    }
        fprintf(fp_help,"   --min-recalib-events INT   minimum number of events to recalbrate (decrease if your reads are very short and could not calibrate) [%d]\n",opt.min_num_events_to_rescale);

#ifdef HAVE_CUDA
        fprintf(fp_help,"   --cuda-mem-frac FLOAT      Fraction of free GPU memory to allocate [0.9 (0.7 for tegra)]\n");
        //fprintf(fp_help,"   --cuda-block-size\n");
#endif
        fprintf(fp_help,"\ndeveloper options:\n");
        fprintf(fp_help,"   --print-events=yes|no      prints the event table\n");
        fprintf(fp_help,"   --print-banded-aln=yes|no  prints the event alignment\n");
        fprintf(fp_help,"   --print-scaling=yes|no     prints the estimated scalings\n");
        fprintf(fp_help,"   --print-raw=yes|no         prints the raw signal\n");
        fprintf(fp_help,"   --debug-break INT          break after processing the specified no. of batches\n");
        fprintf(fp_help,"   --profile-cpu=yes|no       process section by section (used for profiling on CPU)\n");
        fprintf(fp_help,"   --write-dump=yes|no        write the fast5 dump to a file or not\n");
        fprintf(fp_help,"   --read-dump=yes|no         read from a fast5 dump file or not\n");
        fprintf(fp_help,"\nSee the manual page for details (`man ./docs/f5c.1' or https://f5c.bioinf.science/man).\n");
        if(fp_help == stdout){
            exit(EXIT_SUCCESS);
        }
        exit(EXIT_FAILURE);
    }

    //initialise the core data structure
    core_t* core = init_core(bamfilename, fastafile, fastqfile, tmpfile, opt,realtime0,mode,eventalignsummary, slow5file);

    #ifdef ESL_LOG_SUM
        p7_FLogsumInit();
    #endif

    //print the header
    if(mode==0){
        if(core->opt.meth_out_version==1){
            fprintf(stdout, "chromosome\tstart\tend\tread_name\t"
                                 "log_lik_ratio\tlog_lik_methylated\tlog_lik_unmethylated\t"
                                 "num_calling_strands\tnum_cpgs\tsequence\n");
        }
        else if (core->opt.meth_out_version==2){
            fprintf(stdout, "chromosome\tstrand\tstart\tend\tread_name\t"
                                 "log_lik_ratio\tlog_lik_methylated\tlog_lik_unmethylated\t"
                                 "num_calling_strands\tnum_motifs\tsequence\n");
        }
    }
    else if(mode==1){
        if(core->event_summary_fp!=NULL){
            fprintf(core->event_summary_fp,"read_index\tread_name\tfast5_path\tmodel_name\tstrand\tnum_events\t");
            fprintf(core->event_summary_fp,"num_steps\tnum_skips\tnum_stays\ttotal_duration\tshift\tscale\tdrift\tvar\n");
        }
        int8_t print_read_names = (core->opt.flag & F5C_PRINT_RNAME) ? 1 : 0 ;
        int8_t write_samples = (core->opt.flag & F5C_PRINT_SAMPLES) ? 1 : 0 ;
        int8_t write_signal_index = (core->opt.flag & F5C_PRINT_SIGNAL_INDEX) ? 1 : 0;
        int8_t sam_output =  (core->opt.flag & F5C_SAM) ? 1 : 0 ;
        int8_t paf_output =  (core->opt.flag & F5C_PAF) ? 1 : 0 ;
        int8_t m6anet_output =  (core->opt.flag & F5C_M6ANET) ? 1 : 0 ;

        if(sam_output && paf_output){
            ERROR("%s","-c and --sam cannot be used together");
            exit(EXIT_FAILURE);
        }
        if(sam_output && m6anet_output){
            ERROR("%s","-c and --m6anet cannot be used together");
            exit(EXIT_FAILURE);
        }
        if(paf_output && m6anet_output){
            ERROR("%s","--paf and --m6anet cannot be used together");
            exit(EXIT_FAILURE);
        }

        if(sam_output){
            emit_sam_header(core->sam_output, core->m_hdr);
        } else if (paf_output){
            //none
        } else if (m6anet_output){
            emit_event_alignment_tsv_m6anet_header(stdout, print_read_names, write_signal_index);
        } else{
            emit_event_alignment_tsv_header(stdout, print_read_names, write_samples, write_signal_index);
        }

    }
    int32_t counter=0;

 #ifdef IO_PROC_NO_INTERLEAVE   //If input, processing and output are not interleaved (serial mode)

    //initialise a databatch
    db_t* db = init_db(core);

    ret_status_t status = {core->opt.batch_size,core->opt.batch_size_bases};
    while (status.num_reads >= core->opt.batch_size || status.num_bases>=core->opt.batch_size_bases) {

        //load a databatch
        status = load_db(core, db);

        fprintf(stderr, "[%s::%.3f*%.2f] %d Entries (%.1fM bases) loaded\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                status.num_reads,status.num_bases/(1000.0*1000.0));

        //process a databatch
        process_db(core, db);

        fprintf(stderr, "[%s::%.3f*%.2f] %d Entries (%.1fM bases) processed\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                status.num_reads,status.num_bases/(1000.0*1000.0));

        //output print
        output_db(core, db);

        //free temporary
        free_db_tmp(db);

        slow_fast5_warn(core);

        if(opt.debug_break==counter){
            break;
        }
        counter++;
    }

    //free the databatch
    free_db(db);

#else   //input, processing and output are interleaved (default)

    ret_status_t status = {core->opt.batch_size,core->opt.batch_size_bases};
    int8_t first_flag_p=0;
    int8_t first_flag_pp=0;
    pthread_t tid_p; //process thread
    pthread_t tid_pp; //post-process thread


    while (status.num_reads >= core->opt.batch_size || status.num_bases>=core->opt.batch_size_bases) {

        //init and load a databatch
        db_t* db = init_db(core);
        status = load_db(core, db);

        fprintf(stderr, "[%s::%.3f*%.2f] %d Entries (%.1fM bases) loaded\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                status.num_reads,status.num_bases/(1000.0*1000.0));

        if(first_flag_p){ //if not the first time of the "process" wait for the previous "process"
            int ret = pthread_join(tid_p, NULL);
            NEG_CHK(ret);
            if(opt.verbosity>1){
                fprintf(stderr, "[%s::%.3f*%.2f] Joined to processor thread %ld\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                (long)tid_p);
            }
            slow_fast5_warn(core);
        }
        first_flag_p=1;

        //set up args
        pthread_arg2_t *pt_arg = (pthread_arg2_t*)malloc(sizeof(pthread_arg2_t));
        pt_arg->core=core;
        pt_arg->db=db;
        pt_arg->cond = PTHREAD_COND_INITIALIZER;
        pt_arg->mutex = PTHREAD_MUTEX_INITIALIZER;
        pt_arg->finished = 0;

        //process thread launch
        int ret = pthread_create(&tid_p, NULL, pthread_processor,
                                (void*)(pt_arg));
        NEG_CHK(ret);
        if(opt.verbosity>1){
            fprintf(stderr, "[%s::%.3f*%.2f] Spawned processor thread %ld\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                (long)tid_p);
        }

        if(first_flag_pp){ //if not the first time of the post-process wait for the previous post-process
            int ret = pthread_join(tid_pp, NULL);
            NEG_CHK(ret);
            if(opt.verbosity>1){
                fprintf(stderr, "[%s::%.3f*%.2f] Joined to post-processor thread %ld\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                (long)tid_pp);
            }
        }
        first_flag_pp=1;

        //post-process thread launch (output and freeing thread)
        ret = pthread_create(&tid_pp, NULL, pthread_post_processor,
                                (void*)(pt_arg));
        NEG_CHK(ret);
        if(opt.verbosity>1){
            fprintf(stderr, "[%s::%.3f*%.2f] Spawned post-processor thread %ld\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                (long)tid_pp);
        }

        if(opt.debug_break==counter){
            break;
        }
        counter++;
    }

    //final round
    int ret = pthread_join(tid_p, NULL);
    NEG_CHK(ret);
    if(opt.verbosity>1){
        fprintf(stderr, "[%s::%.3f*%.2f] Joined to last processor thread %ld\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                (long)tid_p);
    }
    ret = pthread_join(tid_pp, NULL);
    NEG_CHK(ret);
    if(opt.verbosity>1){
    fprintf(stderr, "[%s::%.3f*%.2f] Joined to last post-processor thread %ld\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                (long)tid_pp);
    }


#endif

    fprintf(stderr, "[%s] skipped unmapped: %ld, skipped secondary: %ld, skipped low_mapq: %ld\n", __func__,(long)core->unmapped_reads, (long)core->skip_sec_reads, (long)core->skip_mapq_reads);
    fprintf(stderr, "[%s] total entries: %ld, qc fail: %ld, could not calibrate: %ld, no alignment: %ld, bad reads: %ld",
             __func__,(long)core->total_reads, (long)core->qc_fail_reads, (long)core->failed_calibration_reads, (long)core->failed_alignment_reads, (long)core->bad_fast5_file);
    fprintf(stderr,"\n[%s] total bases: %.1f Mbases",__func__,core->sum_bases/(float)(1000*1000));

    fprintf(stderr, "\n[%s] Data loading time: %.3f sec", __func__,core->load_db_time);
    fprintf(stderr, "\n[%s]     - bam load time: %.3f sec", __func__, core->db_bam_time);
    fprintf(stderr, "\n[%s]     - fasta load time: %.3f sec", __func__, core->db_fasta_time);
    if(core->opt.flag&F5C_RD_SLOW5){
        fprintf(stderr, "\n[%s]     - slow5 load time: %.3f sec", __func__, core->db_fast5_time);
    } else{
        fprintf(stderr, "\n[%s]     - fast5 load time: %.3f sec", __func__, core->db_fast5_time);
        fprintf(stderr, "\n[%s]         - fast5 open time: %.3f sec", __func__, core->db_fast5_open_time);
        fprintf(stderr, "\n[%s]         - fast5 read time: %.3f sec", __func__, core->db_fast5_read_time);
    }

    fprintf(stderr, "\n[%s] Data processing time: %.3f sec", __func__,core->process_db_time);

    if((core->opt.flag&F5C_SEC_PROF) || (!(core->opt.flag & F5C_DISABLE_CUDA))){
        fprintf(stderr, "\n[%s]     - Events time: %.3f sec",
                __func__, core->event_time);
        fprintf(stderr, "\n[%s]     - Alignment time: %.3f sec",
                __func__, core->align_time);
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
        fprintf(stderr, "\n[%s]     - Estimate scaling time: %.3f sec",
                __func__, core->est_scale_time);
        fprintf(stderr, "\n[%s]     - HMM time: %.3f sec",
                __func__, core->meth_time);

    }
    fprintf(stderr, "\n[%s] Data output time: %.3f sec", __func__,core->output_time);

    fprintf(stderr,"\n");

    if(core->ultra_long_skipped>0){
        assert(tmpfile!=NULL);
        WARNING("%ld ultra long reads (>%.1f kbases) were skipped.",(long)core->ultra_long_skipped,core->opt.ultra_thresh/1000.0);
        fprintf(stderr," Please run samtools index on '%s' followed by f5c with a larger -B on the CPU.\n",tmpfile);
    }


#ifndef IO_PROC_NO_INTERLEAVE
    if((core->load_db_time - core->process_db_time) > (core->process_db_time*0.2) ){
        INFO("Performance bounded by file I/O. File I/O took %.3f sec than processing",core->load_db_time - core->process_db_time);
    }
#endif

    #ifdef HAVE_CUDA
        if(opt.verbosity>0){
            fprintf(stderr, "[%s] max-lf: %.2f, avg-epk: %.2f, max-epk: %.2f, K: %d, B: %ld, T: %d, Ultra: %ld, Align: %.3f, Diff: %.3f\n", __func__,
            opt.cuda_max_readlen,opt.cuda_avg_events_per_kmer,opt.cuda_max_avg_events_per_kmer,opt.batch_size,opt.batch_size_bases,opt.num_thread,opt.ultra_thresh,
            core->align_time,core->extra_load_cpu);
        }
    #endif

    if(core->skip_mapq_reads > core->total_reads){
        WARNING("Skipped %ld reads with MAPQ < %d. Use --min-mapq to change the threshold.",(long)core->skip_mapq_reads,opt.min_mapq);
    }
    if(core->qc_fail_reads + core->failed_calibration_reads + core->failed_alignment_reads > core->total_reads/2){
        WARNING("%ld out of %ld reads failed. Double check --pore and --rna options.",(long)core->qc_fail_reads + core->failed_calibration_reads + core->failed_alignment_reads,(long)core->total_reads);
    }
    if(core->bad_fast5_file > core->total_reads/10.0){
        WARNING("%ld out of %ld reads missing in FAST5/SLOW5.",(long)core->bad_fast5_file,(long)core->total_reads);
    }
    if(core->qc_fail_reads + core->failed_calibration_reads + core->failed_alignment_reads + core->bad_fast5_file >= core->total_reads){
        ERROR("all %ld out of %ld reads failed. Check the input files.",(long)core->total_reads,(long)core->total_reads);
        exit(EXIT_FAILURE);
    }
    if(core->total_reads==0){
        ERROR("%s","0 reads processed. Check the input files.");
        exit(EXIT_FAILURE);
    }
    //free the core data structure
    free_core(core,opt);

    return 0;
}
