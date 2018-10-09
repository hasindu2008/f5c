#include "f5c.h"
#include "f5cmisc.h"
#include <assert.h>
#include <execinfo.h>
#include <getopt.h>
#include <signal.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

//double realtime0 = 0;

static struct option long_options[] = {
    {"reads", required_argument, 0, 'r'},          //0
    {"bam", required_argument, 0, 'b'},            //1
    {"genome", required_argument, 0, 'g'},         //2
    {"threads", required_argument, 0, 't'},        //3
    {"batchsize", required_argument, 0, 'K'},      //4
    {"print", no_argument, 0, 'p'},                //5
    {"aaaaaa", no_argument, 0, 0},                 //6
    {"help", no_argument, 0, 'h'},                 //7
    {"version", no_argument, 0, 'V'},              //8
    {"min-mapq", required_argument, 0, 0},         //9
    {"secondary", required_argument, 0, 0},        //10
    {"kmer-model", required_argument, 0, 0},       //11
    {"skip-unreadable", required_argument, 0, 0},  //12
    {"print-events", required_argument, 0, 0},     //13
    {"print-banded-aln", required_argument, 0, 0}, //14
    {"print-scaling", required_argument, 0, 0},    //15
    {"print-raw", required_argument, 0, 0},        //16
    {"disable-cuda", required_argument, 0, 0},     //17
    {"cuda-block-size",required_argument, 0, 0},   //18
    {"debug-break",required_argument, 0, 0},   //19
    {0, 0, 0, 0}};

void sig_handler(int sig) {
    void* array[100];
    size_t size = backtrace(array, 100);
    ERROR("I regret to inform that a segmentation fault occurred. But at least "
          "it is better than a wrong answer%s",
          ".");
    fprintf(stderr,
            "[%s::DEBUG]\033[1;35m Here is the backtrace in case it is of any "
            "use:\n",
            __func__);
    backtrace_symbols_fd(&array[2], size - 1, STDERR_FILENO);
    fprintf(stderr, "\033[0m\n");
    exit(EXIT_FAILURE);
}

static inline int64_t mm_parse_num(const char* str) //taken from minimap2
{
    double x;
    char* p;
    x = strtod(str, &p);
    if (*p == 'G' || *p == 'g')
        x *= 1e9;
    else if (*p == 'M' || *p == 'm')
        x *= 1e6;
    else if (*p == 'K' || *p == 'k')
        x *= 1e3;
    return (int64_t)(x + .499);
}

static inline void yes_or_no(opt_t* opt, uint64_t flag, int long_idx,
                             const char* arg,
                             int yes_to_set) //taken from minimap2
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



void* pthread_processor(void* voidargs) {
    pthread_arg2_t* args = (pthread_arg2_t*)voidargs;
    db_t* db = args->db;
    core_t* core = args->core;
    double realtime0=core->realtime0;

    process_db(core, db);

    fprintf(stderr, "[%s::%.3f*%.2f] %d Entries processed\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                db->n_bam_rec);

    //need to say that we processed stuff
    pthread_mutex_lock(&args->mutex);
    pthread_cond_signal(&args->cond);
    pthread_mutex_unlock(&args->mutex); 

    fprintf(stderr, "[%s::%.3f*%.2f] Signal sent!\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0));

    pthread_exit(0);
}


void* pthread_post_processor(void* voidargs){
    pthread_arg2_t* args = (pthread_arg2_t*)voidargs;
    db_t* db = args->db;
    core_t* core = args->core;
    double realtime0=core->realtime0;

    pthread_mutex_lock(&args->mutex);
    pthread_cond_wait(&args->cond, &args->mutex);
    pthread_mutex_unlock(&args->mutex);

    fprintf(stderr, "[%s::%.3f*%.2f] Signal got!\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0));

    output_db(core, db);
    free_db_tmp(db);
    free_db(db);    
    free(args);
    pthread_exit(0);    
}

int main(int argc, char* argv[]) {
    
    double realtime0 = realtime();

    signal(SIGSEGV, sig_handler);

    const char* optstring = "r:b:g:t:K:hvp";
    int longindex = 0;
    int32_t c = -1;

    char* bamfilename = NULL;
    char* fastafile = NULL;
    char* fastqfile = NULL;

    opt_t opt;
    init_opt(&opt);

    while ((c = getopt_long(argc, argv, optstring, long_options, &longindex)) >=
           0) {
        if (c == 'r') {
            fastqfile = optarg;
        } else if (c == 'b') {
            bamfilename = optarg;
        } else if (c == 'g') {
            fastafile = optarg;
        } else if (c == 'p') {
            opt.flag |= F5C_PRINT_RAW;
        } else if (c == 'K') {
            opt.batch_size = atoi(optarg);
            if (opt.batch_size < 1) {
                ERROR("Batch size should larger than 0. You entered %d",
                      opt.batch_size);
                exit(EXIT_FAILURE);
            }
        } else if (c == 't') {
            opt.num_thread = atoi(optarg);
            if (opt.num_thread < 1) {
                ERROR("Number of threads should larger than 0. You entered %d",
                      opt.num_thread);
                exit(EXIT_FAILURE);
            }
        } else if (c == 0 && longindex == 9) {
            opt.min_mapq =
                atoi(optarg); //check whether this is between 0 and 60
        } else if (c == 0 && longindex == 10) { //consider secondary
            yes_or_no(&opt, F5C_SECONDARY_YES, longindex, optarg, 1);
        } else if (c == 0 && longindex == 11) {
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
            WARNING("%s",
                    "disable-cuda has no effect when compiled for the CPU");
#endif
        } else if(c == 0 && longindex == 18){
            opt.cuda_block_size = atoi(optarg); //todo : warnining for cpu only mode, check limits
        }else if(c == 0 && longindex == 19){
            yes_or_no(&opt, F5C_DEBUG_BRK, longindex, optarg, 1);
        }

    }

    if (fastqfile == NULL || bamfilename == NULL || fastafile == NULL) {
        fprintf(
            stderr,
            "Usage: %s [OPTIONS] -r reads.fa -b alignments.bam -g genome.fa\n",
            argv[0]);
        exit(EXIT_FAILURE);
    }

    core_t* core = init_core(bamfilename, fastafile, fastqfile, opt,realtime0);

 #ifndef IO_PROC_INTERLEAVE   
    db_t* db = init_db(core);

    int32_t status = db->capacity_bam_rec;
    while (status >= db->capacity_bam_rec) {
        status = load_db(core, db);

        fprintf(stderr, "[%s::%.3f*%.2f] %d Entries loaded\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                status);

        process_db(core, db);

        fprintf(stderr, "[%s::%.3f*%.2f] %d Entries processed\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                status);
        output_db(core, db);
        free_db_tmp(db);
        if(opt.flag & F5C_DEBUG_BRK){
            break;
        }
    }
    free_db(db);

#else
    int32_t status = core->opt.batch_size;
    int8_t first_flag_p=0;
    int8_t first_flag_pp=0;
    pthread_t tid_p;
    pthread_t tid_pp;


    while (status >= core->opt.batch_size) {
        db_t* db = init_db(core);
        status = load_db(core, db);

        fprintf(stderr, "[%s::%.3f*%.2f] %d Entries loaded\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                status);

        if(first_flag_p){
            int ret = pthread_join(tid_p, NULL);
            fprintf(stderr, "[%s::%.3f*%.2f] Joined to thread %lu\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                tid_p);
            NEG_CHK(ret);
        }
        first_flag_p=1;

        pthread_arg2_t *pt_arg = (pthread_arg2_t*)malloc(sizeof(pthread_arg2_t));
        pt_arg->core=core;
        pt_arg->db=db;
        pt_arg->cond = PTHREAD_COND_INITIALIZER;
        pt_arg->mutex = PTHREAD_MUTEX_INITIALIZER;
        int ret = pthread_create(&tid_p, NULL, pthread_processor,
                                (void*)(pt_arg));
        NEG_CHK(ret);
        fprintf(stderr, "[%s::%.3f*%.2f] Spawned thread %lu\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                tid_p);

        if(first_flag_pp){
            int ret = pthread_join(tid_pp, NULL);
            fprintf(stderr, "[%s::%.3f*%.2f] Joined to thread %lu\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                tid_pp);
            NEG_CHK(ret);
        }
        first_flag_pp=1;


        ret = pthread_create(&tid_pp, NULL, pthread_post_processor,
                                (void*)(pt_arg));
        NEG_CHK(ret);
        fprintf(stderr, "[%s::%.3f*%.2f] Spawned thread %lu\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                tid_pp);

        if(opt.flag & F5C_DEBUG_BRK){
            break;
        }
    }
    int ret = pthread_join(tid_p, NULL);
    NEG_CHK(ret);
    fprintf(stderr, "[%s::%.3f*%.2f] Joined to thread %lu\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                tid_p);
    ret = pthread_join(tid_pp, NULL);
    NEG_CHK(ret);
    fprintf(stderr, "[%s::%.3f*%.2f] Joined to thread %lu\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                tid_pp);


#endif


    fprintf(stderr, "[%s] CMD:", __func__);
    for (int i = 0; i < argc; ++i) {
        fprintf(stderr, " %s", argv[i]);
    }

#ifdef SECTIONAL_BENCHMARK 
    fprintf(stderr, "\n\n[%s] Alignment time: %.3f sec",
            __func__, core->align_time);
#endif

    free_core(core);
    fprintf(stderr, "\n[%s] Real time: %.3f sec; CPU time: %.3f sec\n",
            __func__, realtime() - realtime0, cputime());
    // }

    // else {
    // fprintf(stderr,"Usage: %s [OPTIONS] -r reads.fa -b alignments.bam -g
    // genome.fa\n",argv[0]);
    // exit(EXIT_FAILURE);
    // }


    return 0;
}
