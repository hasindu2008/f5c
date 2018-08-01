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

double realtime0 = realtime();

static struct option long_options[] = {
    {"reads", required_argument, 0, 'r'},     //0
    {"bam", required_argument, 0, 'b'},       //1
    {"genome", required_argument, 0, 'g'},    //2
    {"threads", required_argument, 0, 'g'},   //3
    {"batchsize", required_argument, 0, 'K'}, //4
    {"print", no_argument, 0, 'p'},           //5
    {"aaaaaa", no_argument, 0, 0},            //6
    {"help", no_argument, 0, 'h'},            //7
    {"version", no_argument, 0, 'V'},         //8
    {"min-mapq", required_argument, 0, 0},    //9
    {"secondary", required_argument, 0, 0},   //10
    {"kmer-model", required_argument, 0, 0},  //11
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

int main(int argc, char* argv[]) {
    signal(SIGSEGV, sig_handler);

    const char* optstring = "r:b:g:hvp";
    int longindex = 0;
    int c = -1;

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
            opt.print_raw = 1;
        } else if (c == 0 && longindex == 11) {
            opt.model_file = optarg;
        }
    }

    if (fastqfile == NULL || bamfilename == NULL || fastafile == NULL) {
        fprintf(
            stderr,
            "Usage: %s [OPTIONS] -r reads.fa -b alignments.bam -g genome.fa\n",
            argv[0]);
        exit(EXIT_FAILURE);
    }

    core_t* core = init_core(bamfilename, fastafile, fastqfile, opt);
    db_t* db = init_db();

    int32_t status = db->capacity_bam_rec;
    while (status >= db->capacity_bam_rec) {
        status = load_db(core, db);

        fprintf(stderr, "[%s::%.3f*%.2f] %d Entries loaded\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                status);

#ifdef HAVE_CUDA
        //this is very unoptimal just for testing
#else
        process_db(core, db);
#endif
        fprintf(stderr, "[%s::%.3f*%.2f] %d Entries processed\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),
                status);

        free_db_tmp(db);
    }
    free_db(db);
    free_core(core);
    fprintf(stderr, "[%s] CMD:", __func__);
    for (int i = 0; i < argc; ++i) {
        fprintf(stderr, " %s", argv[i]);
    }

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
