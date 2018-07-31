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

static struct option long_options[] = {{"reads", required_argument, 0, 'r'},
                                       {"bam", required_argument, 0, 'b'},
                                       {"genome", required_argument, 0, 'g'},
                                       {"threads", required_argument, 0, 'g'},
                                       {"batchsize", required_argument, 0, 'K'},
                                       {"print", no_argument, 0, 'p'},
                                       {"aaaaaa", no_argument, 0, 0},
                                       {"help", no_argument, 0, 'h'},
                                       {"version", no_argument, 0, 'V'},
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
        switch (c) {
            case 'r':
                fastqfile = optarg;
                break;
                // if(fastqfile==NULL)
                //    fprintf(stderr,ERR "Fastq cannot be null" CEND);
            case 'b':
                bamfilename = optarg;
                break;
            case 'g':
                fastafile = optarg;
                break;
            case 'p':
                opt.print_raw = 1;
                break;
                /*default:
      //puts("dsfgsdgsd");
      //fprintf(stderr,ERR "usage : " CEND);
      //exit(EXIT_FAILURE);
  */
        }
    }

    // if (strcpy(argv[optind],"view")==0){
    // }
    // else if ((strcpy(argv[optind],"call-methylation")==0)){

    if (fastqfile == NULL || bamfilename == NULL || fastafile == NULL) {
        fprintf(
            stderr,
            "Usage: %s [OPTIONS] -r reads.fa -b alignments.bam -g genome.fa\n",
            argv[0]);
        exit(EXIT_FAILURE);
    }

    // const char *bamfilename="/mnt/f/share/778Nanopore/fastq/740475-67.bam";
    // const char *fastafile="/mnt/f/share/reference/hg38noAlt.fa";
    // const char *fastqfile="/mnt/f/share/778Nanopore/fastq/740475-67.fastq";

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
