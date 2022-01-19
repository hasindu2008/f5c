/* @file freq.c
**
** calculate methylation frequency using per-read methylation calls in tsv format
** @author: Thomas Lam
** @author: Hasindu Gamaarachchi (hasindu@unsw.edu.au)
** @@
******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include "error.h"
#include "khash.h"
#include "f5c.h"

#define KEY_SIZE 3

#define ALLOC(type, size) (type *)safe_malloc((size) * sizeof(type))

// TODO : a number of inefficient mallocs are done in this code, which can be removed
// for example getline can use a onetime buffer
// need to check for empty newline chars
// todo : can improve peak memory by disk based sorting - futuristic next level

static const char usage[] = "Usage: %s [OPTIONS] -i methcalls.tsv\n"
                            "\n"
                            "  -c FLOAT        call threshold for the log likelihood ratio. Default is 2.5. \n"
                            "  -i FILE         input file. Read from stdin if not specified\n"
                            "  -o FILE         output file. Write to stdout if not specified\n"
                            "  -s              split groups\n"
                            "  -h              help\n"
                            "  --version       print version\n"
                            "\nSee the manual page for details (`man ./docs/f5c.1' or https://f5c.page.link/man).\n"
                            ;

struct site_stats {
    int num_reads;
    int posterior_methylated;
    int called_sites;
    int called_sites_methylated;
    int group_size;
    char* sequence;
};

struct tsv_record {
    char* chromosome;
    int start;
    int end;
    double log_lik_ratio;
    int num_cpgs;
    char *sequence;
};

KHASH_MAP_INIT_STR(str, struct site_stats);

void* safe_malloc(size_t size) {
    void* p = malloc(size);
    if (!p) {
        perror("Cannot allocated memory");
        exit(EXIT_FAILURE);
    }
    return p;
}

/* Format chromosome, start and end to tab-delimited string. */
char* make_key(char* chromosome, int start, int end) {
    int start_strlen = snprintf(NULL, 0, "%d", start);
    int end_strlen = snprintf(NULL, 0, "%d", end);
    int key_strlen = strlen(chromosome) + start_strlen + end_strlen + 3;
    char* key = ALLOC(char, key_strlen);
    snprintf(key, key_strlen, "%s\t%d\t%d", chromosome, start, end);
    return key;
}

/* Split tab-delimited keys */
char** split_key(char* key, int size) {
    char** tok = ALLOC(char*, size);
    char* cpy = strdup(key);
    char* cpy_start = cpy;

    char* t = strtok(cpy, "\t");
    if (t) {
        tok[0] = strdup(t);
    }

    for (int i = 1; i < size; i++) {
        char* t = strtok(NULL, "\t");
        if (t) {
            tok[i] = strdup(t);
        }
    }

    free(cpy_start);

    return tok;
}

int cmp_key(const void* a, const void* b) {
    char* key_a = *(char**)a;
    char* key_b = *(char**)b;

    char** toks_a = split_key(key_a, KEY_SIZE);
    char** toks_b = split_key(key_b, KEY_SIZE);

    int chromosome_neq = strcmp(toks_a[0], toks_b[0]);

    if (chromosome_neq) {
        for (int i = 0; i < KEY_SIZE; i++) {
            free(toks_a[i]);
            free(toks_b[i]);
        }
        free(toks_a);
        free(toks_b);
        return chromosome_neq;
    }

    int start_a = atoi(toks_a[1]);
    int start_b = atoi(toks_b[1]);
    int end_a = atoi(toks_a[2]);
    int end_b = atoi(toks_b[2]);

    for (int i = 0; i < KEY_SIZE; i++) {
        free(toks_a[i]);
        free(toks_b[i]);
    }
    free(toks_a);
    free(toks_b);

    if (start_a == start_b) {
        return end_a - end_b;
    }

    return start_a - start_b;
}

void update_call_stats(khash_t(str)* sites, char* key, int num_called_cpg_sites,
        bool is_methylated, char* sequence) {
    struct site_stats ss = {
        .num_reads = 0,
        .posterior_methylated = 0,
        .called_sites = 0,
        .called_sites_methylated = 0,
        .group_size = num_called_cpg_sites,
        .sequence = strdup(sequence)
    };

    int absent;
    khint_t k = kh_put(str, sites, key, &absent);

    if (absent == -1) {
        fprintf(stderr, "Failed to insert key: %s\n", key);
        exit(EXIT_FAILURE);
    } else if (absent > 0) {
        kh_key(sites,k)=key;
        kh_value(sites, k) = ss;
    }
    else{
        free(ss.sequence);
        free(key);
    }

    kh_value(sites, k).num_reads++;
    kh_value(sites, k).called_sites += num_called_cpg_sites;

    if (is_methylated) {
        kh_value(sites, k).called_sites_methylated += num_called_cpg_sites;
    }
}

static inline void strtok_null_check(char *tmp, int64_t line_num){
    if(tmp == NULL){
        fprintf(stderr,"Corrupted tsv file? Offending line is line number %ld (1-based index)\n",(long)line_num);
        exit(EXIT_FAILURE);
    }
}

struct tsv_record* get_tsv_line(FILE* fp, int8_t meth_out_version, int64_t line_num) {
    struct tsv_record* record = ALLOC(struct tsv_record, 1);
    char* buf = NULL;
    size_t buf_size = 0;
    char* tmp;

    if (getline(&buf, &buf_size, fp) == -1) {
        free(record);
        if(buf_size>0){
            free(buf);
        }
        return NULL;

    }

    //chromosome
    tmp = strtok(buf, "\t");
    strtok_null_check(tmp,line_num);
    record->chromosome = strdup(tmp);

    //strand
    if(meth_out_version==2){
        tmp = strtok(NULL, "\t"); //ignored
        strtok_null_check(tmp,line_num);;
    }

    //start
    tmp = strtok(NULL, "\t");
    strtok_null_check(tmp,line_num);;
    record->start = atoi(tmp);

    //end
    tmp = strtok(NULL, "\t");
    strtok_null_check(tmp,line_num);;
    record->end = atoi(tmp);

    if(record->start<0 || record->end<0){
        fprintf(stderr,"chromosome coordinates cannot be negative. Offending line is line number %ld (1-based index)\n",(long)line_num);
        exit(EXIT_FAILURE);
    }

    tmp = strtok(NULL, "\t"); //readname
    strtok_null_check(tmp,line_num);;

    //log_lik_ratio
    tmp = strtok(NULL, "\t");
    strtok_null_check(tmp,line_num);;
    record->log_lik_ratio = strtod(tmp, NULL);

    tmp = strtok(NULL, "\t"); //log_lik_methylated
    tmp = strtok(NULL, "\t"); //log_lik_unmethylated
    tmp = strtok(NULL, "\t"); //num_calling_strands
    strtok_null_check(tmp,line_num);;

    //num_cpgs
    tmp = strtok(NULL, "\t");
    strtok_null_check(tmp,line_num);;
    record->num_cpgs = atoi(tmp);
    if(record->num_cpgs<0){
        fprintf(stderr,"num_cpgs cannot be negative. Offending line is line number %ld (1-based index)\n",(long)line_num);
        exit(EXIT_FAILURE);
    }

    //sequence
    tmp = strtok(NULL, "\t\n");
    strtok_null_check(tmp,line_num);;
    record->sequence = strdup(tmp);

    if(strtok(NULL, "\t\n")!=NULL){
        fprintf(stderr,"encountered more columns than expected. Offending line is line number %ld (1-based index)\n",(long)line_num);
        exit(EXIT_FAILURE);
    }

    free(buf);

    return record;
}

int freq_main(int argc, char **argv) {
    FILE* input = stdin;
    double call_threshold = 2.5;
    bool split_groups = false;

    int c;
    int8_t meth_out_version=0; //the output from methylation call - if strand column exists meth_out_version=2  or else meth_out_version=1
    int8_t meth_out_cpg_or_motif=0; //the output from methylation call - if the columns are named num_cpg meth_out_cpg_or_motif =1, if num_motifs then meth_out_cpg_or_motif=2
    int64_t line_num=0;

    int longindex = 0;
    const char *short_opts = "c:i:o:shV";
    //extern char* optarg;
    //extern int optind, optopt;

    static struct option long_opts[] = {
        {"call-threshold", required_argument, 0, 'c'},
        {"input", required_argument, 0, 'i'},
        {"output", required_argument, 0, 'o'},
        {"split-groups", no_argument, 0, 's'},
        {"help",no_argument, 0, 'h'},
        {"version",no_argument, 0, 'V'},
        {0, 0, 0, 0}};

    while ((c = getopt_long(argc, argv, short_opts, long_opts, &longindex)) != -1) {
        switch(c) {
            case 'c':
                /* TODO handle overflow when calling strtod */
                call_threshold = strtod(optarg, NULL);
                break;
            case 'i':
            {
                FILE* in = fopen(optarg, "r");
                F_CHK(in,optarg)
                input = in;
                break;
            }
            case 's':
                split_groups = true;
                break;
            case ':':
                fprintf(stderr, "Option -%c requires an operand\n", optopt);
                fprintf(stderr, usage, argv[0]);
                exit(EXIT_FAILURE);
            case '?':
                fprintf(stderr, "Unrecognized option: -%c\n", optopt);
                fprintf(stderr, usage, argv[0]);
                exit(EXIT_FAILURE);
            case 'o':
                if (strcmp(optarg, "-") != 0) {
                    if (freopen(optarg, "wb", stdout) == NULL) {
                        ERROR("failed to write the output to file %s : %s",optarg, strerror(errno));
                        exit(EXIT_FAILURE);
                    }
                }
                break;
            case 'h':
                fprintf(stdout, usage, argv[0]);
                exit(EXIT_SUCCESS);
            case 'V':
                fprintf(stdout,"F5C %s\n",F5C_VERSION);
                exit(EXIT_SUCCESS);
            default:
                break;
        }
    }

    if(input==stdin){
        fprintf(stderr,"Scanning the input from stdin ...\n");
    }

    khash_t(str)* sites = kh_init(str);
    struct tsv_record* record;
    /* check the header */
    char *tmp = NULL;
    size_t tmp_size = 0;
    ssize_t ret = getline(&tmp, &tmp_size, input);
    if(ret==-1){
        fprintf(stderr,"Bad file format with no header?\n");
        exit(EXIT_FAILURE);
    }
    if(strcmp(tmp,"chromosome\tstart\tend\tread_name\tlog_lik_ratio\tlog_lik_methylated\tlog_lik_unmethylated\tnum_calling_strands\tnum_cpgs\tsequence\n")==0){
        meth_out_version=1;
        meth_out_cpg_or_motif=1;
    }
    else if(strcmp(tmp,"chromosome\tstart\tend\tread_name\tlog_lik_ratio\tlog_lik_methylated\tlog_lik_unmethylated\tnum_calling_strands\tnum_motifs\tsequence\n")==0){
        meth_out_version=1;
        meth_out_cpg_or_motif=2;
    }
    else if(strcmp(tmp,"chromosome\tstrand\tstart\tend\tread_name\tlog_lik_ratio\tlog_lik_methylated\tlog_lik_unmethylated\tnum_calling_strands\tnum_cpgs\tsequence\n")==0){
        meth_out_version=2;
        meth_out_cpg_or_motif=1;
    }
    else if(strcmp(tmp,"chromosome\tstrand\tstart\tend\tread_name\tlog_lik_ratio\tlog_lik_methylated\tlog_lik_unmethylated\tnum_calling_strands\tnum_motifs\tsequence\n")==0){
        meth_out_version=2;
        meth_out_cpg_or_motif=2;
    }
    else{
        fprintf(stderr, "Incorrect header: %s\n", tmp);
        exit(EXIT_FAILURE);
    }
    line_num++;
    line_num++;

    while ((record = get_tsv_line(input,meth_out_version,line_num)) != NULL) {
        line_num++;
        int num_sites = record->num_cpgs;
        double llr = record->log_lik_ratio;

        if (fabs(llr) < call_threshold) {
            free(record->sequence);
            free(record->chromosome);
            free(record);
            continue;
        }

        char* sequence = record->sequence;
        bool is_methylated = llr > 0;

        if (split_groups && num_sites > 1) {
            char* c = record->chromosome;
            int s = record->start;
            /* TODO following variable is unused in original Python script */
            /* int e = records[i].end; */
            char* substring = strstr(sequence, "CG");
            /* str.find() is equivalent to the offset of the needle pointer location */
            /* and the haystack pointer location */
            long cg_pos = substring == NULL ? -1 : substring - sequence;
            long first_cg_pos = cg_pos;

            while (cg_pos != -1) {
                char* key = make_key(c, s + cg_pos - first_cg_pos, s + cg_pos - first_cg_pos);
                const char* sg = "split-group";
                update_call_stats(sites, key, 1, is_methylated, (char*)sg);
                char* substring = strstr(sequence + cg_pos + 1, "CG");
                cg_pos = substring == NULL ? -1 : substring - sequence;
            }
        } else {
            char* key = make_key(record->chromosome, record->start, record->end);
            update_call_stats(sites, key, num_sites, is_methylated, sequence);

        }
        free(record->sequence);
        free(record->chromosome);
        free(record);
    }

    if(meth_out_cpg_or_motif==1){
        printf("chromosome\tstart\tend\tnum_cpgs_in_group\tcalled_sites\tcalled_sites_methylated\tmethylated_frequency\tgroup_sequence\n");
    }
    else if(meth_out_cpg_or_motif==2){
        printf("chromosome\tstart\tend\tnum_motifs_in_group\tcalled_sites\tcalled_sites_methylated\tmethylated_frequency\tgroup_sequence\n");
    }

    char** sorted_keys = ALLOC(char*, kh_size(sites));
    int size = 0;

    for (khint_t k = kh_begin(sites); k != kh_end(sites); k++) {
        if (kh_exist(sites, k)) {
            sorted_keys[size++] = (char*)kh_key(sites, k);
        }
    }

    qsort(sorted_keys, size, sizeof(char*), cmp_key);

    for (int i = 0; i < size; i++) {
        khint_t k = kh_get(str, sites, sorted_keys[i]);
        if (kh_value(sites, k).called_sites > 0) {
            char** toks = split_key((char*)kh_key(sites, k), KEY_SIZE);
            char* c = toks[0];
            char* s = toks[1];
            char* e = toks[2];
            struct site_stats site = kh_value(sites, k);
            double f = (double)site.called_sites_methylated / site.called_sites;
            printf("%s\t%s\t%s\t%d\t%d\t%d\t%.3lf\t%s\n", c, s, e, site.group_size, site.called_sites, site.called_sites_methylated, f, site.sequence);
            free(site.sequence);
            free(c);
            free(s);
            free(e);
            free(toks);
            //free((char*)kh_key(sites, k));
        }
    }


    for (khint_t k = kh_begin(sites); k != kh_end(sites); k++) {
        if (kh_exist(sites, k)) {
            free((char*)kh_key(sites, k));
        }
    }

    fclose(input);
    kh_destroy(str, sites);
    free(sorted_keys);
    free(tmp);
    return 0;
}
