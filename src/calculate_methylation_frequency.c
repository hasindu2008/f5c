#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include "khash.h"
#include "ksort.h"

#define KEY_SIZE 3

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

/* TODO needs to fine tune the comparsion logic for tuples */
#define key_lt(a, b) (strcmp(a[0], b[0]) || atoi(a[1]) < atoi(b[1]))

/* Format chromosome, start and end to colon-delimited string. */
char* make_key(char* chromosome, int start, int end) {
    int start_strlen = snprintf(NULL, 0, "%d", start);
    int end_strlen = snprintf(NULL, 0, "%d", end);
    int key_strlen = start_strlen + end_strlen + 3;
    char* key = (char*)malloc(key_strlen);
    snprintf(key, key_strlen, "%s:%d:%d", chromosome, start, end);
    return key;
}

/* Split colon-delimited keys */
char** split_key(char* key) {
    /* TODO to be free'd somewhere */
    char** tok = (char **)malloc(KEY_SIZE * sizeof(char*));
    for (int i = 0; i < KEY_SIZE; i++) {
        tok[i] = strsep(&key, ":");
    }
    return tok;
}

void update_call_stats(khash_t(str)* sites, char* key, int num_called_cpg_sites,
        bool is_methylated, char* sequence) {
    struct site_stats ss = {
        .num_reads = 0,
        .posterior_methylated = 0,
        .called_sites = 0,
        .called_sites_methylated = 0,
        .group_size = num_called_cpg_sites,
        .sequence = sequence
    };

    int absent;
    khint_t k = kh_put(str, sites, key, &absent);

    if (absent == -1) {
        fprintf(stderr, "Failed to insert key: %s\n", key);
        exit(EXIT_FAILURE);
    } else if (absent > 0) {
        kh_value(sites, k) = ss;
    }

    kh_value(sites, k).num_reads++;
    kh_value(sites, k).called_sites += num_called_cpg_sites;

    if (is_methylated) {
        kh_value(sites, k).called_sites_methylated += num_called_cpg_sites;
    }
}

struct tsv_record* get_tsv_line(FILE* fp) {
    struct tsv_record* record = (struct tsv_record*)malloc(sizeof(struct tsv_record));
    char* buf = NULL;
    size_t buf_size = 0;

    if (getline(&buf, &buf_size, fp) == -1) {
        return NULL;
    }

    record->chromosome = strsep(&buf, "\t");
    record->start = atoi(strsep(&buf, "\t"));
    record->end = atoi(strsep(&buf, "\t"));
    strsep(&buf, "\t");
    record->log_lik_ratio = strtod(strsep(&buf, "\t"), NULL);
    strsep(&buf, "\t");
    strsep(&buf, "\t");
    strsep(&buf, "\t");
    strsep(&buf, "\t");
    record->num_cpgs = atoi(strsep(&buf, "\t"));
    record->sequence = strsep(&buf, "\t");

    return record;
}

int main(int argc, char **argv) {
    FILE* input = stdin;
    double call_threshold = 2.5;
    bool split_groups = false;

    int c;

    extern char* optarg;
    extern int optind, optopt;

    while ((c = getopt(argc, argv, "c:i:s")) != -1) {
        switch(c) {
            case 'c':
                /* TODO handle overflow when calling strtod */
                call_threshold = strtod(optarg, NULL);
                break;
            case 'i': {
                          FILE* in = fopen(optarg, "r");
                          if (in == NULL) {
                              perror("Incorrect input");
                              exit(EXIT_FAILURE);
                          }
                          input = in;
                          break;
                      }
            case 's':
                      split_groups = true;
                      break;
            case ':':
                      fprintf(stderr, "Option -%c requires an operand\n", optopt);
                      exit(EXIT_FAILURE);
            case '?':
                      fprintf(stderr, "Unrecognized option: -%c\n", optopt);
                      exit(EXIT_FAILURE);
            default:
                      break;
        }
    }

    khash_t(str)* sites = kh_init(str);
    struct tsv_record* record;
    fscanf(input, "%*s\n");

    while ((record = get_tsv_line(input)) != NULL) {
        int num_sites = record->num_cpgs;
        double llr = record->log_lik_ratio;

        if (fabs(llr) < call_threshold) {
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
                const char* sg = "split_groups";
                update_call_stats(sites, key, 1, is_methylated, (char*)sg);
                char* substring = strstr(sequence + cg_pos + 1, "CG");
                cg_pos = substring == NULL ? -1 : substring - sequence;
            }
        } else {
            char* key = make_key(record->chromosome, record->start, record->end);
            update_call_stats(sites, key, num_sites, is_methylated, sequence);
        }
        free(record);
    }

    printf("chromosome\tstart\tend\tnum_cpgs_in_group\tcalled_sites\tcalled_sites_methylated\tmethylated_frequency\tgroup_sequence\n");

    char** sorted_keys = (char**)malloc(kh_size(sites) * sizeof(char*));
    /* TODO sort keys here */
    int size = 0;
    for (khint_t k = kh_begin(sites); k != kh_end(sites); k++) {
        sorted_keys[size++] = (char*)kh_key(sites, k);
    }
    for (int i = 0; i < size; i++) {
        khint_t k = kh_get(str, sites, sorted_keys[i]);
        if (kh_value(sites, k).called_sites > 0) {
            /*
             (c, s, e) = key.split(":")
             f = float(sites[key].called_sites_methylated) / sites[key].called_sites
             print("%s\t%s\t%s\t%d\t%d\t%d\t%.3f\t%s" % (c, s, e, sites[key].group_size, sites[key].called_sites, sites[key].called_sites_methylated, f, sites[key].sequence))
             */
        }
    }

    fclose(input);
    kh_destroy(str, sites);
    return 0;
}
