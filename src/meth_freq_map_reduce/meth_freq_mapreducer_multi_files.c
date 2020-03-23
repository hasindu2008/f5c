#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <unistd.h>
#include <assert.h>

#define TSV_HEADER_LENGTH 200
#define FILE_NAME_LENGTH 20

FILE *fout;

char (*buf)[FILE_NAME_LENGTH];

int read_line(int file_no,FILE* fp){
    strcpy(buf[file_no], NULL);
    size_t buf_size = 0;
    if (getline(&buf[file_no], &buf_size, fp) == -1) {
        if(buf_size>0){
            free(buf[file_no]);
        }
        return -1;
    }
    if(strcmp(buf[file_no],"\n")==0)return -1;
    return 1;
}

static inline void strtok_null_check(char *tmp, int64_t line_num){
    // TODO: the file no is passed. have to change
    if(tmp == NULL){
        fprintf(stderr,"Corrupted tsv file? Offending line is line number %ld (1-based index)\n",line_num);
        exit(EXIT_FAILURE);
    }
}

struct tsv_record {
    char* chromosome;
    char *group_sequence;
    int start;
    int end;
    int num_cpgs_in_group;
    int called_sites;
    int called_sites_methylated;
};

void get_tsv_line(struct tsv_record* record, int file_no, int64_t line_num) {
    char* tmp;
    char * buffer = strdup(buf[file_no]);
    //chromosome
    tmp = strtok(buffer, "\t");
    strtok_null_check(tmp,line_num);
    if(record->chromosome){
        assert(strcmp(record->chromosome,tmp) == 0); // should have the same group sequence
    } else{
        record->chromosome = strdup(tmp);
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
        fprintf(stderr,"chromosome coordinates cannot be negative. Offending line is line number %ld (1-based index)\n",line_num);
        free(buffer);
        exit(EXIT_FAILURE);
    }

    //num_cpgs_in_group
    tmp = strtok(NULL, "\t");
    strtok_null_check(tmp,line_num);;
    record->num_cpgs_in_group = atoi(tmp);
    if(record->num_cpgs_in_group<0){
        fprintf(stderr,"num_cpgs_in_group cannot be negative. Offending line is line number %ld (1-based index)\n",line_num);
        free(buffer);
        exit(EXIT_FAILURE);
    }

    //called_sites
    tmp = strtok(NULL, "\t");
    strtok_null_check(tmp,line_num);;
    record->called_sites = atoi(tmp);

    //called_sites_methylated
    tmp = strtok(NULL, "\t");
    strtok_null_check(tmp,line_num);;
    record->called_sites_methylated = atoi(tmp);

    //methylated_frequency
    tmp = strtok(NULL, "\t");
    strtok_null_check(tmp,line_num);;

    //group_sequence
    tmp = strtok(NULL, "\t\n");
    strtok_null_check(tmp,line_num);
    if(record->group_sequence){
        assert(strcmp(record->group_sequence,tmp) == 0); // should have the same group sequence
    } else{
        record->group_sequence = strdup(tmp);
    }

    if(strtok(NULL, "\t\n")!=NULL){
        fprintf(stderr,"encountered more columns than expected. Offending line is line number %ld (1-based index)\n",line_num);
        free(buffer);
        exit(EXIT_FAILURE);
    }
    free(buffer);
}

static const char *MAP_REDUCE_USAGE_MESSAGE =
        "Usage: f5c map-reduce -o [OUTPUT_FILE_NAME] -n [NO_OF_INPUT_FILES] -f [INPUT_FILE_1] [INPUT_FILE_2] ...\n";

int main(int argc, char **argv) {

    char (*inputfileNames)[FILE_NAME_LENGTH];

    char *outputFileName = NULL;
    int no_of_files = 0;
    int index;
    int c;

    while ((c = getopt (argc, argv, "o:n:f")) != -1)
        switch (c) {
            case 'o':
                outputFileName = optarg;
                break;
            case 'n':
                no_of_files = atoi(optarg);
                inputfileNames = malloc(sizeof * inputfileNames * no_of_files);
                break;
            case 'f':
                for(index = 0; optind < argc && *argv[optind] != '-'; optind++, index++){
                    strcpy(inputfileNames[index], argv[optind]);
                }
                if (index != no_of_files){
                    fprintf (stderr, "Number of input files and input file names mismatch\n");
                    return 1;
                }
                break;
            default:
                fprintf (stderr, "%s", MAP_REDUCE_USAGE_MESSAGE);
                return 1;
        }

    for (int i = 0; i < no_of_files; i++){
        printf("file %d = %s\n", i, inputfileNames[i]);
    }

    FILE* file_pointers [no_of_files];
    char tmp[TSV_HEADER_LENGTH];
    char header[] = {"chromosome\tstart\tend\tnum_cpgs_in_group\tcalled_sites\tcalled_sites_methylated\tmethylated_frequency\tgroup_sequence\n"};
    for(int i = 0; i<no_of_files; i++){
        file_pointers[i] = fopen(inputfileNames[i], "r"); // read mode
        if (file_pointers[i] == NULL){
            perror("Error while opening the file.\n");
            exit(EXIT_FAILURE);
        }
        /* ignore header */
        char *ret0=fgets(tmp, TSV_HEADER_LENGTH, file_pointers[i]);
        if(ret0==-1 || strcmp(tmp,header)!=0){
            fprintf(stderr,"Bad file format with no header?\n");
            exit(1);
        }
    }
    fout = fopen(outputFileName, "w"); // read mode
    if (fout == NULL){
        perror("Error while opening the file.\n");
        exit(EXIT_FAILURE);
    }

    fprintf(fout,"%s",header);

    int starts[no_of_files];
    int lines[no_of_files];
    char * chromosome [no_of_files];
    int active_file = -1;
    for(int i = 0; i<no_of_files; i++){
        lines[i] = read_line(i,file_pointers[i]);
    }
    for(int i = 0; i<no_of_files; i++){
        if(lines[i] == 1){
            active_file = i;
            break;
        }
    }
    while(active_file!=-1){
        for(int i = 0; i<no_of_files; i++){
            if(lines[i]!=-1){
//                sscanf( buf[i], "%s\t%d", chromosome[i], &starts[i]);
                char* tmp;
                char * buffer = strdup(buf[i]);
                //chromosome
                tmp = strtok(buffer, "\t");
                strtok_null_check(tmp,i);
                chromosome[i] = strdup(tmp);

                //start
                tmp = strtok(NULL, "\t");
                strtok_null_check(tmp,i);;
                starts[i] = atoi(tmp);

                free(buffer);

            }
        }
        int lexi_idx [no_of_files];
        for(int i = 0; i<no_of_files; i++){
            lexi_idx[i] = -1;
        }
        int min_lexi_idx = active_file;
        for(int i = min_lexi_idx+1; i<no_of_files; i++){
            if(lines[i]==-1){
                continue;
            }
            int lexi = strcmp(chromosome[min_lexi_idx], chromosome[i]);
            if(lexi > 0){
                min_lexi_idx = i;
            } else if (lexi == 0){
                int temp = min_lexi_idx;
                min_lexi_idx = i;
                lexi_idx[min_lexi_idx] = temp;
            }
        }
        if(lexi_idx[min_lexi_idx] == -1){
            fprintf(fout,"%s",buf[min_lexi_idx]);
            lines[min_lexi_idx] = read_line(min_lexi_idx,file_pointers[min_lexi_idx]);
        } else{
            assert(lexi_idx[min_lexi_idx] < no_of_files);
            int start_comp[no_of_files];
            for(int i = 0; i<no_of_files; i++){
                start_comp[i] = -1;
            }
            int next = min_lexi_idx;
            int min_start_idx = next;
            int min_start = starts[min_start_idx];
            next = lexi_idx[next];
            while(next!=-1){
                if(min_start > starts[next]){
                    min_start = starts[next];
                    min_start_idx = next;
                } else if(min_start == starts[next]){
                    start_comp[next] = min_start_idx;
                    min_start_idx = next;
                }
                next = lexi_idx[next];
            }
            if(start_comp[min_start_idx] == -1){
                fprintf(fout,"%s",buf[min_start_idx]);
                lines[min_start_idx] = read_line(min_start_idx,file_pointers[min_start_idx]);
            } else{
                struct tsv_record* rcd = malloc(sizeof(struct tsv_record));
                rcd->chromosome = NULL;
                rcd->group_sequence = NULL;
                if(rcd == NULL){
                    perror("malloc failed\n");
                    exit(EXIT_FAILURE);
                }
                double methylated_frequency;
                int called_sites = 0;
                int called_sites_methylated = 0;
                next = min_start_idx;
                assert(start_comp[next] != -1); //there should be at least two same start values
                while(next!=-1){
                    assert(starts[next] == min_start);
                    get_tsv_line(rcd,next,next);
//                    sscanf( buf[next], "%s\t%d\t%d\t%d\t%d\t%d\t%f\t%s", chromosome[next], &starts[next], &end, &num_cpgs_in_group, &temp_called_sites, &temp_called_sites_methylated, &methylated_frequency, group_sequence);
                    called_sites += rcd->called_sites;
                    called_sites_methylated += rcd->called_sites_methylated;
                    lines[next] = read_line(next,file_pointers[next]);
                    next = start_comp[next];
                }
                methylated_frequency = (double)(called_sites_methylated)/(called_sites);
                assert(chromosome[min_start_idx == rcd->chromosome]);
                assert(starts[min_start_idx == rcd->start]);
                fprintf(fout,"%s\t%d\t%d\t%d\t%d\t%d\t%.3lf\t%s\n",rcd->chromosome,rcd->start,rcd->end,rcd->num_cpgs_in_group,called_sites,called_sites_methylated,methylated_frequency,rcd->group_sequence);
                free(rcd->chromosome);
                free(rcd->group_sequence);
                free(rcd);
            }
        }
        active_file = -1;
        for(int i = 0; i<no_of_files; i++){
            if(lines[i] == 1){
                active_file = i;
                break;
            }
        }
    }
    for(int i = 0; i<no_of_files; i++){
        fclose(file_pointers[i]);
    }
//    fclose(fout);
    return 0;
}
