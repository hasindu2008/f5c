#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <unistd.h>
#define TSV_HEADER_LENGTH 200

FILE *f0,*f1,*fout;
char file_name0[] = "/home/shan/CLionProjects/meth_freq_mapreducer/freq00.tsv";
char file_name1[] = "/home/shan/CLionProjects/meth_freq_mapreducer/freq11.tsv";
char *files[] = {file_name0,file_name1};
//char file_name0[] = "/media/shan/OS/Hiruna/FYP/methcall_freq23/2meth-freq-split.tsv";
//char file_name1[] = "/media/shan/OS/Hiruna/FYP/methcall_freq23/3meth-freq-split.tsv";


char file_name_out[] = "/home/shan/CLionProjects/meth_freq_mapreducer/merged-split.txt";

char* buf[] = {NULL,NULL};


int read_line(int file_no,FILE* fp){
//    FILE * fp = (file_no==0)?f0:f1;
    buf[file_no] = NULL;
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


int main() {
    int no_of_files = 2;
    FILE* file_pointers [no_of_files];
    char tmp[TSV_HEADER_LENGTH];
    for(int i = 0; i<no_of_files; i++){
        file_pointers[i] = fopen(files[i], "r"); // read mode
        if (file_pointers[i] == NULL){
            perror("Error while opening the file.\n");
            exit(EXIT_FAILURE);
        }
        /* ignore header */
        char *ret0=fgets(tmp, TSV_HEADER_LENGTH, file_pointers[i]);
        if(ret0==NULL){
            fprintf(stderr,"Bad file format with no header?\n");
            exit(1);
        }
    }
    fout = stdout;//fopen(file_name_out, "w"); // read mode
    if (fout == NULL){
        perror("Error while opening the file.\n");
        exit(EXIT_FAILURE);
    }

    fprintf(fout,"chromosome\tstart\tend\tnum_cpgs_in_group\tcalled_sites\tcalled_sites_methylated\tmethylated_frequency\tgroup_sequence\n");

    int starts[no_of_files];
    int lines[no_of_files];
    char chromosome [no_of_files][100];
    int max_line = -1;
    for(int i = 0; i<no_of_files; i++){
        lines[i] = read_line(i,file_pointers[i]);
        if(max_line<lines[i]){
            max_line = lines[i];
        }
    }
//    int line0 = read_line(0);
//    int line1 = read_line(1);

    while(max_line!=-1){
        for(int i = 0; i<no_of_files; i++){
            if(lines[i]){
                sscanf( buf[i], "%s %d", chromosome[i], &starts[i]);
            }
        }
        int lexi_idx [no_of_files];
        for(int i = 0; i<no_of_files; i++){
            lexi_idx[i] = -1;
        }
        int min_lexi_idx = 0;
        for(int i = 1; i<no_of_files; i++){
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
                char group_sequence[100];
                int end,num_cpgs_in_group;
                double methylated_frequency;
                int called_sites = 0;
                int called_sites_methylated = 0;
                next = min_start_idx;

                while(next!=-1){
                    int temp_called_sites,temp_called_sites_methylated;
                    sscanf( buf[next], "%s %d %d %d %d %d %f %s", chromosome[next], &starts[next], &end, &num_cpgs_in_group, &temp_called_sites, &temp_called_sites_methylated, &methylated_frequency, group_sequence);
                    called_sites += temp_called_sites;
                    called_sites_methylated += temp_called_sites_methylated;
                    lines[next] = read_line(next,file_pointers[next]);
                    next = start_comp[next];
                }
                methylated_frequency = (double)(called_sites_methylated)/(called_sites);
                fprintf(fout,"%s\t%d\t%d\t%d\t%d\t%d\t%.3lf\t%s\n",chromosome[min_start_idx],starts[min_start_idx],end,num_cpgs_in_group,called_sites,called_sites_methylated,methylated_frequency,group_sequence);
            }
        }
        max_line = -1;
        for(int i = 0; i<no_of_files; i++){
            if(max_line<lines[i]){
                max_line = lines[i];
            }
        }
    }
    for(int i = 0; i<no_of_files; i++){
        fclose(file_pointers[i]);
    }
//    fclose(fout);
    return 0;
}
