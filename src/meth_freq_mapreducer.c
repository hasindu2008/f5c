#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <unistd.h>
#define TSV_HEADER_LENGTH 200

FILE *f0,*f1,*fout;
//char file_name0[] = "/home/shan/CLionProjects/meth_freq_mapreducer/freq00.tsv";
//char file_name1[] = "/home/shan/CLionProjects/meth_freq_mapreducer/freq11.tsv";

char file_name0[] = "/media/shan/OS/Hiruna/FYP/methcall_freq23/2meth-freq.tsv";
char file_name1[] = "/media/shan/OS/Hiruna/FYP/methcall_freq23/3meth-freq.tsv";


char file_name_out[] = "/home/shan/CLionProjects/meth_freq_mapreducer/merged.txt";

char* buf[] = {NULL,NULL};


int read_line(int file_no){
    FILE * fp = (file_no==0)?f0:f1;
    buf[file_no] = NULL;
    size_t buf_size = 0;
    if (getline(&buf[file_no], &buf_size, fp) == -1) {
        if(buf_size>0){
            free(buf[file_no]);
        }
        return -1;
    }
    return 1;
}


int main() {
    f0 = fopen(file_name0, "r"); // read mode
    f1 = fopen(file_name1, "r"); // read mode
    fout = fopen(file_name_out, "w"); // read mode
    if (f0 == NULL || f1 == NULL || fout == NULL){
        perror("Error while opening the file.\n");
        exit(EXIT_FAILURE);
    }
    /* ignore header */
    char tmp[TSV_HEADER_LENGTH];
    char *ret0=fgets(tmp, TSV_HEADER_LENGTH, f0);
    char *ret1=fgets(tmp, TSV_HEADER_LENGTH, f1);
    if(ret0==NULL || ret1==NULL){
        fprintf(stderr,"Bad file format with no header?\n");
        exit(1);
    }

    fprintf(fout,"chromosome\tstart\tend\tnum_cpgs_in_group\tcalled_sites\tcalled_sites_methylated\tmethylated_frequency\tgroup_sequence\n");

    int start[] = {0,0};
    char chromosome [2][100];
    int line0 = read_line(0);
    int line1 = read_line(1);
    while(line0 >0 && line1 > 0){
        sscanf( buf[0], "chr%s %d", chromosome[0], &start[0]);
        sscanf( buf[1], "chr%s %d", chromosome[1], &start[1]);
        int lexi = strcmp(chromosome[0], chromosome[1]);
        if(lexi < 0){
            fprintf(fout,"%s",buf[0]);
            line0 = read_line(0);
        } else if (lexi > 0){
            fprintf(fout,"%s",buf[1]);
            line1 = read_line(1);
        } else{
            if(start[0]<start[1]){
                fprintf(fout,"%s",buf[0]);
                line0 = read_line(0);
            } else if (start[0]>start[1]){
                fprintf(fout,"%s",buf[1]);
                line1 = read_line(1);
            } else{
                char group_sequence[100];
                int end,num_cpgs_in_group;
                double methylated_frequency;
                int called_sites_buf,called_sites_max_buf,called_sites_methylated_buf,called_sites_methylated_max_buf;
                sscanf( buf[0], "chr%s %d %d %d %d %d %f %s", chromosome[0], &start[0], &end, &num_cpgs_in_group, &called_sites_buf, &called_sites_methylated_buf, &methylated_frequency, group_sequence);
                sscanf( buf[1], "chr%s %d %d %d %d %d %f %s", chromosome[1], &start[1], &end, &num_cpgs_in_group, &called_sites_max_buf, &called_sites_methylated_max_buf, &methylated_frequency, group_sequence);
                methylated_frequency = (double)(called_sites_methylated_buf+called_sites_methylated_max_buf)/(called_sites_buf+called_sites_max_buf);
                fprintf(fout,"chr%s\t%d\t%d\t%d\t%d\t%d\t%.3lf\t%s\n",chromosome[0],start[0],end,num_cpgs_in_group,called_sites_buf+called_sites_max_buf,called_sites_methylated_buf+called_sites_methylated_max_buf,methylated_frequency,group_sequence);
                line0 = read_line(0);
                line1 = read_line(1);

            }
        }

    }
//    fprintf(fout,"file break\n");
    while (line0 > 0){
        fprintf(fout,"%s",buf[0]);
        line0 = read_line(0);
    }
    while (line1 > 0){
        fprintf(fout,"%s",buf[1]);
        line1 = read_line(1);
    }

    fclose(f0);
    fclose(f1);
    fclose(fout);
    return 0;
}
