#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TSV_HEADER_LENGTH 200

FILE *f0,*f1,*fout;
char file_name0[] = "/home/shan/CLionProjects/meth_freq_mapreducer/freq0.tsv";
char file_name1[] = "/home/shan/CLionProjects/meth_freq_mapreducer/freq1.tsv";

char file_name_out[] = "/home/shan/CLionProjects/meth_freq_mapreducer/merged.txt";

char* buf = NULL;
char* max_buf = NULL;
int max,line,file;


int read_file(int file_no){
    FILE * fp = (file_no==0)?f0:f1;
    buf = NULL;
    size_t buf_size = 0;
    if (getline(&buf, &buf_size, fp) == -1) {
        if(buf_size>0){
            free(buf);
        }
        return -1;
    }
//    fprintf(stdout,"buf %s",buf);
    char chromosome[100];
    int start;
    sscanf( buf, "%s %d", chromosome, &start);
    return start;
}

void decide(){
    int line0 = read_file(0);
    max_buf = malloc((strlen(buf) + 1) * sizeof(char));
    strcpy(max_buf,buf);
    int line1 = read_file(1);
    file = (line0>line1)?1:0;
    max = (line0>line1)?line0:line1;
    line = (line0>line1)?line1:line0;
    if(line0<line1){
        char* temp_buf = malloc((strlen(buf) + 1) * sizeof(char));
        strcpy(temp_buf,buf);
        free(buf);
        buf = malloc((strlen(max_buf) + 1) * sizeof(char));
        strcpy(buf,max_buf);
        free(max_buf);
        max_buf = malloc((strlen(temp_buf) + 1) * sizeof(char));
        strcpy(max_buf,temp_buf);
        free(temp_buf);
    }
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
    decide();

    fprintf(fout,"chromosome\tstart\tend\tnum_cpgs_in_group\tcalled_sites\tcalled_sites_methylated\tmethylated_frequency\tgroup_sequence\n");
    while(1){
        while (line>0 && max > line ){
            fprintf(fout,"%s",buf);
//            fprintf(fout,"%d %s\n",line,buf);
            line = read_file(file);
        }
        if(max!=line){
            if(line < 0){
                fprintf(fout,"%s","\n"); //end of a file. a new line is should be added.
            }
            fprintf(fout,"%s",max_buf);
            max = line;
            free(max_buf);
            max_buf = malloc((strlen(buf) + 1) * sizeof(char));
            strcpy(max_buf,buf);

            file = (file==0)?1:0;
            if(line < 0){
                break;
            }
            line = read_file(file);
        }
        else{
            char chromosome[100],group_sequence[100];
            int start,end,num_cpgs_in_group;
            double methylated_frequency;
            int called_sites_buf,called_sites_max_buf,called_sites_methylated_buf,called_sites_methylated_max_buf;
            sscanf( buf, "%s %d %d %d %d %d %f %s", chromosome, &start, &end, &num_cpgs_in_group, &called_sites_buf, &called_sites_methylated_buf, &methylated_frequency, group_sequence);
            sscanf( max_buf, "%s %d %d %d %d %d %f %s", chromosome, &start, &end, &num_cpgs_in_group, &called_sites_max_buf, &called_sites_methylated_max_buf, &methylated_frequency, group_sequence);
            methylated_frequency = (double)(called_sites_methylated_buf+called_sites_methylated_max_buf)/(called_sites_buf+called_sites_max_buf);
            fprintf(stdout,"max == line = %d and meth_freq = %.3lf\n",line,methylated_frequency);
            fprintf(fout,"%s\t%d\t%d\t%d\t%d\t%d\t%.3lf\t%s\n",chromosome,start,end,num_cpgs_in_group,called_sites_buf+called_sites_max_buf,called_sites_methylated_buf+called_sites_methylated_max_buf,methylated_frequency,group_sequence);
            decide();
            if(line < 0){
                break;
            }
        }

    }
    line = read_file(file);
    while(line>0){
        fprintf(fout,"%s",buf);
        line = read_file(file);
    }
    free(max_buf);
    fclose(f0);
    fclose(f1);
    fclose(fout);
    return 0;
}
