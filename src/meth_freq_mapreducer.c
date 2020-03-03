#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <unistd.h>
#define TSV_HEADER_LENGTH 200

//#define HAVE_EXECINFO_H

#ifdef HAVE_EXECINFO_H
#include <execinfo.h>
#endif

//make the segmentation faults a bit cool
void sig_handler(int sig) {
#ifdef HAVE_EXECINFO_H
    void* array[100];
    size_t size = backtrace(array, 100);
    fprintf(stderr,"I regret to inform that a segmentation fault occurred. But at least "
          "it is better than a wrong answer%s",
          ".");
    fprintf(stderr,
            "[%s::DEBUG]\033[1;35m Here is the backtrace in case it is of any "
            "use:\n",
            __func__);
    backtrace_symbols_fd(&array[2], size - 1, STDERR_FILENO);
    fprintf(stderr, "\033[0m\n");
#else
    fprintf(stderr,"I regret to inform that a segmentation fault occurred. But at least "
          "it is better than a wrong answer%s",
          ".");
#endif
    exit(EXIT_FAILURE);
}


FILE *f0,*f1,*fout;
char file_name0[] = "/home/shan/CLionProjects/meth_freq_mapreducer/freq00.tsv";
char file_name1[] = "/home/shan/CLionProjects/meth_freq_mapreducer/freq11.tsv";

//char file_name0[] = "/media/shan/OS/Hiruna/FYP/methcall_freq23/2meth-freq.tsv";
//char file_name1[] = "/media/shan/OS/Hiruna/FYP/methcall_freq23/3meth-freq.tsv";


char file_name_out[] = "/home/shan/CLionProjects/meth_freq_mapreducer/merged.txt";

char* buf = NULL;
char* max_buf = NULL;
int max,line,file;

char * sorted_chromosomes[] = {
        "chr1","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
        "chr2","chr20","chr21","chr22",
        "chr3","chr4","chr5","chr6","chr8","chr9",
        "chrEBV",
        "chrM",
        "chrUn_GL",
        "chrUn_KI",
        "chrX"
};
char chromosome [100];
char prev_chromosome [100];
int do_compare = 1;
char* memory_buf0 = NULL;
char* memory_buf1 = NULL;

int read_file(int file_no,int from){
    FILE * fp = (file_no==0)?f0:f1;
    if(buf!=NULL){
        int temp;
        sscanf( buf, "chr%s %d", prev_chromosome,&temp);
    }
    buf = NULL;
    size_t buf_size = 0;
    if(from == 0)buf = memory_buf0;
    else if (from == 1)buf = memory_buf1;
    else if (getline(&buf, &buf_size, fp) == -1) {
        if(buf_size>0){
            free(buf);
        }
        return -2;
    }
//    fprintf(stdout,"buf %s",buf);
    int start;
    sscanf( buf, "chr%s %d", chromosome, &start);
    int strcmp_value = strcmp(chromosome,prev_chromosome);
    return (strcmp_value==0 || do_compare==0)?start:-1;
}

void store_buf(int file){
    if(!file){
        if(memory_buf0)free(memory_buf0);
        memory_buf0 = malloc((strlen(buf) + 1) * sizeof(char));
        strcpy(memory_buf0,buf);
    }
    else{
        if(memory_buf1)free(memory_buf1);
        memory_buf1 = malloc((strlen(buf) + 1) * sizeof(char));
        strcpy(memory_buf1,buf);
    }
}

void decide(int memory){
    int line0 = read_file(0,memory+1);
    if(line0 < 0){
        line = line0;
        return;
    }
    if(max_buf!=NULL)free(max_buf);
    max_buf = malloc((strlen(buf) + 1) * sizeof(char));
    strcpy(max_buf,buf);
    int line1 = read_file(1,memory+2);
    if(line1 < 0){
        line = line1;
        return;
    }
//    if(strcmp(chromosome,prev_chromosome)){
//
//    }
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

int for_each_chromosome(){
    while(1){
        while (line>0 && max > line ){
            fprintf(fout,"%s",buf);
//            fprintf(fout,"%d %s\n",line,buf);
            line = read_file(file,-1);
        }
        if(line == -1){
            fprintf(fout,"%s",max_buf);
            file = (file==0)?1:0;
            break;
        }
        if(max!=line){
            file = (file==0)?1:0;
            if(line == -2){
                fprintf(fout,"%s","\n"); //end of a file. a new line is should be added.
            }
            fprintf(fout,"%s",max_buf);
            if(line < 0){
                break;
            }
            max = line;
            free(max_buf);
            max_buf = malloc((strlen(buf) + 1) * sizeof(char));
            strcpy(max_buf,buf);
            line = read_file(file,-1);
        }
        else{
            char group_sequence[100];
            int chromosome,start,end,num_cpgs_in_group;
            double methylated_frequency;
            int called_sites_buf,called_sites_max_buf,called_sites_methylated_buf,called_sites_methylated_max_buf;
            sscanf( buf, "chr%d %d %d %d %d %d %f %s", &chromosome, &start, &end, &num_cpgs_in_group, &called_sites_buf, &called_sites_methylated_buf, &methylated_frequency, group_sequence);
            sscanf( max_buf, "chr%d %d %d %d %d %d %f %s", &chromosome, &start, &end, &num_cpgs_in_group, &called_sites_max_buf, &called_sites_methylated_max_buf, &methylated_frequency, group_sequence);
            methylated_frequency = (double)(called_sites_methylated_buf+called_sites_methylated_max_buf)/(called_sites_buf+called_sites_max_buf);
//            fprintf(stdout,"max == line = %d and meth_freq = %.3lf\n",line,methylated_frequency);
            fprintf(fout,"chr%d\t%d\t%d\t%d\t%d\t%d\t%.3lf\t%s\n",chromosome,start,end,num_cpgs_in_group,called_sites_buf+called_sites_max_buf,called_sites_methylated_buf+called_sites_methylated_max_buf,methylated_frequency,group_sequence);
            decide(1);
            if(line < 0){
                file = (file==0)?1:0;
                break;
            }
        }

    }
    if(line != -2){
        do_compare = 0;
        store_buf((file+1)%2);
        line = read_file(file,-1);
        do_compare = 1; // not the best place to put this
        if(line != -2 && strcmp(chromosome,prev_chromosome)){
            while(line>0){
                fprintf(fout,"%s",buf);
                line = read_file(file,-1);
            }
        }
        if(line != -2){
            store_buf((file+2)%2);
            return -1;
        }
        read_file(0,(file+1)%2);
        fprintf(fout,"%s",buf);
        file = (file==0)?1:0;
    }
    do_compare = 0;
    line = read_file(file,-1);
    while(line != -2){
        fprintf(fout,"%s",buf);
        line = read_file(file,-1);
    }
    return line;
}

int main() {
    signal(SIGSEGV, sig_handler);
    f0 = fopen(file_name0, "r"); // read mode
    f1 = fopen(file_name1, "r"); // read mode
    fout = stdout;//fopen(file_name_out, "w"); // read mode
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
    int ret = 1;
    do{
        do_compare =0;
        decide(ret);
        do_compare =1;
        ret = for_each_chromosome();
//        fprintf(fout,"\nret = %d\n",ret);
    }while (ret == -1);

    free(max_buf);
    fclose(f0);
    fclose(f1);
//    fclose(fout);
    return 0;
}
