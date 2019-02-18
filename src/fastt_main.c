#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "fast5lite.h"
#include "f5cmisc.h"
#include "error.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <dirent.h>
#include <getopt.h>
#include <unistd.h>
#include "ftidx.h"

static double realtime0=0;
static int64_t bad_fast5_file=0;
static int64_t total_reads=0;

// ref: http://stackoverflow.com/a/612176/717706
// return true if the given name is a directory
static bool is_directory(const std::string& file_name){
    auto dir = opendir(file_name.c_str());
    if(not dir) {
        return false;
    }
    closedir(dir);
    return true;
}

//from nanopolish idnex
static std::vector< std::string > list_directory(const std::string& file_name)
{
    std::vector< std::string > res;
    DIR* dir;
    struct dirent *ent;

    dir = opendir(file_name.c_str());
    if(not dir) {
        return res;
    }
    while((ent = readdir(dir)) != nullptr) {
        res.push_back(ent->d_name);
    }
    closedir(dir);
    return res;
}


int fast5_to_fastt(std::string fast5_path_str){

    total_reads++;

    char* fast5_path =
        (char*)malloc(fast5_path_str.size() + 10); // is +10 needed? do errorcheck
    strcpy(fast5_path, fast5_path_str.c_str());

    fast5_t f5;

    hid_t hdf5_file = fast5_open(fast5_path);
    if (hdf5_file >= 0) {
        int32_t ret=fast5_read(hdf5_file, &f5);
        if(ret<0){
            WARNING("Fast5 file [%s] is unreadable and will be skipped", fast5_path_str.c_str());
            bad_fast5_file++;
            fast5_close(hdf5_file);
            free(fast5_path);
            return 0;
        }
        std::string read_id = fast5_get_read_id(hdf5_file);
        if (read_id==""){
            WARNING("Fast5 file [%s] does not have a read ID and will be skipped", fast5_path_str.c_str());
            bad_fast5_file++;
            fast5_close(hdf5_file);
            free(fast5_path);   
            return 0;         
        }
        fast5_close(hdf5_file);

        //printf("@read_id\tn_samples\tdigitisation\toffset\trange\tsample_rate\traw_signal\tnum_bases\tsequence\nfast5_path");

        printf("%s\t%lld\t%.1f\t%.1f\t%.1f\t%.1f\t", read_id.c_str(), 
                f5.nsample,f5.digitisation, f5.offset, f5.range, f5.sample_rate);
        uint32_t j = 0;
        for (j = 0; j < f5.nsample-1; j++) {
            printf("%d,", (int)f5.rawptr[j]);
        }
        if(j<f5.nsample){
            printf("%d", (int)f5.rawptr[j]);
            j++;
        }
        printf("\t%d\t%s\t%s\n",0,".",fast5_path);
    }
    else{
        WARNING("Fast5 file [%s] is unreadable and will be skipped", fast5_path_str.c_str());
        bad_fast5_file++;
        free(fast5_path);
        return 0;     
    }
    free(f5.rawptr);
    free(fast5_path);
    return 1;

}

// adapted from nanopolish index
// recurse through a directory and extract all fast5 on the way
void recurse_dir(const std::string& path)
{
    fprintf(stderr, "[%s::%.3f*%.2f] Extracting fast5 from %s\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0),path.c_str());

    if (is_directory(path)) {
        auto dir_list = list_directory(path);
        for (const auto& fn : dir_list) {
            if(fn == "." or fn == "..") {
                continue;
            }

            std::string full_fn = path + "/" + fn;
            if(is_directory(full_fn)) {
                // recurse
                recurse_dir(full_fn);
            } else if (full_fn.find(".fast5") !=  std::string::npos) {
                //write fast5
                fast5_to_fastt(full_fn);    
            }
        }
    }
} 


int fastt_main(int argc, char** argv){


    realtime0 = realtime();
    FILE *fp_help = stderr;
    char *tsvfile=NULL;

    int c;
    while ((c = getopt(argc, argv, "i:h")) != -1) {
        switch(c) {
            case 'h':
                fp_help=stdout;
                break;
            case 'i': 
                tsvfile=optarg;
                if(tsvfile==NULL){
                    fprintf(stderr,"File name cannot be empty\n");
                    exit(EXIT_FAILURE);
                }
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


    if(tsvfile!=NULL){
        if(argc-optind >= 1){
            ftidx_t *ftidx= fti_load(tsvfile);
            if (ftidx == NULL) {
                fprintf(stderr, "Error loading ftidx for %s\n",
                        tsvfile);
                exit(EXIT_FAILURE);
            }
            int i;
            for(i=optind;i<argc;i++){
                int len=0;
                fprintf(stderr, "Fetching %s\n",
                        argv[i]);
                char *record = fti_fetch(ftidx, argv[i], &len);
                if(record==NULL || len <0){
                    fprintf(stderr, "Error locating %s\n",
                        argv[i]);
                    exit(EXIT_FAILURE);  
                }
                printf("%s\n",record);
            }
            fti_destroy(ftidx);        
        }
        else{
            // build ftidx
            int ret = fti_build(tsvfile);
            if (ret != 0) {
                fprintf(stderr, "Error running ftidx_build on %s\n",
                        tsvfile);
                exit(EXIT_FAILURE);
            }
        }
    }
    else{

        if (argc-optind < 1 || fp_help==stdout) {
            fprintf(fp_help,
                "Usage: f5c fastt [OPTIONS] [fast5_dir] [...]\n");        
            if(fp_help == stdout){
                exit(EXIT_SUCCESS);
            }
            exit(EXIT_FAILURE);
        }

        printf("#read_id\tn_samples\tdigitisation\toffset\trange\tsample_rate\traw_signal\tnum_bases\tsequence\tfast5_path\n");

        int i;
        for(i=optind;i<argc;i++){
            std::string path = argv[i];
            recurse_dir(path);
        }

        fprintf(stderr, "\n[%s] total reads: %ld, bad fast5: %ld",
                __func__,total_reads,bad_fast5_file);
        //fprintf(stderr,"\n[%s] total bases: %.1f Mbases",__func__,core->sum_bases/(float)(1000*1000));
    }


    return 0;
}