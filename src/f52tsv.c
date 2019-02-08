#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "f5c.h"
#include "f5cmisc.h"
#include "error.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <dirent.h>




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


int fast5_to_tsv(std::string fast5_path_str){

    char* fast5_path =
        (char*)malloc(fast5_path_str.size() + 10); // is +10 needed? do errorcheck
    strcpy(fast5_path, fast5_path_str.c_str());

    fast5_t f5;

    hid_t hdf5_file = fast5_open(fast5_path);
    if (hdf5_file >= 0) {
        int32_t ret=fast5_read(hdf5_file, &f5);
        if(ret<0){
            WARNING("Fast5 file [%s] is unreadable and will be skipped", fast5_path_str.c_str());
            free(fast5_path);
            return 0;
        }
        fast5_close(hdf5_file);

        printf(">%s\tPATH:%s\tLN:%llu\n", "aaaa", fast5_path,
                f5.nsample);
            uint32_t j = 0;
        for (j = 0; j < f5.nsample; j++) {
            printf("%d\t", (int)f5.rawptr[j]);
        }
        printf("\n");
    }
    else{
        WARNING("Fast5 file [%s] is unreadable and will be skipped", fast5_path_str.c_str());
        free(fast5_path);
        return 0;     
    }
    free(f5.rawptr);
    free(fast5_path);
    return 1;

}

//adapted from nanopolish index
// recurse through a directory and handle all fast5 on the way
void recurse_dir(const std::string& path)
{
    fprintf(stderr, "[readdb] indexing %s\n", path.c_str());
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
                fast5_to_tsv(full_fn);    
            }
        }
    }
} 


int f52tsv_main(int argc, char** argv){

    //bool die = false;

    if (argc < 1) {
        std::cerr << "[f5c f52tsv] at least one folder is needed arguments\n";
        return 1;
        //die = true;
    }

    int i=0;
    for(i=0;i<argc;i++){
        std::string path = argv[i];
        recurse_dir(path);
    }

    // if (die)
    // {
    //     std::cout << "\n" << INDEX_USAGE_MESSAGE;
    //     exit(EXIT_FAILURE);
    // }    
    
    return 0;
}