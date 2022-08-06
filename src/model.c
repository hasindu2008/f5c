
/* @file model.c
**
** load a pore model from file or memory
** @author: Hasindu Gamaarachchi (hasindu@unsw.edu.au)
** @@
******************************************************************************/

#include "model.h"
#include "f5c.h"
#include "f5cmisc.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

//#define DEBUG_MODEL_PRINT 1

uint32_t eval_num_kmer(uint32_t kmer_size,uint32_t type){

    uint32_t num_kmer = 0;
    if(type==MODEL_TYPE_NUCLEOTIDE){
        num_kmer = (uint32_t)(1 << 2*kmer_size); //num_kmer should be 4^kmer_size
        assert(num_kmer <= MAX_NUM_KMER);
    }
    else if (type==MODEL_TYPE_METH){
        num_kmer = (uint32_t)pow(5,kmer_size); //num_kmer should be 5^kmer_size
        assert(num_kmer <= MAX_NUM_KMER_METH);
    }
    else {
        assert(0);
    }
    return num_kmer;
}


uint32_t read_model(model_t* model, const char* file, uint32_t type) {

    uint32_t kmer_size = MAX_KMER_SIZE;
    uint32_t num_kmer = eval_num_kmer(kmer_size,type);


    WARNING("Reading the model from file %s. This is an experimental feature. Use with caution", file);

    FILE* fp = fopen(file, "r");
    F_CHK(fp, file);

    //these two are discarded from the model. dummy vars
    char kmer[MAX_KMER_SIZE+4];
    float weight;

    //buffers for getline
    char* buffer = (char*)malloc(sizeof(char) * (100)); //READ+newline+nullcharacter
    MALLOC_CHK(buffer);
    size_t bufferSize = 100;
    ssize_t readlinebytes = 0;

    uint32_t num_k = 0;
    uint32_t line_no = 0;


    while ((readlinebytes = getline(&buffer, &bufferSize, fp)) != -1) {
        line_no++;
        if (buffer[0] == '#' ||
            strcmp(buffer, "kmer\tlevel_mean\tlevel_stdv\tsd_mean\tsd_stdv\tweight\n") ==0 ||
            strcmp(buffer, "kmer\tlevel_mean\tlevel_stdv\tsd_mean\tsd_stdv\n") ==0 ||
            strcmp(buffer,"kmer\tlevel_mean\tlevel_stdv\tsd_mean\tsd_stdv\tig_lambda\tweight\n") == 0 ||
            buffer[0] == '\n' || buffer[0] == '\r') { //comments and header
            //todo : (make generic)
            //fprintf(stderr, "%s\n", buffer);
            char key[1000];
            int32_t val=0;
            int ret=sscanf(buffer,"%s\t%d",key, &val);
            if(ret==2 && strcmp(key,"#k")==0){
                if(val<=0){
                    ERROR("k-mer size (#k\t%d) in file %s is invalid.",val,file);
                    exit(EXIT_FAILURE);
                }
                else if(val>MAX_KMER_SIZE){
                    ERROR("k-mer size (#k\t%d) in file %s larger than MAX_KMER_SIZE (%d).",val,file,MAX_KMER_SIZE);
                    exit(EXIT_FAILURE);
                }
                else{
                    INFO("k-mer size in file %s is %d",file,val);
                    kmer_size=val;
                    num_kmer = eval_num_kmer(kmer_size,type);
                }
            }
            continue;
        } else {
            //as sd_mean and sd_stdv seems not to be used just read to the dummy weight
            #ifdef LOAD_SD_MEANSSTDV
                int32_t ret =
                    sscanf(buffer, "%s\t%f\t%f\t%f\t%f\t%f", kmer,
                        &model[num_k].level_mean, &model[num_k].level_stdv,
                        &model[num_k].sd_mean, &model[num_k].sd_stdv, &weight);
            #else
                int32_t ret =
                    sscanf(buffer, "%s\t%f\t%f\t%f\t%f", kmer,
                        &model[num_k].level_mean, &model[num_k].level_stdv,
                        &weight, &weight);
            #endif
            #ifdef CACHED_LOG
                model[num_k].level_log_stdv=log(model[num_k].level_stdv);
            #endif
            num_k++;
            if (ret != 5) {
                ERROR("File %s is corrupted at line %d. Does the format adhere to examples at test/r9-models?", file, line_no);
            }
            if (num_k > num_kmer) {
                ERROR("File %s has too many entries. Expected %d kmers in the "
                      "model, but file had more than that",
                      file, num_kmer);
                exit(EXIT_FAILURE);
            }
        }

    }

    if (num_k != num_kmer) {
        ERROR("File %s prematurely ended. Expected %d kmers in the model, but "
              "file had only %d",
              file, num_kmer, num_k);
        exit(EXIT_FAILURE);
    }

#ifdef DEBUG_MODEL_PRINT
    uint32_t i = 0;
    fprintf(stderr, "level_mean\tlevel_stdv\tsd_mean\tsd_stdv\n");
    for (i = 0; i < num_kmer; i++) {
        fprintf(stderr, "%f\t%f\t%f\t%f\n", model[i].level_mean,
                model[i].level_stdv, 0.0,0.0);
    }
#endif

    free(buffer);
    fclose(fp);

    return kmer_size;
}


uint32_t set_model(model_t* model, uint32_t model_id) {

    uint32_t kmer_size=0;
    uint32_t num_kmer=0;
    float *inbuilt_model=NULL;

    if(model_id==MODEL_ID_DNA_NUCLEOTIDE){
        kmer_size=6;
        num_kmer=4096;
        inbuilt_model=r9_4_450bps_nucleotide_6mer_template_model_builtin_data;
        assert(num_kmer == (uint32_t)(1 << 2*kmer_size)); //num_kmer should be 4^kmer_size
    }
    else if(model_id==MODEL_ID_DNA_CPG){
        kmer_size=6;
        num_kmer=15625;
        inbuilt_model=r9_4_450bps_cpg_6mer_template_model_builtin_data;
        assert(num_kmer == (uint32_t)pow(5,kmer_size)); //num_kmer should be 5^kmer_size
    }
    else if(model_id==MODEL_ID_RNA_NUCLEOTIDE){
        kmer_size=5;
        num_kmer=1024;
        inbuilt_model=r9_4_70bps_u_to_t_rna_5mer_template_model_builtin_data;
        assert(num_kmer == (uint32_t)(1 << 2*kmer_size)); //num_kmer should be 4^kmer_size
    }
    else{
        assert(0);
    }

    uint32_t i = 0;
    for (i = 0; i < num_kmer; i++) {
        model[i].level_mean = inbuilt_model[i * 4 + 0];
        model[i].level_stdv = inbuilt_model[i * 4 + 1];
    #ifdef LOAD_SD_MEANSSTDV
        model[i].sd_mean = inbuilt_model[i * 4 + 2];
        model[i].sd_stdv = inbuilt_model[i * 4 + 3];
    #endif
    #ifdef CACHED_LOG
        model[i].level_log_stdv=log(model[i].level_stdv);
    #endif
    }

#ifdef DEBUG_MODEL_PRINT
    i = 0;
    fprintf(stderr, "level_mean\tlevel_stdv\tsd_mean\tsd_stdv\n");
    for (i = 0; i < num_kmer; i++) {
        fprintf(stderr, "%f\t%f\t%f\t%f\n", model[i].level_mean,
                model[i].level_stdv, 0.0, 0.0);
    }
#endif

    return kmer_size;
}
