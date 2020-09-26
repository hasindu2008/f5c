/* @file profiles.c
**
** implementation of profiles for various systems for better performance
** @author: David Hyland
** @author: Hasindu Gamaarachchi (hasindu@unsw.edu.au)
** @@
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "f5c.h"
#include "f5cmisc.h"
#include "profiles.h"


//helper function to set profile to a preset value.
void set_opt_profile(opt_t *opt, parameters machine){
    opt->cuda_max_readlen = machine.cuda_max_readlen;
    opt->cuda_avg_events_per_kmer = machine.cuda_avg_events_per_kmer;
    opt->cuda_max_avg_events_per_kmer = machine.cuda_max_events_per_kmer;
    opt->batch_size = machine.batch_size;
    opt->batch_size_bases = machine.batch_size_bases;
    opt->num_thread = machine.num_thread;
    opt->ultra_thresh = machine.ultra_thresh;
    opt->num_iop = machine.num_iop;
    if(machine.disable_cuda){
        opt->flag |= F5C_DISABLE_CUDA;
    }

    fprintf(stderr, "[%s] max-lf: %.1f, avg-epk: %.1f, max-epk: %.1f, K: %d, B: %.1fM, t: %d, ultra-thresh: %.1fk, iop: %d\n", __func__,
        opt->cuda_max_readlen,opt->cuda_avg_events_per_kmer,opt->cuda_max_avg_events_per_kmer,opt->batch_size,opt->batch_size_bases/(1000.0*1000.0),opt->num_thread,
        opt->ultra_thresh/(1000.0),opt->num_iop);

}

//helper function to set user specified options. Pass -1 to corresponding arg if using default parameter value.
void set_opts(opt_t *opt, int32_t batch_size, int64_t batch_size_bases, int32_t num_thread, int64_t ultra_thresh, float cuda_max_readlen, float cuda_avg_events_per_kmer, float cuda_max_avg_events_per_kmer){
    if( batch_size != -1) opt->batch_size = batch_size;
    if( batch_size_bases != -1) opt->batch_size_bases = batch_size_bases;
    if( num_thread != -1) opt->num_thread = num_thread;
    if( ultra_thresh != -1) opt->ultra_thresh = ultra_thresh;
    if( cuda_max_readlen != -1) opt->cuda_max_readlen = cuda_max_readlen;
    if( cuda_avg_events_per_kmer != -1) opt->cuda_avg_events_per_kmer = cuda_avg_events_per_kmer;
    if( cuda_max_avg_events_per_kmer != -1) opt->cuda_max_avg_events_per_kmer = cuda_max_avg_events_per_kmer;
}


//function that sets the parameter values based on specified profile name or config file.
int set_profile(char *profile, opt_t *opt){
    //Load preset value if argument passed is the name of a machine for which there is a default profile
    if(strcmp(profile,"jetson-nano") == 0){
        set_opt_profile(opt,jetson_nano);
    }else if(strcmp(profile,"jetson-tx2") == 0){
        set_opt_profile(opt,jetson_tx2);
    }else if(strcmp(profile,"jetson-xavier") == 0){
        set_opt_profile(opt,jetson_xavier);

    }else if(strcmp(profile,"laptop-low") == 0){
        set_opt_profile(opt,laptop_low);
    }else if(strcmp(profile,"laptop-mid") == 0 || strcmp(profile,"laptop") == 0 ){
        set_opt_profile(opt,laptop_mid);
    }else if(strcmp(profile,"laptop-high") == 0){
        set_opt_profile(opt,laptop_high);

    }else if(strcmp(profile,"desktop-low") == 0){
        set_opt_profile(opt,desktop_low);
    }else if(strcmp(profile,"desktop-mid") == 0 || strcmp(profile,"desktop") == 0){
        set_opt_profile(opt,desktop_mid);
    }else if(strcmp(profile,"desktop-high") == 0){
        set_opt_profile(opt,desktop_high);

    }else if(strcmp(profile,"hpc-low") == 0){
        set_opt_profile(opt,hpc_low);
    }else if(strcmp(profile,"hpc-mid") == 0 || strcmp(profile,"hpc") == 0){
        set_opt_profile(opt,hpc_mid);
    }else if(strcmp(profile,"hpc-high") == 0){
        set_opt_profile(opt,hpc_high);


    }else if(strcmp(profile,"hpc-cpu") == 0){
        set_opt_profile(opt,hpc_cpu);
    }else if(strcmp(profile,"hpc-gpu") == 0){
        set_opt_profile(opt,hpc_gpu);

    }else if(strcmp(profile,"nci-gadi") == 0){
        set_opt_profile(opt,nci_gadi);
    }else{
        ERROR("Unknown profile %s. Trying to read profile from file.",profile);
        //Try to read from .profile file
        //profile name specifies a file from which to read values from.
        FILE *fptr = fopen(profile, "r");
        F_CHK(fptr,profile);
        int32_t batch_size, num_thread;
        int64_t batch_size_bases, ultra_thresh;
        float cuda_max_readlen,cuda_avg_events_per_kmer,cuda_max_avg_events_per_kmer;

        //read file and set parameter values
        int result = fscanf(fptr, "%f %f %f %d %" PRId64 " %d %" PRId64,
        &cuda_max_readlen,&cuda_avg_events_per_kmer,&cuda_max_avg_events_per_kmer,
        &batch_size,&batch_size_bases,&num_thread,&ultra_thresh);

        fprintf(stderr,"[%s] Profile loaded\n",__func__);
        fprintf(stderr,"[%s] batch_size: %d\nbatch_size_bases: %ld\nnum_thread: %d\nultra_thresh: %ld\ncuda_max_readlen: %f\ncuda_avg_events_per_kmer: %.2f\ncuda_max_avg_events_per_kmer: %.2f\n",
        __func__,batch_size,(long)batch_size_bases,num_thread,(long)ultra_thresh,cuda_max_readlen,cuda_avg_events_per_kmer,cuda_max_avg_events_per_kmer);

        if(result < 7){
            ERROR("%s","Malformed profile config file.");
            exit(EXIT_FAILURE);
        }

        set_opts(opt,batch_size,batch_size_bases,num_thread,ultra_thresh,cuda_max_readlen,
        cuda_avg_events_per_kmer,cuda_max_avg_events_per_kmer);
    }
    return 0;
}
