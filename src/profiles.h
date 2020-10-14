/* @file profiles.h
**
** profiles for various systems for better performance
** @author: David Hyland
** @author: Hasindu Gamaarachchi (hasindu@unsw.edu.au)
** @@
******************************************************************************/

#ifndef PROFILES_H
#define PROFILES_H

#include <stdint.h>

typedef struct{
    float cuda_max_readlen; // max-lf
    float cuda_avg_events_per_kmer; // avg-epk
    float cuda_max_events_per_kmer; // max-epk
    int32_t batch_size; // K
    int64_t batch_size_bases; // B
    int32_t num_thread; // t
    int64_t ultra_thresh; // ultra-thresh
    int32_t num_iop; //iop
    int8_t disable_cuda;
} parameters;

//6 core CPU, 8 GB integrated RAM
parameters jetson_tx2 = {
    .cuda_max_readlen = 3.0,
    .cuda_avg_events_per_kmer = 2.0,
    .cuda_max_events_per_kmer = 5.0,
    .batch_size = 512,
    .batch_size_bases = 2350000,
    .num_thread = 6,
    .ultra_thresh = 100000,
    .num_iop = 1,
    .disable_cuda=0
};

//4 core CPU, 4 GB integrated RAM
parameters jetson_nano = {
    .cuda_max_readlen = 3.0,
    .cuda_avg_events_per_kmer = 2.0,
    .cuda_max_events_per_kmer = 5.0,
    .batch_size = 200,
    .batch_size_bases = 1400000,
    .num_thread = 4,
    .ultra_thresh = 100000,
    .num_iop = 1,
    .disable_cuda=0
};

//8 core CPU, 16 GB integrated RAM
parameters jetson_xavier = {
    .cuda_max_readlen = 3.0,
    .cuda_avg_events_per_kmer = 2.0,
    .cuda_max_events_per_kmer = 6.25,
    .batch_size = 1024,
    .batch_size_bases = 4700000,
    .num_thread = 8,
    .ultra_thresh = 100000,
    .num_iop = 2,
    .disable_cuda=0
};

//GPU with 4GB RAM, 12 core CPU, 16GB RAM
parameters laptop_high = {
    .cuda_max_readlen = 5.0,
    .cuda_avg_events_per_kmer = 2.0,
    .cuda_max_events_per_kmer = 5.0,
    .batch_size = 512,
    .batch_size_bases = 2500000,
    .num_thread = 12,
    .ultra_thresh = 100000,
    .num_iop = 2,
    .disable_cuda=0
};

//GPU with 3GB RAM, 8 core CPU, 8GB RAM
parameters laptop_mid = {
    .cuda_max_readlen = 5.0,
    .cuda_avg_events_per_kmer = 2.0,
    .cuda_max_events_per_kmer = 5.0,
    .batch_size = 350,
    .batch_size_bases = 2000000,
    .num_thread = 8,
    .ultra_thresh = 100000,
    .num_iop = 2,
    .disable_cuda=0
};

//GPU with 2GB RAM, 4 core CPU, 4GB RAM
parameters laptop_low = {
    .cuda_max_readlen = 5.0,
    .cuda_avg_events_per_kmer = 2.0,
    .cuda_max_events_per_kmer = 5.0,
    .batch_size = 256,
    .batch_size_bases = 1500000,
    .num_thread = 4,
    .ultra_thresh = 100000,
    .num_iop = 1,
    .disable_cuda=0
};

//GPU with 12GB RAM, 16 core CPU, 64 GB RAM
parameters desktop_high = {
    .cuda_max_readlen = 5.0,
    .cuda_avg_events_per_kmer = 2.0,
    .cuda_max_events_per_kmer = 5.0,
    .batch_size = 1024,
    .batch_size_bases = 7500000,
    .num_thread = 16,
    .ultra_thresh = 100000,
    .num_iop = 6,
    .disable_cuda=0
};

//GPU with 10GB RAM, 12 core CPU, 32 GB RAM
parameters desktop_mid = {
    .cuda_max_readlen = 5.0,
    .cuda_avg_events_per_kmer = 2.0,
    .cuda_max_events_per_kmer = 5.0,
    .batch_size = 768,
    .batch_size_bases = 6250000,
    .num_thread = 12,
    .ultra_thresh = 100000,
    .num_iop = 4,
    .disable_cuda=0
};

//GPU with 8GB RAM, 8 core CPU, 32 GB RAM
parameters desktop_low = {
    .cuda_max_readlen = 5.0,
    .cuda_avg_events_per_kmer = 2.0,
    .cuda_max_events_per_kmer = 5.0,
    .batch_size = 512,
    .batch_size_bases = 5000000,
    .num_thread = 8,
    .ultra_thresh = 100000,
    .num_iop = 2,
    .disable_cuda=0
};


//GPU wuth 40GB RAM (eg: Tesla A100), 64 core CPU, 384 GB RAM
parameters hpc_high = {
    .cuda_max_readlen = 5.0,
    .cuda_avg_events_per_kmer = 2.0,
    .cuda_max_events_per_kmer = 5.0,
    .batch_size = 2560,
    .batch_size_bases = 25000000,
    .num_thread = 64,
    .ultra_thresh = 100000,
    .num_iop = 64,
    .disable_cuda=0
};

//GPU with 32GB RAM (eg: Tesla V100), 48 core CPU, 256 GB RAM
parameters hpc_mid = {
    .cuda_max_readlen = 5.0,
    .cuda_avg_events_per_kmer = 2.0,
    .cuda_max_events_per_kmer = 5.0,
    .batch_size = 2048,
    .batch_size_bases = 20000000,
    .num_thread = 48,
    .ultra_thresh = 100000,
    .num_iop = 64,
    .disable_cuda=0
};

//GPU with 16GB RAM (eg: Tesla V100), 32 core CPU, 128 GB RAM
parameters hpc_low = {
    .cuda_max_readlen = 5.0,
    .cuda_avg_events_per_kmer = 2.0,
    .cuda_max_events_per_kmer = 5.0,
    .batch_size = 1024,
    .batch_size_bases = 10000000,
    .num_thread = 32,
    .ultra_thresh = 100000,
    .num_iop = 64,
    .disable_cuda=0
};

//GPU with 16GB RAM, 32 core CPU, 128 GB RAM
parameters hpc_gpu = {
    .cuda_max_readlen = 5.0,
    .cuda_avg_events_per_kmer = 2.0,
    .cuda_max_events_per_kmer = 5.0,
    .batch_size = 1024,
    .batch_size_bases = 10000000,
    .num_thread = 32,
    .ultra_thresh = 100000,
    .num_iop = 32,
    .disable_cuda=0
};

//32 core CPU, 128 GB RAM
parameters hpc_cpu = {
    .cuda_max_readlen = 5.0,
    .cuda_avg_events_per_kmer = 2.0,
    .cuda_max_events_per_kmer = 5.0,
    .batch_size = 4096,
    .batch_size_bases = 50000000,
    .num_thread = 32,
    .ultra_thresh = 100000,
    .num_iop = 32,
    .disable_cuda=1
};

//Tesla V100-32GB. fast5 on scratch
parameters nci_gadi = {
    .cuda_max_readlen = 5.0,
    .cuda_avg_events_per_kmer = 2.0,
    .cuda_max_events_per_kmer = 5.0,
    .batch_size = 2048,
    .batch_size_bases = 20000000,
    .num_thread = 12,
    .ultra_thresh = 100000,
    .num_iop = 64,
    .disable_cuda = 0
};


#endif
