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

parameters jetson_xavier = {
    .cuda_max_readlen = 3.0,
    .cuda_avg_events_per_kmer = 2.0,
    .cuda_max_events_per_kmer = 6.25,
    .batch_size = 1024,
    .batch_size_bases = 4700000,
    .num_thread = 8,
    .ultra_thresh = 100000,
    .num_iop = 1,
    .disable_cuda=0
};

parameters laptop_high = {
    .cuda_max_readlen = 5.0,
    .cuda_avg_events_per_kmer = 2.0,
    .cuda_max_events_per_kmer = 5.0,
    .batch_size = 512,
    .batch_size_bases = 2500000,
    .num_thread = 12,
    .ultra_thresh = 100000,
    .num_iop = 1,
    .disable_cuda=0
};

parameters laptop_mid = {
    .cuda_max_readlen = 5.0,
    .cuda_avg_events_per_kmer = 2.0,
    .cuda_max_events_per_kmer = 5.0,
    .batch_size = 350,
    .batch_size_bases = 2000000,
    .num_thread = 8,
    .ultra_thresh = 100000,
    .num_iop = 1,
    .disable_cuda=0
};


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


// parameters laptop_gpu = {
//     .cuda_max_readlen = 5.0,
//     .cuda_avg_events_per_kmer = 2.0,
//     .cuda_max_events_per_kmer = 5.0,
//     .batch_size = 512,
//     .batch_size_bases = 2.5*1000*1000,
//     .num_thread = 12,
//     .ultra_thresh = 100000,
//     .num_iop = 1
// };

void set_opt_profile(opt_t *opt, parameters machine);
#endif
