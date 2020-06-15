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
} parameters;

parameters JetsonTx2 = {
    .cuda_max_readlen = 3.0,
    .cuda_avg_events_per_kmer = 2.0,
    .cuda_max_events_per_kmer = 5.0,
    .batch_size = 512,
    .batch_size_bases = 2350000,
    .num_thread = 6,
    .ultra_thresh = 100000
};

parameters Nanojet = {
    .cuda_max_readlen = 3.0,
    .cuda_avg_events_per_kmer = 2.0,
    .cuda_max_events_per_kmer = 5.0,
    .batch_size = 200,
    .batch_size_bases = 1400000,
    .num_thread = 4,
    .ultra_thresh = 100000
};

parameters Xavier = {
    .cuda_max_readlen = 3.0,
    .cuda_avg_events_per_kmer = 2.0,
    .cuda_max_events_per_kmer = 6.25,
    .batch_size = 1024,
    .batch_size_bases = 4700000,
    .num_thread = 8,
    .ultra_thresh = 100000
};

void set_opt_profile(opt_t *opt, parameters machine);
#endif
