#ifndef F5CMISC_CUH
#define F5CMISC_CUH

#include <stdint.h>

/* check whether the last CUDA function or CUDA kernel launch is erroneous and if yes an error message will be printed
and then the program will be aborted*/
#define CUDA_CHK()                                                             \
    { gpu_assert(__FILE__, __LINE__); }



__global__ void align_kernel(AlignedPair* event_align_pairs,
        int32_t* n_event_align_pairs, char* read,
        int32_t* read_len, int32_t* read_ptr,
        event_t* event_table, int32_t* n_events,
        int32_t* event_ptr, model_t* model,
        scalings_t* scalings, int32_t n_bam_rec,size_t* kmer_ranks,float *bands,uint8_t *trace, EventKmerPair* band_lower_left);






static inline void gpu_assert(const char* file, uint64_t line) {
    cudaError_t code = cudaGetLastError();
    if (code != cudaSuccess) {
        fprintf(stderr, "Cuda error: %s \n in file : %s line number : %lu\n",
                cudaGetErrorString(code), file, line);
        exit(-1);
    }
}

static inline int32_t cuda_exists() {
    //check cuda devices
    int32_t nDevices;
    cudaGetDeviceCount(&nDevices);
    if (nDevices == 0) {
        fprintf(stderr, "No CUDA device found. Use the CPU version\n");
        exit(1);
    }

    return nDevices;
}

static inline uint64_t cuda_freemem(int32_t devicenum) {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, devicenum);
    fprintf(stderr, "Device name: %s\n", prop.name);
    uint64_t golabalmem = prop.totalGlobalMem;
    fprintf(stderr, "Total global memory: %lf GB\n",
            (golabalmem / double(1024 * 1024 * 1024)));
    uint64_t freemem, total;
    cudaMemGetInfo(&freemem, &total);
    fprintf(stderr, "%lf GB free of total %lf GB\n",
            freemem / double(1024 * 1024 * 1024),
            total / double(1024 * 1024 * 1024));

    return freemem;
}


#endif