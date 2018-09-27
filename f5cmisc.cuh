#ifndef F5CMISC_CUH
#define F5CMISC_CUH

#include <stdint.h>

#define CUDA_DEBUG 1
//#define CONST_MEM 1
//#define DYNAMIC_PARALLELISM 1
#define DYNAMIC_THRESH 63
#define DYNAMIC_BLOCK_LEN 64

#define ALIGN_KERNEL_SLICED 1 //not to be used with CONST_MEM defined
#define WARP_HACK 1
#define TWODIM_KERNEL 1

#define MY_KERNEL_MAX_THREADS 128
#define MY_KERNEL_MIN_BLOCKS 64

//for 2d kernel
#define BLOCK_LEN_Y 1
#define BLOCK_LEN_X 128

/* check whether the last CUDA function or CUDA kernel launch is erroneous and if yes an error message will be printed
and then the program will be aborted*/
#define CUDA_CHK()                                                             \
    { gpu_assert(__FILE__, __LINE__); }


#ifndef CONST_MEM

    #ifndef ALIGN_KERNEL_SLICED 
        __global__ void align_kernel(AlignedPair* event_align_pairs,
            int32_t* n_event_align_pairs, char* read,
            int32_t* read_len, int32_t* read_ptr,
            event_t* event_table, int32_t* n_events,
            int32_t* event_ptr, model_t* model,
            scalings_t* scalings, int32_t n_bam_rec,int32_t* kmer_ranks,float *bands,uint8_t *trace, EventKmerPair* band_lower_left);
    #else
        __global__ void align_kernel_pre(AlignedPair* event_align_pairs,
        int32_t* n_event_align_pairs, char* read,
        int32_t* read_len, int32_t* read_ptr,
        event_t* event_table, int32_t* n_events,
        int32_t* event_ptr, model_t* model,
        scalings_t* scalings, int32_t n_bam_rec,int32_t* kmer_ranks,float *bands,uint8_t *trace, EventKmerPair* band_lower_left);
    
        __global__ void align_kernel_core(AlignedPair* event_align_pairs,
            int32_t* n_event_align_pairs, char* read,
            int32_t* read_len, int32_t* read_ptr,
            event_t* event_table, int32_t* n_events,
            int32_t* event_ptr, model_t* model,
            scalings_t* scalings, int32_t n_bam_rec,int32_t* kmer_ranks,float *bands,uint8_t *trace, EventKmerPair* band_lower_left);
    
        __global__ void align_kernel_post(AlignedPair* event_align_pairs,
            int32_t* n_event_align_pairs, char* read,
            int32_t* read_len, int32_t* read_ptr,
            event_t* event_table, int32_t* n_events,
            int32_t* event_ptr, model_t* model,
            scalings_t* scalings, int32_t n_bam_rec,int32_t* kmer_ranks,float *bands,uint8_t *trace, EventKmerPair* band_lower_left);
    #endif


    __global__ void 
    //__launch_bounds__(MY_KERNEL_MAX_THREADS, MY_KERNEL_MIN_BLOCKS)
    align_kernel_core_2d(AlignedPair* event_align_pairs,
        int32_t* n_event_align_pairs, char* read,
        int32_t* read_len, int32_t* read_ptr,
        event_t* event_table, int32_t* n_events1,
        int32_t* event_ptr, model_t* models,
        scalings_t* scalings, int32_t n_bam_rec,int32_t* kmer_rank,float *band,uint8_t *traces, EventKmerPair* band_lower_lefts) ;


#else
    __global__ void align_kernel(AlignedPair* event_align_pairs,
            int32_t* n_event_align_pairs, char* read,
            int32_t* read_len, int32_t* read_ptr,
            event_t* event_table, int32_t* n_events,
            int32_t* event_ptr,
            scalings_t* scalings, int32_t n_bam_rec,int32_t* kmer_ranks,float *bands,uint8_t *trace, EventKmerPair* band_lower_left);

#endif






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
