#ifndef F5CMISC_CUH
#define F5CMISC_CUH

#include <stdint.h>

#define CUDA_DEBUG 1 //whether perform CUDA_device_synchronise or not

#define WARP_HACK 1 //whether the kernels are  performed in 1D with a warp hack (effective only  if specific TWODIM_ALIGN is not defined)

//align-core-kernel options

#define BLOCK_LEN_READS 1 //the block size along y axis (the number of reads)
#define BLOCK_LEN_BANDWIDTH 128 //the block size along the x axis, should be >= ALN_BANDWIDTH
#define ALIGN_KERNEL_FLOAT 1 //(for 2d kernel only)

//align-pre-kernel options
#define BLOCK_LEN_NUMBAND 16    //the block size along the x axis (BANDWDITH)
#define BLOCK_LEN_READS2 16 // //the block size along y axis (the number of reads)



/* check whether the last CUDA function or CUDA kernel launch is erroneous and if yes an error message will be printed
and then the program will be aborted*/
#define CUDA_CHK()                                                             \
    { gpu_assert(__FILE__, __LINE__); }



__global__ void 
//__launch_bounds__(MY_KERNEL_MAX_THREADS, MY_KERNEL_MIN_BLOCKS)
align_kernel_core_2d_shm(int32_t* read_len, int32_t* read_ptr,
    event_t* event_table, int32_t* n_events1,
    int32_t* event_ptr, 
    scalings_t* scalings, int32_t n_bam_rec,model_t* model_kmer_caches,float *band,uint8_t *traces, EventKmerPair* band_lower_lefts) ;

__global__ void align_kernel_pre_2d(char* read,
    int32_t* read_len, int32_t* read_ptr,
    int32_t* n_events,
    int32_t* event_ptr, model_t* models,
    int32_t n_bam_rec,model_t* model_kmer_caches,float *bands1,uint8_t *trace1, EventKmerPair* band_lower_left1) ;


__global__ void align_kernel_post(AlignedPair* event_align_pairs,
    int32_t* n_event_align_pairs, 
    int32_t* read_len, int32_t* read_ptr,
    event_t* event_table, int32_t* n_events,
    int32_t* event_ptr, 
    scalings_t* scalings, int32_t n_bam_rec,model_t* model_kmer_caches,float *bands1,uint8_t *trace1, EventKmerPair* band_lower_left1);




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
