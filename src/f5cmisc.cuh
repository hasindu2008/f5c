#ifndef F5CMISC_CUH
#define F5CMISC_CUH

#include <stdint.h>

#define CUDA_DEBUG 1 //whether perform CUDA_device_synchronise or not
#define CUDA_PRE_MALLOC 1 //whether static cuda/cpu arrays are preallocated
#define CPU_GPU_PROC 1 //CUDA_PRE_MALLOC should be always 1 if this is set
//#define CUDA_DYNAMIC_MALLOC 1 //only effective with CPU_GPU_PROC (whether big dynamic loops are statically preallocated)
#define WARP_HACK 1 //whether the kernels are  performed in 1D with a warp hack (effective only  if specific TWODIM_ALIGN is not defined)

//align-core-kernel options

#define BLOCK_LEN_READS 1 //the block size along y axis (the number of reads)
#define BLOCK_LEN_BANDWIDTH 128 //the block size along the x axis, should be >= ALN_BANDWIDTH
#define ALIGN_KERNEL_FLOAT 1 //(for 2d kernel only)

//align-pre-kernel options
#define BLOCK_LEN_NUMBAND 16    //the block size along the x axis (BANDWDITH)
#define BLOCK_LEN_READS2 16 // //the block size along y axis (the number of reads)

#define AVG_EVENTS_PER_KMER 2.5f // the average number of events per base/k-mer
//AVG_EVENTS_PER_KMER is used to pre-allocate arrays on GPU that depends on the number of events

//if avverage events per base of a read < AVG_EVENTS_PER_KMER_GPU_THRESH process on GPU
//else go for the CPU
#define AVG_EVENTS_PER_KMER_GPU_THRESH 5.0f

#define AVG_EVENTS_PER_KMER_MAX 15.0f

#define TEGRA_MEM_FACTOR 0.8f //in tegra we cannot grab all 
//the free memory as we have to reserve some space for RAM as well
//TEGRA_MEM_FACTOR is the factor of the free memory allocated for the gpu

#define MEM_FACTOR 0.9f

#define REVERSAL_ON_CPU 1

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
    CUDA_CHK();
    if (nDevices == 0) {
        fprintf(stderr, "No CUDA device found. Run with --disable-cuda=yes to run on the CPU\n");
        exit(1);
    }

    return nDevices;
}

// static inline uint64_t cuda_freemem(int32_t devicenum) {
//     cudaDeviceProp prop;
//     cudaGetDeviceProperties(&prop, devicenum);
//     fprintf(stderr, "Device name: %s\n", prop.name);
//     uint64_t golabalmem = prop.totalGlobalMem;
//     fprintf(stderr, "Total global memory: %lf GB\n",
//             (golabalmem / double(1024 * 1024 * 1024)));
//     uint64_t freemem, total;
//     cudaMemGetInfo(&freemem, &total);
//     fprintf(stderr, "%lf GB free of total %lf GB\n",
//             freemem / double(1024 * 1024 * 1024),
//             total / double(1024 * 1024 * 1024));

//     return freemem;
// }

static inline uint64_t cuda_freemem(int32_t devicenum) {

    uint64_t freemem, total;
    cudaMemGetInfo(&freemem, &total); 
    CUDA_CHK();
    fprintf(stderr, "%lf GB free of total %lf GB global memory\n",
            freemem / double(1024 * 1024 * 1024),
            total / double(1024 * 1024 * 1024));

    return freemem;
}

#endif
