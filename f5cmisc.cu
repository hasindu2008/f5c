#include "f5c.cuh"
#include <stdint.h>

void gpu_assert(const char* file, uint64_t line) {
    cudaError_t code = cudaGetLastError();
    if (code != cudaSuccess) {
        fprintf(stderr, "Cuda error: %s \n in file : %s line number : %d\n",
                cudaGetErrorString(code), file, line);
        exit(-1);
    }
}

int32_t cuda_exists() {
    //check cuda devices
    int32_t nDevices;
    cudaGetDeviceCount(&nDevices);
    if (nDevices == 0) {
        fprintf(stderr, "No CUDA device found. Use the CPU version\n");
        exit(1);
    }

    return nDevices;
}

uint64_t cuda_freemem(int32_t devicenum) {
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
