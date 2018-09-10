#ifndef F5CMISC_CUH
#define F5CMISC_CUH

#include <stdint.h>

/* check whether the last CUDA function or CUDA kernel launch is erroneous and if yes an error message will be printed
and then the program will be aborted*/

#define CUDA_CHK()                                                             \
    { gpu_assert(__FILE__, __LINE__); }

void gpu_assert(const char* file, uint64_t line);

#endif