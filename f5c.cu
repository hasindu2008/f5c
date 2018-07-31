#include "f5c.cuh"

void* process_db_cuda(core_t* core, db_t* db) {
    event_table* et = (event_table*)malloc(sizeof(event_table) * db->n_bam_rec);
    MALLOC_CHK(et);

    int32_t i;
    for (i = 0; i < db->n_bam_rec; i++) {
        float* rawptr = db->f5[i]->rawptr;
        float range = db->f5[i]->range;
        float digitisation = db->f5[i]->digitisation;
        float offset = db->f5[i]->offset;
        int32_t nsample = db->f5[i]->nsample;

        // convert to pA
        float raw_unit = range / digitisation;
        for (int32_t j = 0; j < nsample; j++) {
            rawptr[j] = (rawptr[j] + offset) * raw_unit;
        }
        et[i] = getevents(db->f5[i]->nsample, rawptr);
    }

    return (void*)et;
}


// 	cudaMalloc((void**)&dev_sample, SAMPLES*KEYBYTES*sizeof(unsigned int));		checkCudaError();
// 	cudaMalloc((void**)&dev_corelation, KEYS*KEYBYTES*sizeof(double));			checkCudaError();
// 	cudaMalloc((void**)&dev_hammingArray, KEYS*KEYBYTES*SAMPLES*sizeof(byte));	checkCudaError();
// 	cudaMalloc((void**)&dev_wavestat, 2*WAVELENGTH*sizeof(double));				checkCudaError();
// 	cudaMalloc((void**)&dev_wavestat2, KEYS*KEYBYTES*WAVELENGTH*sizeof(double));checkCudaError();
// 	cudaMalloc((void**)&dev_hammingstat, 2*KEYS*KEYBYTES*sizeof(double));		checkCudaError();
	
// 	//make all correlation values 0 at the beginning
// 	cudaMemset(dev_corelation,0, KEYS*KEYBYTES*sizeof(double));		checkCudaError();
// 	//copy plain text samples to GPU
// 	cudaMemcpy(dev_sample,sample, SAMPLES*KEYBYTES*sizeof(unsigned int),cudaMemcpyHostToDevice);	checkCudaError();
	
// 	//cuda kernel configuraion parameters
// 	dim3 grid(KEYBYTES/16,KEYS/16);
// 	dim3 block(16,16);

// 	//find hamming statistics
// 	hammingkernel<<<grid,block>>>(dev_sample,dev_hammingArray,dev_hammingstat,SAMPLES);
// #ifdef DEBUG		
// 	cudaDeviceSynchronize(); checkCudaError();
// #endif