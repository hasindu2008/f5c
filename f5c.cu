#include "f5c.cuh"

#define CUDA_DEBUG 1
#define BLOCK_LEN 128

__device__ int32_t align(AlignedPair* out_2, char* sequence, int32_t sequence_len,
    event_table events, model_t* models, scalings_t scaling,
    float sample_rate);

__global__ void align_kernel(){

	int i=blockDim.x*blockIdx.x+threadIdx.x;

	if (i < n_bam_rec){ 
        n_event_align_pairs[i]=align(event_align_pairs[i],read[i], read_len[i],et[i], model, scalings[i],
            sample_rate);
    }

}

void *align_cuda(core_t* core, db_t* db){
 
    int32_t n_bam_rec = db->n_bam_rec;

	//pointers for GPU
	AlignedPair* event_align_pairs;
    int32_t* n_event_align_pairs;

    char* read;
    int32_t* read_len;

    int32_t* n_events;
    event_t* event_table;

    model_t* model;
    scalings_t* scalings;

    //loop through and sum the total size for event_align_pairs and read and event_table
    int64_t sum_read_len;
    int64_t sum_n_events;

	//memory allocation in cuda
    //outputs
	cudaMalloc((void**)&event_align_pairs, 2*sum_n_events*sizeof(AlignedPair)); 	//todo : need better huristic	
    CUDA_CHK();
	cudaMalloc((void**)&n_event_align_pairs, n_bam_rec*sizeof(int32_t));		
    CUDA_CHK();

    //inputs
	//read sequences : need flattening
    
    cudaMalloc((void**)&read, sum_read_len*sizeof(char)+ n_bam_rec*sizeof(char));	//with null char
    CUDA_CHK();

    //read length : already linear
	cudaMalloc((void**)&read_len, n_bam_rec*sizeof(int32_t));	
    CUDA_CHK();
    cudaMemcpy(read_len,db->read_len,n_bam_rec*sizeof(int32_t),cudaMemcpyHostToDevice);	
    CUDA_CHK();

    //num events : need flattening
	cudaMalloc((void**)&n_events, n_bam_rec*sizeof(int32_t));				
    CUDA_CHK();
    //event table : need flattening
    cudaMemcpy(n_events,db->read_len, SAMPLES*KEYBYTES*n_bam_rec*sizeof(int32_t),cudaMemcpyHostToDevice);	
    CUDA_CHK();
	cudaMalloc((void**)&event_table, sum_n_events*sizeof(event_table));
    CUDA_CHK();

    //model : already linear
	cudaMalloc((void**)&model, NUM_KMER*sizeof(model_t));		//constant memory
    CUDA_CHK();


	cudaMalloc((void**)&scalings, n_bam_rec*sizeof(scalings_t));		
    CUDA_CHK();

	//make all correlation values 0 at the beginning
	cudaMemset(dev_corelation,0, KEYS*KEYBYTES*sizeof(double));		CUDA_CHK();
	//copy plain text samples to GPU
	cudaMemcpy(dev_sample,sample, SAMPLES*KEYBYTES*sizeof(unsigned int),cudaMemcpyHostToDevice);	CUDA_CHK();

	//cuda kernel configuraion parameters
	dim3 grid((db->n_bam_rec+BLOCK_LEN-1)/BLOCK_LEN);
    dim3 block(BLOCK_LEN);
        
    db->n_event_align_pairs[i] = 
        align_kernel<<<grid,block>>>(db->event_align_pairs[i],db->read[i], db->read_len[i],db->et[i], core->model, db->scalings[i],
                db->f5[i]->sample_rate);

    //fprintf(stderr,"readlen %d,n_events %d\n",db->read_len[i],n_event_align_pairs);    

    #ifdef CUDA_DEBUG		
        cudaDeviceSynchronize(); CUDA_CHK();
    #endif    

}


// void* process_db_cuda(core_t* core, db_t* db) {
//     event_table* et = (event_table*)malloc(sizeof(event_table) * db->n_bam_rec);
//     MALLOC_CHK(et);

//     int32_t i;
//     for (i = 0; i < db->n_bam_rec; i++) {
//         float* rawptr = db->f5[i]->rawptr;
//         float range = db->f5[i]->range;
//         float digitisation = db->f5[i]->digitisation;
//         float offset = db->f5[i]->offset;
//         int32_t nsample = db->f5[i]->nsample;

//         // convert to pA
//         float raw_unit = range / digitisation;
//         for (int32_t j = 0; j < nsample; j++) {
//             rawptr[j] = (rawptr[j] + offset) * raw_unit;
//         }
//         et[i] = getevents(db->f5[i]->nsample, rawptr);
//     }

//     return (void*)et;
// }

// 	cudaMalloc((void**)&dev_sample, SAMPLES*KEYBYTES*sizeof(unsigned int));		CUDA_CHK();
// 	cudaMalloc((void**)&dev_corelation, KEYS*KEYBYTES*sizeof(double));			CUDA_CHK();
// 	cudaMalloc((void**)&dev_hammingArray, KEYS*KEYBYTES*SAMPLES*sizeof(byte));	CUDA_CHK();
// 	cudaMalloc((void**)&dev_wavestat, 2*WAVELENGTH*sizeof(double));				CUDA_CHK();
// 	cudaMalloc((void**)&dev_wavestat2, KEYS*KEYBYTES*WAVELENGTH*sizeof(double));CUDA_CHK();
// 	cudaMalloc((void**)&dev_hammingstat, 2*KEYS*KEYBYTES*sizeof(double));		CUDA_CHK();

// 	//make all correlation values 0 at the beginning
// 	cudaMemset(dev_corelation,0, KEYS*KEYBYTES*sizeof(double));		CUDA_CHK();
// 	//copy plain text samples to GPU
// 	cudaMemcpy(dev_sample,sample, SAMPLES*KEYBYTES*sizeof(unsigned int),cudaMemcpyHostToDevice);	CUDA_CHK();

// 	//cuda kernel configuraion parameters
// 	dim3 grid(KEYBYTES/16,KEYS/16);
// 	dim3 block(16,16);

// 	//find hamming statistics
// 	hammingkernel<<<grid,block>>>(dev_sample,dev_hammingArray,dev_hammingstat,SAMPLES);
// #ifdef DEBUG
// 	cudaDeviceSynchronize(); CUDA_CHK();
// #endif