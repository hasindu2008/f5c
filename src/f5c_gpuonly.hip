/* @file f5c_gpuonly.cu
**
** implementation of the f5c GPU-only framework (opposed to the CPU-GPU hybrid approach in f5c.cu)
** not compiled by default
** @author: Hasindu Gamaarachchi (hasindu@unsw.edu.au)
** @@
******************************************************************************/

#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "error.h"
#include "f5c.h"
#include "f5cmisc.cuh"
#include "f5cmisc.h"


#ifndef CPU_GPU_PROC

/* if defined, static cuda/cpu arrays (model and arrays dependent on K)
   are preallocated at the beginning of the program, rather than repeatedly doing so inside a loop */
#define CUDA_PRE_MALLOC 1

void align_cuda(core_t* core, db_t* db) {
    int32_t i;
    int32_t n_bam_rec = db->n_bam_rec;
    double realtime1;

    /**cuda pointers*/
    char* read;        //flattened reads sequences
    ptr_t* read_ptr; //index pointer for flattedned "reads"
    int32_t* read_len;
    int64_t sum_read_len;
    int32_t* n_events;
    event_t* event_table;
    ptr_t* event_ptr;
    int64_t sum_n_events;
    scalings_t* scalings;
    AlignedPair* event_align_pairs;
    int32_t* n_event_align_pairs;
    float *bands;
    uint8_t *trace;
    EventKmerPair* band_lower_left;

realtime1 = realtime();

    int32_t cuda_device_num = core->opt.cuda_dev_id;
    cudaSetDevice(cuda_device_num);
    CUDA_CHK();

#ifdef CUDA_PRE_MALLOC
    ptr_t* read_ptr_host = core->cuda->read_ptr_host;
#else
    //get the total size and create the pointers
    ptr_t* read_ptr_host = (ptr_t*)malloc(sizeof(ptr_t) * n_bam_rec);
    MALLOC_CHK(read_ptr_host);
#endif
    sum_read_len = 0;

    //read sequences : needflattening
    for (i = 0; i < n_bam_rec; i++) {
        read_ptr_host[i] = sum_read_len;
        sum_read_len += (db->read_len[i] + 1); //with null term
    }
    //form the temporary flattened array on host
    char* read_host = (char*)malloc(sizeof(char) * sum_read_len);
    MALLOC_CHK(read_host);
    for (i = 0; i < n_bam_rec; i++) {
        ptr_t idx = read_ptr_host[i];
        strcpy(&read_host[idx], db->read[i]);
    }

    //now the events : need flattening
    //num events : need flattening
    //get the total size and create the pointers
#ifdef CUDA_PRE_MALLOC
    int32_t* n_events_host = core->cuda->n_events_host;
    ptr_t* event_ptr_host = core->cuda->event_ptr_host;
#else
    int32_t* n_events_host = (int32_t*)malloc(sizeof(int32_t) * n_bam_rec);
    MALLOC_CHK(n_events_host);
    ptr_t* event_ptr_host = (ptr_t*)malloc(sizeof(ptr_t) * n_bam_rec);
    MALLOC_CHK(event_ptr_host);
#endif

    sum_n_events = 0;
    for (i = 0; i < n_bam_rec; i++) {
        n_events_host[i] = db->et[i].n;
        event_ptr_host[i] = sum_n_events;
        sum_n_events += db->et[i].n;
    }

    //event table flatten
    //form the temporary flattened array on host
    event_t* event_table_host =
        (event_t*)malloc(sizeof(event_t) * sum_n_events);
    MALLOC_CHK(event_table_host);
    for (i = 0; i < n_bam_rec; i++) {
        ptr_t idx = event_ptr_host[i];
        memcpy(&event_table_host[idx], db->et[i].event,
               sizeof(event_t) * db->et[i].n);
    }

    AlignedPair* event_align_pairs_host =
        (AlignedPair*)malloc(2 * sum_n_events * sizeof(AlignedPair));
    MALLOC_CHK(event_align_pairs_host);

core->align_cuda_preprocess += (realtime() - realtime1);

    /** Start GPU mallocs**/
realtime1 = realtime();

#ifdef CUDA_PRE_MALLOC
    read_ptr =core->cuda->read_ptr;
    read_len=core->cuda->read_len;
    n_events=core->cuda->n_events;
    event_ptr=core->cuda->event_ptr;
    scalings=core->cuda->scalings;
    model_t* model = core->cuda->model;
#else

    if(core->opt.verbosity>1) print_size("read_ptr array",n_bam_rec * sizeof(ptr_t));
    cudaMalloc((void**)&read_ptr, n_bam_rec * sizeof(ptr_t));
    CUDA_CHK();

    if(core->opt.verbosity>1) print_size("read_lens",n_bam_rec * sizeof(int32_t));
    cudaMalloc((void**)&read_len, n_bam_rec * sizeof(int32_t));
    CUDA_CHK();
    //n_events
    if(core->opt.verbosity>1) print_size("n_events",n_bam_rec * sizeof(int32_t));
    cudaMalloc((void**)&n_events, n_bam_rec * sizeof(int32_t));
    CUDA_CHK();
    //event ptr
    if(core->opt.verbosity>1) print_size("event ptr",n_bam_rec * sizeof(ptr_t));
    cudaMalloc((void**)&event_ptr, n_bam_rec * sizeof(ptr_t));
    CUDA_CHK();
    //scalings : already linear
    if(core->opt.verbosity>1) print_size("Scalings",n_bam_rec * sizeof(scalings_t));
    cudaMalloc((void**)&scalings, n_bam_rec * sizeof(scalings_t));
    CUDA_CHK();
    //model : already linear
    model_t* model;
    cudaMalloc((void**)&model,
            MAX_NUM_KMER * sizeof(model_t));
    CUDA_CHK();
#endif


    if(core->opt.verbosity>1) print_size("read array",sum_read_len * sizeof(char));
    cudaMalloc((void**)&read, sum_read_len * sizeof(char)); //with null char
    CUDA_CHK();
    if(core->opt.verbosity>1) print_size("event table",sum_n_events * sizeof(event_t));
    cudaMalloc((void**)&event_table, sum_n_events * sizeof(event_t));
    CUDA_CHK();
    model_t* model_kmer_cache;
    cudaMalloc((void**)&model_kmer_cache, sum_read_len * sizeof(model_t));
    CUDA_CHK();

    /**allocate output arrays for cuda**/
    if(core->opt.verbosity>1) print_size("event align pairs",2 * sum_n_events *sizeof(AlignedPair));
    cudaMalloc((void**)&event_align_pairs,
            2 * sum_n_events *
                sizeof(AlignedPair)); //todo : need better huristic
    CUDA_CHK();
#ifdef CUDA_PRE_MALLOC
    n_event_align_pairs=core->cuda->n_event_align_pairs;
#else
    if(core->opt.verbosity>1) print_size("n_event_align_pairs",n_bam_rec * sizeof(int32_t));
    cudaMalloc((void**)&n_event_align_pairs, n_bam_rec * sizeof(int32_t));
    CUDA_CHK();
#endif
    //scratch arrays
    size_t sum_n_bands = sum_n_events + sum_read_len; //todo : can be optimised
    if(core->opt.verbosity>1) print_size("bands",sizeof(float) * sum_n_bands * ALN_BANDWIDTH);
    cudaMalloc((void**)&bands,sizeof(float) * sum_n_bands * ALN_BANDWIDTH);
    CUDA_CHK();
    if(core->opt.verbosity>1) print_size("trace",sizeof(uint8_t) * sum_n_bands * ALN_BANDWIDTH);
    cudaMalloc((void**)&trace, sizeof(uint8_t) * sum_n_bands * ALN_BANDWIDTH);
    CUDA_CHK();
    cudaMemset(trace,0,sizeof(uint8_t) * sum_n_bands * ALN_BANDWIDTH); //initialise the trace array to 0
    if(core->opt.verbosity>1) print_size("band_lower_left",sizeof(EventKmerPair)* sum_n_bands);
    cudaMalloc((void**)&band_lower_left, sizeof(EventKmerPair)* sum_n_bands);
    CUDA_CHK();
core->align_cuda_malloc += (realtime() - realtime1);

    /* cuda mem copys*/
realtime1 =realtime();
    cudaMemcpy(read_ptr, read_ptr_host, n_bam_rec * sizeof(ptr_t),
               cudaMemcpyHostToDevice);
    CUDA_CHK();
    cudaMemcpy(read, read_host, sum_read_len * sizeof(char),
               cudaMemcpyHostToDevice);
    CUDA_CHK();
    //read length : already linear hence direct copy
    cudaMemcpy(read_len, db->read_len, n_bam_rec * sizeof(int32_t),
               cudaMemcpyHostToDevice);
    CUDA_CHK();
    cudaMemcpy(n_events, n_events_host, n_bam_rec * sizeof(int32_t),
               cudaMemcpyHostToDevice);
    CUDA_CHK();
    cudaMemcpy(event_ptr, event_ptr_host, n_bam_rec * sizeof(ptr_t),
               cudaMemcpyHostToDevice);
    CUDA_CHK();
    cudaMemcpy(event_table, event_table_host, sizeof(event_t) * sum_n_events,
               cudaMemcpyHostToDevice);
    CUDA_CHK();

#ifndef CUDA_PRE_MALLOC
//model : already linear //move to cuda_init
    cudaMemcpy(model, core->model, MAX_NUM_KMER * sizeof(model_t),
            cudaMemcpyHostToDevice);
    CUDA_CHK();
#endif
    //can be interleaved
    cudaMemcpy(scalings, db->scalings, sizeof(scalings_t) * n_bam_rec,
               cudaMemcpyHostToDevice);
    CUDA_CHK();
core->align_cuda_memcpy += (realtime() - realtime1);

    uint32_t kmer_size = core->kmer_size;

realtime1 = realtime();

    /*pre kernel*/
    assert(BLOCK_LEN_BANDWIDTH>=ALN_BANDWIDTH);
    dim3 gridpre(1,(db->n_bam_rec + BLOCK_LEN_READS - 1) / BLOCK_LEN_READS);
    dim3 blockpre(BLOCK_LEN_BANDWIDTH,BLOCK_LEN_READS);
	if(core->opt.verbosity>1) fprintf(stderr,"grid %d,%d, block %d,%d\n",gridpre.x,gridpre.y, blockpre.x,blockpre.y);

    align_kernel_pre_2d<<<gridpre, blockpre>>>( read,
        read_len, read_ptr, n_events,
        event_ptr, model, kmer_size, n_bam_rec, model_kmer_cache,bands,trace,band_lower_left);

    cudaDeviceSynchronize();CUDA_CHK();
    if(core->opt.verbosity>1) fprintf(stderr, "[%s::%.3f*%.2f] align-pre kernel done\n", __func__,
            realtime() - realtime1, cputime() / (realtime() - realtime1));
core->align_kernel_time += (realtime() - realtime1);
core->align_pre_kernel_time += (realtime() - realtime1);

realtime1 = realtime();

    /* core kernel*/
    assert(BLOCK_LEN_BANDWIDTH>=ALN_BANDWIDTH);
    dim3 grid1(1,(db->n_bam_rec + BLOCK_LEN_READS - 1) / BLOCK_LEN_READS);
    dim3 block1(BLOCK_LEN_BANDWIDTH,BLOCK_LEN_READS);
    align_kernel_core_2d_shm<<<grid1, block1>>>(read_len, read_ptr, event_table, n_events,
            event_ptr, scalings, n_bam_rec, model_kmer_cache, kmer_size, bands,trace,band_lower_left );

    cudaDeviceSynchronize();CUDA_CHK();
    if(core->opt.verbosity>1) fprintf(stderr, "[%s::%.3f*%.2f] align-core kernel done\n", __func__,
    realtime() - realtime1, cputime() / (realtime() - realtime1));
    core->align_kernel_time += (realtime() - realtime1);
core->align_core_kernel_time += (realtime() - realtime1);

realtime1 = realtime();

    /*post kernel*/
    int32_t BLOCK_LEN = core->opt.cuda_block_size;
    dim3 gridpost((db->n_bam_rec + BLOCK_LEN - 1) / BLOCK_LEN);
    dim3 blockpost(BLOCK_LEN);
    #ifndef WARP_HACK
        align_kernel_post<<<gridpost, blockpost>>>(event_align_pairs, n_event_align_pairs,
            read_len, read_ptr, event_table, n_events,
            event_ptr,scalings, n_bam_rec, model_kmer_cache, kmer_size, bands,trace,band_lower_left );

    #else
        assert(BLOCK_LEN>=32);
        dim3 grid1post((db->n_bam_rec + (BLOCK_LEN/32) - 1) / (BLOCK_LEN/32));
        if(core->opt.verbosity>1) fprintf(stderr,"grid new %d\n",grid1post.x);
        align_kernel_post<<<grid1post, blockpost>>>(event_align_pairs, n_event_align_pairs,
            read_len, read_ptr, event_table, n_events,
            event_ptr, scalings, n_bam_rec, model_kmer_cache, kmer_size, bands,trace,band_lower_left );
    #endif
    cudaDeviceSynchronize();CUDA_CHK();
    if(core->opt.verbosity>1) fprintf(stderr, "[%s::%.3f*%.2f] align-post kernel done\n", __func__,
            realtime() - realtime1, cputime() / (realtime() - realtime1));
    core->align_kernel_time += (realtime() - realtime1);
core->align_post_kernel_time += (realtime() - realtime1);


    //fprintf(stderr,"readlen %d,n_events %d\n",db->read_len[i],n_event_align_pairs);

#ifdef CUDA_DEBUG

    cudaDeviceSynchronize();
    CUDA_CHK();

#endif

    /** copyback ans**/
realtime1 =  realtime();
    cudaMemcpy(db->n_event_align_pairs, n_event_align_pairs,
               n_bam_rec * sizeof(int32_t), cudaMemcpyDeviceToHost);
    CUDA_CHK();

    cudaMemcpy(event_align_pairs_host, event_align_pairs,
               2 * sum_n_events * sizeof(AlignedPair), cudaMemcpyDeviceToHost);
    CUDA_CHK();
core->align_cuda_memcpy += (realtime() - realtime1);

realtime1 =  realtime();
#ifndef CUDA_PRE_MALLOC
    cudaFree(read_ptr);
    cudaFree(read_len);
    cudaFree(n_events);
    cudaFree(event_ptr);
    cudaFree(model); //constant memory
    cudaFree(scalings);
    cudaFree(n_event_align_pairs);
#endif
    cudaFree(read); //with null char
    cudaFree(event_table);
    cudaFree(event_align_pairs);
    cudaFree(bands);
    cudaFree(trace);
    cudaFree(band_lower_left);
    cudaFree(model_kmer_cache);

core->align_cuda_malloc += (realtime() - realtime1);

    /** post work**/
realtime1 =  realtime();
    //copy back
    for (i = 0; i < n_bam_rec; i++) {
        ptr_t idx = event_ptr_host[i];
        memcpy(db->event_align_pairs[i], &event_align_pairs_host[idx * 2],
               sizeof(AlignedPair) * db->n_event_align_pairs[i]);
    }

    //free the temp arrays on host
#ifndef CUDA_PRE_MALLOC
    free(read_ptr_host);
    free(n_events_host);
    free(event_ptr_host);
#endif
    free(read_host);
    free(event_table_host);
    free(event_align_pairs_host);


core->align_cuda_postprocess += (realtime() - realtime1);

}


#endif
