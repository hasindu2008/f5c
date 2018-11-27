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

void init_cuda(core_t* core){

    core->cuda = (cuda_data_t*)malloc(sizeof(cuda_data_t)); 
    MALLOC_CHK(core->cuda);

    core->align_kernel_time=0;
    core->align_pre_kernel_time=0;
    core->align_core_kernel_time=0;
    core->align_post_kernel_time=0;
    core->align_cuda_malloc=0;
    core->extra_load_cpu=0;
    core->align_cuda_memcpy=0;
    core->align_cuda_postprocess=0;
    core->align_cuda_preprocess=0;

#ifdef CUDA_PRE_MALLOC

    int32_t n_bam_rec = core->opt.batch_size;
    //cpu arrays
    core->cuda->read_ptr_host = (int32_t*)malloc(sizeof(int32_t) * n_bam_rec);
    MALLOC_CHK(core->cuda->read_ptr_host);
    core->cuda->n_events_host = (int32_t*)malloc(sizeof(int32_t) * n_bam_rec);
    MALLOC_CHK(core->cuda->n_events_host);
    core->cuda->event_ptr_host = (int32_t*)malloc(sizeof(int32_t) * n_bam_rec);
    MALLOC_CHK(core->cuda->event_ptr_host);

    core->cuda->read_len_host = (int32_t*)malloc(sizeof(int32_t) * n_bam_rec);
    MALLOC_CHK(core->cuda->read_len_host);
    core->cuda->scalings_host = (scalings_t*)malloc(sizeof(scalings_t) * n_bam_rec);
    MALLOC_CHK(core->cuda->scalings_host);
    core->cuda->n_event_align_pairs_host = (int32_t*)malloc(sizeof(int32_t) * n_bam_rec);
    MALLOC_CHK(core->cuda->n_event_align_pairs_host);

    //cuda arrays
    print_size("read_ptr array",n_bam_rec * sizeof(int32_t));
    cudaMalloc((void**)&(core->cuda->read_ptr), n_bam_rec * sizeof(int32_t));
    CUDA_CHK();
    print_size("read_lens",n_bam_rec * sizeof(int32_t));
    cudaMalloc((void**)&(core->cuda->read_len), n_bam_rec * sizeof(int32_t));
    CUDA_CHK();
    //n_events
    print_size("n_events",n_bam_rec * sizeof(int32_t));
    cudaMalloc((void**)&(core->cuda->n_events), n_bam_rec * sizeof(int32_t));
    CUDA_CHK();
    //event ptr
    print_size("event ptr",n_bam_rec * sizeof(int32_t));
    cudaMalloc((void**)&(core->cuda->event_ptr), n_bam_rec * sizeof(int32_t));
    CUDA_CHK();
    //scalings : already linear
    print_size("Scalings",n_bam_rec * sizeof(scalings_t));
    cudaMalloc((void**)&(core->cuda->scalings), n_bam_rec * sizeof(scalings_t));
    CUDA_CHK();
    cudaMalloc((void**)&(core->cuda->model),
            NUM_KMER * sizeof(model_t));
    CUDA_CHK();  

    print_size("n_event_align_pairs",n_bam_rec * sizeof(int32_t));
    cudaMalloc((void**)&(core->cuda->n_event_align_pairs), n_bam_rec * sizeof(int32_t));
    CUDA_CHK();

    //model : already linear //move to cuda_init
    cudaMemcpy(core->cuda->model, core->model, NUM_KMER * sizeof(model_t),
    cudaMemcpyHostToDevice);
    CUDA_CHK();

    

#endif

    return;
}

void free_cuda(core_t* core){

#ifdef CUDA_PRE_MALLOC
    free(core->cuda->event_ptr_host);
    free(core->cuda->n_events_host);
    free(core->cuda->read_ptr_host);
    free(core->cuda->read_len_host);
    free(core->cuda->scalings_host);
    free(core->cuda->n_event_align_pairs_host);
    
    cudaFree(core->cuda->read_ptr);
    cudaFree(core->cuda->read_len);
    cudaFree(core->cuda->n_events);
    cudaFree(core->cuda->event_ptr);
    cudaFree(core->cuda->model); //constant memory
    cudaFree(core->cuda->scalings);
    cudaFree(core->cuda->n_event_align_pairs);

#endif

    free(core->cuda);
    return;
}


#ifndef CPU_GPU_PROC

void align_cuda(core_t* core, db_t* db) {
    int32_t i;
    int32_t n_bam_rec = db->n_bam_rec;
    double realtime1;

    /**cuda pointers*/
    char* read;        //flattened reads sequences
    int32_t* read_ptr; //index pointer for flattedned "reads"
    int32_t* read_len;
    int64_t sum_read_len;
    int32_t* n_events;
    event_t* event_table;
    int32_t* event_ptr;
    int64_t sum_n_events;
    scalings_t* scalings;
    AlignedPair* event_align_pairs;
    int32_t* n_event_align_pairs;
    float *bands;
    uint8_t *trace;
    EventKmerPair* band_lower_left;

realtime1 = realtime();

#ifdef CUDA_PRE_MALLOC
    int32_t* read_ptr_host = core->cuda->read_ptr_host;
#else
    //get the total size and create the pointers
    int32_t* read_ptr_host = (int32_t*)malloc(sizeof(int32_t) * n_bam_rec);
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
        int32_t idx = read_ptr_host[i];
        strcpy(&read_host[idx], db->read[i]);
    }

    //now the events : need flattening
    //num events : need flattening
    //get the total size and create the pointers
#ifdef CUDA_PRE_MALLOC
    int32_t* n_events_host = core->cuda->n_events_host;
    int32_t* event_ptr_host = core->cuda->event_ptr_host;
#else
    int32_t* n_events_host = (int32_t*)malloc(sizeof(int32_t) * n_bam_rec);
    MALLOC_CHK(n_events_host);
    int32_t* event_ptr_host = (int32_t*)malloc(sizeof(int32_t) * n_bam_rec);
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
        int32_t idx = event_ptr_host[i];
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

    print_size("read_ptr array",n_bam_rec * sizeof(int32_t));
    cudaMalloc((void**)&read_ptr, n_bam_rec * sizeof(int32_t));
    CUDA_CHK();

    print_size("read_lens",n_bam_rec * sizeof(int32_t));
    cudaMalloc((void**)&read_len, n_bam_rec * sizeof(int32_t));
    CUDA_CHK();
    //n_events
    print_size("n_events",n_bam_rec * sizeof(int32_t));
    cudaMalloc((void**)&n_events, n_bam_rec * sizeof(int32_t));
    CUDA_CHK();
    //event ptr
    print_size("event ptr",n_bam_rec * sizeof(int32_t));
    cudaMalloc((void**)&event_ptr, n_bam_rec * sizeof(int32_t));
    CUDA_CHK();
    //scalings : already linear
    print_size("Scalings",n_bam_rec * sizeof(scalings_t));
    cudaMalloc((void**)&scalings, n_bam_rec * sizeof(scalings_t));
    CUDA_CHK();
    //model : already linear
    model_t* model;
    cudaMalloc((void**)&model,
            NUM_KMER * sizeof(model_t));
    CUDA_CHK();  
#endif


    print_size("read array",sum_read_len * sizeof(char));
    cudaMalloc((void**)&read, sum_read_len * sizeof(char)); //with null char
    CUDA_CHK();
    print_size("event table",sum_n_events * sizeof(event_t));
    cudaMalloc((void**)&event_table, sum_n_events * sizeof(event_t));
    CUDA_CHK();
    model_t* model_kmer_cache;
    cudaMalloc((void**)&model_kmer_cache, sum_read_len * sizeof(model_t)); 
    CUDA_CHK();
 
    /**allocate output arrays for cuda**/
    print_size("event align pairs",2 * sum_n_events *sizeof(AlignedPair));
    cudaMalloc((void**)&event_align_pairs,
            2 * sum_n_events *
                sizeof(AlignedPair)); //todo : need better huristic
    CUDA_CHK();
#ifdef CUDA_PRE_MALLOC
    n_event_align_pairs=core->cuda->n_event_align_pairs;
#else
    print_size("n_event_align_pairs",n_bam_rec * sizeof(int32_t));
    cudaMalloc((void**)&n_event_align_pairs, n_bam_rec * sizeof(int32_t));
    CUDA_CHK();
#endif
    //scratch arrays
    size_t sum_n_bands = sum_n_events + sum_read_len; //todo : can be optimised 
    print_size("bands",sizeof(float) * sum_n_bands * ALN_BANDWIDTH);
    cudaMalloc((void**)&bands,sizeof(float) * sum_n_bands * ALN_BANDWIDTH);
    CUDA_CHK();
    print_size("trace",sizeof(uint8_t) * sum_n_bands * ALN_BANDWIDTH);
    cudaMalloc((void**)&trace, sizeof(uint8_t) * sum_n_bands * ALN_BANDWIDTH);
    CUDA_CHK();
    cudaMemset(trace,0,sizeof(uint8_t) * sum_n_bands * ALN_BANDWIDTH); //initialise the trace array to 0
    print_size("band_lower_left",sizeof(EventKmerPair)* sum_n_bands);
    cudaMalloc((void**)&band_lower_left, sizeof(EventKmerPair)* sum_n_bands);
    CUDA_CHK();   
core->align_cuda_malloc += (realtime() - realtime1);

    /* cuda mem copys*/
realtime1 =realtime();
    cudaMemcpy(read_ptr, read_ptr_host, n_bam_rec * sizeof(int32_t),
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
    cudaMemcpy(event_ptr, event_ptr_host, n_bam_rec * sizeof(int32_t),
               cudaMemcpyHostToDevice);
    CUDA_CHK();
    cudaMemcpy(event_table, event_table_host, sizeof(event_t) * sum_n_events,
               cudaMemcpyHostToDevice);
    CUDA_CHK();

#ifndef CUDA_PRE_MALLOC
//model : already linear //move to cuda_init
    cudaMemcpy(model, core->model, NUM_KMER * sizeof(model_t),
            cudaMemcpyHostToDevice);
    CUDA_CHK();
#endif
    //can be interleaved
    cudaMemcpy(scalings, db->scalings, sizeof(scalings_t) * n_bam_rec,
               cudaMemcpyHostToDevice);
    CUDA_CHK();
core->align_cuda_memcpy += (realtime() - realtime1);



realtime1 = realtime();    
 
    /*pre kernel*/
    assert(BLOCK_LEN_BANDWIDTH>=ALN_BANDWIDTH);
    dim3 gridpre(1,(db->n_bam_rec + BLOCK_LEN_READS - 1) / BLOCK_LEN_READS);
    dim3 blockpre(BLOCK_LEN_BANDWIDTH,BLOCK_LEN_READS);  
	fprintf(stderr,"grid %d,%d, block %d,%d\n",gridpre.x,gridpre.y, blockpre.x,blockpre.y);	

    align_kernel_pre_2d<<<gridpre, blockpre>>>( read,
        read_len, read_ptr, n_events,
        event_ptr, model, n_bam_rec, model_kmer_cache,bands,trace,band_lower_left); 
       
    cudaDeviceSynchronize();CUDA_CHK();
    fprintf(stderr, "[%s::%.3f*%.2f] align pre done\n", __func__,
            realtime() - realtime1, cputime() / (realtime() - realtime1));
core->align_kernel_time += (realtime() - realtime1);        
core->align_pre_kernel_time += (realtime() - realtime1);        
                
realtime1 = realtime();

    /* core kernel*/
    assert(BLOCK_LEN_BANDWIDTH>=ALN_BANDWIDTH);
    dim3 grid1(1,(db->n_bam_rec + BLOCK_LEN_READS - 1) / BLOCK_LEN_READS);
    dim3 block1(BLOCK_LEN_BANDWIDTH,BLOCK_LEN_READS);
    align_kernel_core_2d_shm<<<grid1, block1>>>(read_len, read_ptr, event_table, n_events,
            event_ptr, scalings, n_bam_rec, model_kmer_cache,bands,trace,band_lower_left );
    
    cudaDeviceSynchronize();CUDA_CHK();
    fprintf(stderr, "[%s::%.3f*%.2f] align done\n", __func__,
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
            event_ptr,scalings, n_bam_rec, model_kmer_cache,bands,trace,band_lower_left );     

    #else
        assert(BLOCK_LEN>=32);    
        dim3 grid1post((db->n_bam_rec + (BLOCK_LEN/32) - 1) / (BLOCK_LEN/32)); 
        fprintf(stderr,"grid new %d\n",grid1post.x);   
        align_kernel_post<<<grid1post, blockpost>>>(event_align_pairs, n_event_align_pairs,
            read_len, read_ptr, event_table, n_events,
            event_ptr, scalings, n_bam_rec, model_kmer_cache,bands,trace,band_lower_left );
    #endif
    cudaDeviceSynchronize();CUDA_CHK();
    fprintf(stderr, "[%s::%.3f*%.2f] align post done\n", __func__,
            realtime() - realtime1, cputime() / (realtime() - realtime1));
    core->align_kernel_time += (realtime() - realtime1);        
core->align_post_kernel_time += (realtime() - realtime1);        


    //fprintf(stderr,"readlen %d,n_events %d\n",db->read_len[i],n_event_align_pairs);

#ifdef CUDA_DEBUG

    cudaDeviceSynchronize();
    cudaError_t code = cudaGetLastError();
    //todo : print a message to detect the launch timed out
    if (code == cudaErrorLaunchTimeout) {
        ERROR("%s", "The kernel timed out. You have to first disable the cuda "
                    "time out.");
        fprintf(
            stderr,
            "On Ubuntu do the following\nOpen the file /etc/X11/xorg.conf\nYou "
            "will have a section about your NVIDIA device. Add the following "
            "line to it.\nOption \"Interactive\" \"0\"\nIf you do not have a "
            "section about your NVIDIA device in /etc/X11/xorg.conf or you do "
            "not have a file named /etc/X11/xorg.conf, run the command sudo "
            "nvidia-xconfig to generate a xorg.conf file and do as above.\n\n");
    }
    if (code != cudaSuccess) {
        fprintf(stderr, "Cuda error: %s \n in file : %s line number : %lu\n",
                cudaGetErrorString(code), __FILE__, __LINE__);
        exit(-1);
    }        
    
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
        int32_t idx = event_ptr_host[i];
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



#else



#ifdef WORK_STEAL
static inline int32_t steal_work(pthread_arg_t* all_args, int32_t n_threads)
{

	int32_t i, c_i = -1;
	int32_t k;
	for (i = 0; i < n_threads; ++i){
        pthread_arg_t args = all_args[i];
        //fprintf(stderr,"endi : %d, starti : %d\n",args.endi,args.starti);
		if (args.endi-args.starti > STEAL_THRESH) {
            //fprintf(stderr,"gap : %d\n",args.endi-args.starti);
            c_i = i;
            break;
        }
    }
    if(c_i<0){
        return -1;
    }
	k = __sync_fetch_and_add(&(all_args[c_i].starti), 1);
    //fprintf(stderr,"k : %d, end %d, start %d\n",k,all_args[c_i].endi,all_args[c_i].starti);
	return k >= all_args[c_i].endi ? -1 : k;
}
#endif

void* pthread_cusingle(void* voidargs) {
    int32_t i,j;
    pthread_arg_t* args = (pthread_arg_t*)voidargs;
    db_t* db = args->db;
    core_t* core = args->core;

#ifndef WORK_STEAL
    for (i = args->starti; i < args->endi; i++) {
        j=args->ultra_long_reads[i];
        args->func(core,db,j);
    }
#else
    pthread_arg_t* all_args = (pthread_arg_t*)(args->all_pthread_args);
    //adapted from ktherad
	for (;;) {
		i = __sync_fetch_and_add(&args->starti, 1);
		if (i >= args->endi) {
            break;
        }
        j=args->ultra_long_reads[i];
		args->func(core,db,j);
	}
	while ((i = steal_work(all_args,core->opt.num_thread)) >= 0){
        j=args->ultra_long_reads[i];
		args->func(core,db,j);  
    }  
#endif

    //fprintf(stderr,"Thread %d done\n",(myargs->position)/THREADS);
    pthread_exit(0);
}


void pthread_cudb(core_t* core, db_t* db, int32_t* ultra_long_reads, int32_t  n_ultra_long_reads,void (*func)(core_t*,db_t*,int)){ 
    //create threads
    pthread_t tids[core->opt.num_thread];
    pthread_arg_t pt_args[core->opt.num_thread];
    int32_t t, ret;
    int32_t i = 0;
    int32_t num_thread = core->opt.num_thread;
    int32_t step = (n_ultra_long_reads + num_thread - 1) / num_thread;
    //todo : check for higher num of threads than the data
    //current works but many threads are created despite 

    //set the data structures
    for (t = 0; t < num_thread; t++) {
        pt_args[t].core = core;
        pt_args[t].db = db;
        pt_args[t].starti = i;
        i += step;
        if (i > n_ultra_long_reads) {
            pt_args[t].endi = n_ultra_long_reads;
        } else {
            pt_args[t].endi = i;
        }
        pt_args[t].func=func;
        pt_args[t].ultra_long_reads=ultra_long_reads;
    #ifdef WORK_STEAL    
        pt_args[t].all_pthread_args =  (void *)pt_args;
    #endif
        //fprintf(stderr,"t%d : %d-%d\n",t,pt_args[t].starti,pt_args[t].endi);

    }

    //create threads
    for(t = 0; t < core->opt.num_thread; t++){
        ret = pthread_create(&tids[t], NULL, pthread_cusingle,
                                (void*)(&pt_args[t]));
        NEG_CHK(ret);
    }

    //pthread joining
    for (t = 0; t < core->opt.num_thread; t++) {
        int ret = pthread_join(tids[t], NULL);
        NEG_CHK(ret);
    }
}


void* align_cudb(void* voidargs){

    pthread_arg_t* args = (pthread_arg_t*)voidargs;
    db_t* db = args->db;
    core_t* core = args->core;
    int32_t* ultra_long_reads = args->ultra_long_reads;
    int32_t n_ultra_long_reads = args->endi;
    //fprintf(stderr,"ultra long guys : %d\n",n_ultra_long_reads);
    //fprintf(stderr, "cpu\n");
    if (core->opt.num_thread == 1) {
        int j;
        for(j=0;j<n_ultra_long_reads;j++) {
            int32_t i = ultra_long_reads[j];
            align_single(core, db, i);
            // db->n_event_align_pairs[i] =
            //     align(db->event_align_pairs[i], db->read[i],
            //           db->read_len[i], db->et[i], core->model,
            //           db->scalings[i], db->f5[i]->sample_rate);
            //fprintf(stderr,"readlen %d,n_events %d\n",db->read_len[i],n_event_align_pairs);
        }
    } else {
        pthread_cudb(core, db, ultra_long_reads,n_ultra_long_reads,align_single);
        
    }

    fprintf(stderr,"%d reads (length>%d kb) processed on cpu\n",n_ultra_long_reads,GPU_MAX_READ_LEN/1000);

    return NULL;
}
    
pthread_t align_cudb_async(pthread_arg_t *pt_args,core_t* core, db_t* db, int32_t* ultra_long_reads, int32_t  n_ultra_long_reads) {
    assert(pt_args==NULL);
    pt_args = (pthread_arg_t *)malloc(sizeof(pthread_arg_t));
    MALLOC_CHK(pt_args);
    pt_args->core = core;
    pt_args->db = db;
    pt_args->starti = 0;
    pt_args->endi = n_ultra_long_reads;
    pt_args->ultra_long_reads=ultra_long_reads;

    pthread_t tid;    
    int ret = pthread_create(&tid, NULL, align_cudb,(void*)(pt_args));
    NEG_CHK(ret);

    return tid;
}

void align_cudb_async_join(pthread_arg_t *pt_args, pthread_t tid) {
    int ret = pthread_join(tid, NULL);
    NEG_CHK(ret);
    free(pt_args);


}

void align_cuda(core_t* core, db_t* db) {
    int32_t i,j;
    int32_t n_bam_rec = db->n_bam_rec;
    int32_t n_bam_rec_cuda;
    double realtime1;
    int32_t n_ultra_long_reads=0;
    int32_t ultra_long_reads[n_bam_rec];

    //cpu temp pointers
    int32_t* read_ptr_host;
    int32_t* n_events_host;
    int32_t* event_ptr_host;
    event_t* event_table_host;
    AlignedPair* event_align_pairs_host;
    int32_t* read_len_host;
    scalings_t* scalings_host;
    int32_t* n_event_align_pairs_host;
    char* read_host;

    /**cuda pointers*/
    char* read;        //flattened reads sequences
    int32_t* read_ptr; //index pointer for flattedned "reads"
    int32_t* read_len;
    int64_t sum_read_len;
    int32_t* n_events;
    event_t* event_table;
    int32_t* event_ptr;
    int64_t sum_n_events;
    scalings_t* scalings;
    AlignedPair* event_align_pairs;
    int32_t* n_event_align_pairs;
    float *bands;
    uint8_t *trace;
    EventKmerPair* band_lower_left;
    model_t* model_kmer_cache;
    model_t* model;

realtime1 = realtime();

    read_ptr_host = core->cuda->read_ptr_host;

    sum_read_len = 0;

    //read sequences : needflattening
    for (i = 0,j=0; i < n_bam_rec; i++) {
        if(db->read_len[i]<GPU_MAX_READ_LEN){
            read_ptr_host[j] = sum_read_len;
            sum_read_len += (db->read_len[i] + 1); //with null term
            j++;
        }
        else{
            ultra_long_reads[n_ultra_long_reads]=i;
            n_ultra_long_reads++;
        }
    }
    n_bam_rec_cuda = j;
    

    //can start processing on the ultra long reads on the CPU
    pthread_arg_t *tmparg=NULL;
    pthread_t tid =  align_cudb_async(tmparg,core, db, ultra_long_reads, n_ultra_long_reads);
    

    read_len_host = core->cuda->read_len_host;
    scalings_host = core->cuda->scalings_host;
    n_event_align_pairs_host = core->cuda->n_event_align_pairs_host;

    //form the temporary flattened array on host
    read_host = (char*)malloc(sizeof(char) * sum_read_len);
    MALLOC_CHK(read_host);
    for (i = 0,j=0; i < n_bam_rec; i++) {
        if(db->read_len[i]<GPU_MAX_READ_LEN){
            int32_t idx = read_ptr_host[j];
            strcpy(&read_host[idx], db->read[i]);
            read_len_host[j]=db->read_len[i];
            scalings_host[j]=db->scalings[i];
            j++;
        }

    }


    //now the events : need flattening
    //num events : need flattening
    //get the total size and create the pointers
    n_events_host = core->cuda->n_events_host;
    event_ptr_host = core->cuda->event_ptr_host;

    sum_n_events = 0;
    for (i = 0,j=0; i < n_bam_rec; i++) {
        if(db->read_len[i]<GPU_MAX_READ_LEN){
            n_events_host[j] = db->et[i].n;
            event_ptr_host[j] = sum_n_events;
            sum_n_events += db->et[i].n;
            j++;
        }
    }

    //event table flatten
    //form the temporary flattened array on host
    event_table_host =
        (event_t*)malloc(sizeof(event_t) * sum_n_events);
    MALLOC_CHK(event_table_host);
    for (i = 0,j=0; i < n_bam_rec; i++) {
        if(db->read_len[i]<GPU_MAX_READ_LEN){
            int32_t idx = event_ptr_host[j];
            memcpy(&event_table_host[idx], db->et[i].event,
                sizeof(event_t) * db->et[i].n);
                j++;
            }
    }

    event_align_pairs_host =
        (AlignedPair*)malloc(2 * sum_n_events * sizeof(AlignedPair));
    MALLOC_CHK(event_align_pairs_host);

core->align_cuda_preprocess += (realtime() - realtime1);

    /** Start GPU mallocs**/
realtime1 = realtime();


    read_ptr =core->cuda->read_ptr;
    read_len=core->cuda->read_len;
    n_events=core->cuda->n_events;
    event_ptr=core->cuda->event_ptr;
    scalings=core->cuda->scalings;
    model = core->cuda->model;


    print_size("read array",sum_read_len * sizeof(char));
    cudaMalloc((void**)&read, sum_read_len * sizeof(char)); //with null char
    CUDA_CHK();
    print_size("event table",sum_n_events * sizeof(event_t));
    cudaMalloc((void**)&event_table, sum_n_events * sizeof(event_t));
    CUDA_CHK();
    cudaMalloc((void**)&model_kmer_cache, sum_read_len * sizeof(model_t)); 
    CUDA_CHK();
 
    /**allocate output arrays for cuda**/
    print_size("event align pairs",2 * sum_n_events *sizeof(AlignedPair));
    cudaMalloc((void**)&event_align_pairs,
            2 * sum_n_events *
                sizeof(AlignedPair)); //todo : need better huristic
    CUDA_CHK();
    n_event_align_pairs=core->cuda->n_event_align_pairs;

    //scratch arrays
    size_t sum_n_bands = sum_n_events + sum_read_len; //todo : can be optimised 
    print_size("bands",sizeof(float) * sum_n_bands * ALN_BANDWIDTH);
    cudaMalloc((void**)&bands,sizeof(float) * sum_n_bands * ALN_BANDWIDTH);
    CUDA_CHK();
    print_size("trace",sizeof(uint8_t) * sum_n_bands * ALN_BANDWIDTH);
    cudaMalloc((void**)&trace, sizeof(uint8_t) * sum_n_bands * ALN_BANDWIDTH);
    CUDA_CHK();
    cudaMemset(trace,0,sizeof(uint8_t) * sum_n_bands * ALN_BANDWIDTH); //initialise the trace array to 0
    print_size("band_lower_left",sizeof(EventKmerPair)* sum_n_bands);
    cudaMalloc((void**)&band_lower_left, sizeof(EventKmerPair)* sum_n_bands);
    CUDA_CHK();   
core->align_cuda_malloc += (realtime() - realtime1);

    /* cuda mem copys*/
realtime1 =realtime();
    cudaMemcpy(read_ptr, read_ptr_host, n_bam_rec_cuda * sizeof(int32_t),
               cudaMemcpyHostToDevice);
    CUDA_CHK();
    cudaMemcpy(read, read_host, sum_read_len * sizeof(char),
               cudaMemcpyHostToDevice);
    CUDA_CHK();
    //read length : already linear hence direct copy
    cudaMemcpy(read_len, read_len_host, n_bam_rec_cuda * sizeof(int32_t),
               cudaMemcpyHostToDevice);
    CUDA_CHK();
    cudaMemcpy(n_events, n_events_host, n_bam_rec_cuda * sizeof(int32_t),
               cudaMemcpyHostToDevice);
    CUDA_CHK();
    cudaMemcpy(event_ptr, event_ptr_host, n_bam_rec_cuda * sizeof(int32_t),
               cudaMemcpyHostToDevice);
    CUDA_CHK();
    cudaMemcpy(event_table, event_table_host, sizeof(event_t) * sum_n_events,
               cudaMemcpyHostToDevice);
    CUDA_CHK();


    //can be interleaved
    cudaMemcpy(scalings, scalings_host, sizeof(scalings_t) * n_bam_rec_cuda,
               cudaMemcpyHostToDevice);
    CUDA_CHK();
core->align_cuda_memcpy += (realtime() - realtime1);



realtime1 = realtime();    
 
    /*pre kernel*/
    assert(BLOCK_LEN_BANDWIDTH>=ALN_BANDWIDTH);
    dim3 gridpre(1,(n_bam_rec_cuda + BLOCK_LEN_READS - 1) / BLOCK_LEN_READS);
    dim3 blockpre(BLOCK_LEN_BANDWIDTH,BLOCK_LEN_READS);  
	fprintf(stderr,"grid %d,%d, block %d,%d\n",gridpre.x,gridpre.y, blockpre.x,blockpre.y);	

    align_kernel_pre_2d<<<gridpre, blockpre>>>( read,
        read_len, read_ptr, n_events,
        event_ptr, model, n_bam_rec_cuda, model_kmer_cache,bands,trace,band_lower_left); 
       
    cudaDeviceSynchronize();CUDA_CHK();
    fprintf(stderr, "[%s::%.3f*%.2f] align pre done\n", __func__,
            realtime() - realtime1, cputime() / (realtime() - realtime1));
core->align_kernel_time += (realtime() - realtime1);        
core->align_pre_kernel_time += (realtime() - realtime1);        
                
realtime1 = realtime();

    /* core kernel*/
    assert(BLOCK_LEN_BANDWIDTH>=ALN_BANDWIDTH);
    dim3 grid1(1,(n_bam_rec_cuda + BLOCK_LEN_READS - 1) / BLOCK_LEN_READS);
    dim3 block1(BLOCK_LEN_BANDWIDTH,BLOCK_LEN_READS);
    align_kernel_core_2d_shm<<<grid1, block1>>>(read_len, read_ptr, event_table, n_events,
            event_ptr, scalings, n_bam_rec_cuda, model_kmer_cache,bands,trace,band_lower_left );
    
    cudaDeviceSynchronize();CUDA_CHK();
    fprintf(stderr, "[%s::%.3f*%.2f] align done\n", __func__,
    realtime() - realtime1, cputime() / (realtime() - realtime1));
    core->align_kernel_time += (realtime() - realtime1);
core->align_core_kernel_time += (realtime() - realtime1);

realtime1 = realtime();
    
    /*post kernel*/
    int32_t BLOCK_LEN = core->opt.cuda_block_size;
    dim3 gridpost((n_bam_rec_cuda + BLOCK_LEN - 1) / BLOCK_LEN);
    dim3 blockpost(BLOCK_LEN);
    #ifndef WARP_HACK  
        align_kernel_post<<<gridpost, blockpost>>>(event_align_pairs, n_event_align_pairs,
            read_len, read_ptr, event_table, n_events,
            event_ptr,scalings, n_bam_rec_cuda, model_kmer_cache,bands,trace,band_lower_left );     

    #else
        assert(BLOCK_LEN>=32);    
        dim3 grid1post((n_bam_rec_cuda + (BLOCK_LEN/32) - 1) / (BLOCK_LEN/32)); 
        fprintf(stderr,"grid new %d\n",grid1post.x);   
        align_kernel_post<<<grid1post, blockpost>>>(event_align_pairs, n_event_align_pairs,
            read_len, read_ptr, event_table, n_events,
            event_ptr, scalings, n_bam_rec_cuda, model_kmer_cache,bands,trace,band_lower_left );
    #endif
    cudaDeviceSynchronize();CUDA_CHK();
    fprintf(stderr, "[%s::%.3f*%.2f] align post done\n", __func__,
            realtime() - realtime1, cputime() / (realtime() - realtime1));
    core->align_kernel_time += (realtime() - realtime1);        
core->align_post_kernel_time += (realtime() - realtime1);        


    //fprintf(stderr,"readlen %d,n_events %d\n",db->read_len[i],n_event_align_pairs);

#ifdef CUDA_DEBUG

    cudaDeviceSynchronize();
    cudaError_t code = cudaGetLastError();
    //todo : print a message to detect the launch timed out
    if (code == cudaErrorLaunchTimeout) {
        ERROR("%s", "The kernel timed out. You have to first disable the cuda "
                    "time out.");
        fprintf(
            stderr,
            "On Ubuntu do the following\nOpen the file /etc/X11/xorg.conf\nYou "
            "will have a section about your NVIDIA device. Add the following "
            "line to it.\nOption \"Interactive\" \"0\"\nIf you do not have a "
            "section about your NVIDIA device in /etc/X11/xorg.conf or you do "
            "not have a file named /etc/X11/xorg.conf, run the command sudo "
            "nvidia-xconfig to generate a xorg.conf file and do as above.\n\n");
    }
    if (code != cudaSuccess) {
        fprintf(stderr, "Cuda error: %s \n in file : %s line number : %lu\n",
                cudaGetErrorString(code), __FILE__, __LINE__);
        exit(-1);
    }        
    
#endif

    /** copyback ans**/
realtime1 =  realtime();
    cudaMemcpy(n_event_align_pairs_host, n_event_align_pairs,
               n_bam_rec_cuda * sizeof(int32_t), cudaMemcpyDeviceToHost);
    CUDA_CHK();
     
    cudaMemcpy(event_align_pairs_host, event_align_pairs,
               2 * sum_n_events * sizeof(AlignedPair), cudaMemcpyDeviceToHost);
    CUDA_CHK();
core->align_cuda_memcpy += (realtime() - realtime1);

realtime1 =  realtime();

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
    for (i = 0,j=0; i < n_bam_rec; i++) {
        if(db->read_len[i]<GPU_MAX_READ_LEN){
            int32_t idx = event_ptr_host[j];
            db->n_event_align_pairs[i]=n_event_align_pairs_host[j];
            memcpy(db->event_align_pairs[i], &event_align_pairs_host[idx * 2],
                sizeof(AlignedPair) * db->n_event_align_pairs[i]);
                j++;
        }
    }

    //free the temp arrays on host
    free(read_host);
    free(event_table_host);
    free(event_align_pairs_host);

core->align_cuda_postprocess += (realtime() - realtime1);

realtime1 =  realtime(); 
    align_cudb_async_join(tmparg,tid);  
core->extra_load_cpu += (realtime() - realtime1);


}


#endif