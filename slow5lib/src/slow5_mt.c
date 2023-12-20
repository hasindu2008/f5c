/* @file slow5_mt.c
**
** @@
******************************************************************************/
#ifdef SLOW5_ENABLE_MT

#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <slow5/slow5.h>
#include <slow5/slow5_mt.h>

#define SLOW5_WORK_STEAL 1 //simple work stealing enabled or not (no work stealing mean no load balancing)
#define SLOW5_STEAL_THRESH 1 //stealing threshold

extern enum slow5_log_level_opt  slow5_log_level;
extern enum slow5_exit_condition_opt  slow5_exit_condition;

void *slow5_get_next_mem(size_t *n, const slow5_file_t *s5p);
int slow5_rec_depress_parse(char **mem, size_t *bytes, const char *read_id, slow5_rec_t **read, slow5_file_t *s5p);

#define SLOW5_MALLOC_CHK_LAZY_EXIT(ret) { \
    SLOW5_MALLOC_CHK(ret) \
    if (ret == NULL) { \
        exit(EXIT_FAILURE); \
    } \
}


/* argument wrapper for the multithreaded framework used for data processing */
typedef struct {
    slow5_mt_t* core;
    slow5_batch_t* db;
    int32_t starti;
    int32_t endi;
    void (*func)(slow5_mt_t*,slow5_batch_t*,int);
    int32_t thread_index;
#ifdef SLOW5_WORK_STEAL
    void *all_pthread_args;
#endif
} slow5_pt_arg_t;


/* initialise the core data structure */
slow5_mt_t *slow5_init_mt(int num_thread, slow5_file_t *s5p) {

    slow5_mt_t* core = (slow5_mt_t*)malloc(sizeof(slow5_mt_t));
    SLOW5_MALLOC_CHK_LAZY_EXIT(core);

    core->sf = s5p;
    core->num_thread = num_thread;

    return core;
}

/* free the core data structure */
void slow5_free_mt(slow5_mt_t* core) {
    free(core);
}

/* initialise a data batch */
slow5_batch_t* slow5_init_batch(int batch_capacity){
    slow5_batch_t* db = (slow5_batch_t*)(malloc(sizeof(slow5_batch_t)));
    SLOW5_MALLOC_CHK_LAZY_EXIT(db);

    db->capacity_rec = batch_capacity;
    db->n_rec = 0;

    db->mem_records = (char**)(calloc(db->capacity_rec,sizeof(char*)));
    SLOW5_MALLOC_CHK_LAZY_EXIT(db->mem_records);
    db->mem_bytes = (size_t*)(calloc(db->capacity_rec,sizeof(size_t)));
    SLOW5_MALLOC_CHK_LAZY_EXIT(db->mem_bytes);

    db->slow5_rec = (slow5_rec_t**)calloc(db->capacity_rec,sizeof(slow5_rec_t*));
    SLOW5_MALLOC_CHK_LAZY_EXIT(db->slow5_rec);

    return db;
}

/* load a data batch from disk */
static int slow5_load_db(slow5_mt_t* core, slow5_batch_t* db) {

    int32_t i = 0;
    while (i < db->n_rec) {

        db->mem_records[i] = (char *)slow5_get_next_mem(&(db->mem_bytes[i]), core->sf);

        if (db->mem_records[i] == NULL) {
            if (slow5_errno != SLOW5_ERR_EOF) {
                SLOW5_ERROR("Error reading from SLOW5 file %d\n", slow5_errno);
                exit(EXIT_FAILURE);
            }
            else {
                SLOW5_LOG_DEBUG("%s","Last Batch!\n");
                break;
            }
        }
        else {
            i++;
        }
    }

    return i;
}


static int slow5_write_db(slow5_mt_t* core, slow5_batch_t* db) {

    int32_t i = 0;
    for(i=0;i<db->n_rec;i++) {
        size_t n = fwrite(db->mem_records[i], db->mem_bytes[i], 1, core->sf->fp);
        if (n != 1) {
            SLOW5_ERROR("Writing failed for read id %s!\n", db->slow5_rec[i]->read_id);
        }
    }

    return i;
}


static void slow5_parse_single(slow5_mt_t* core,slow5_batch_t* db, int32_t i){

    assert(db->mem_bytes[i]>0);
    assert(db->mem_records[i]!=NULL);
    int ret=slow5_rec_depress_parse(&db->mem_records[i], &db->mem_bytes[i], NULL, &db->slow5_rec[i], core->sf);
    if(ret!=0){
        SLOW5_ERROR("Error parsing the record %s",db->slow5_rec[i]->read_id);
        exit(EXIT_FAILURE);
    }

}


static void slow5_work_per_single_read(slow5_mt_t* core,slow5_batch_t* db, int32_t i){
    slow5_parse_single(core,db,i);
}

static void slow5_work_per_single_read2(slow5_mt_t* core,slow5_batch_t* db, int32_t i){
    assert(db->rid[i]!=NULL);
    int ret = slow5_get(db->rid[i],&db->slow5_rec[i], core->sf);
    if(ret<0){
        SLOW5_ERROR("Error when fetching the read %s\n",db->rid[i]);
        exit(EXIT_FAILURE);
    }
    db->mem_bytes[i]=ret;

}

static void slow5_work_per_single_read3(slow5_mt_t* core,slow5_batch_t* db, int32_t i){
    assert(db->slow5_rec[i]!=NULL);
    slow5_file_t *sf = core->sf;
    //fprintf(stderr,"Here %d\n",i);
    slow5_press_t *press_ptr = NULL;

    if(sf->compress){
        assert(sf->compress->record_press!=NULL);
        assert(sf->compress->signal_press!=NULL);

        slow5_press_method_t press_out = {sf->compress->record_press->method, sf->compress->signal_press->method};
        press_ptr = slow5_press_init(press_out);
        if(!press_ptr){
            SLOW5_ERROR("Could not initialize the slow5 compression method%s","");
            exit(EXIT_FAILURE);
        }
    }

    //TODO: check if ASCII if press_ptr is still NULL

    db->mem_records[i] = slow5_rec_to_mem(db->slow5_rec[i], sf->header->aux_meta, sf->format, press_ptr, &(db->mem_bytes[i]));
    //fprintf(stderr,"Here 2 %d\n",i);
    slow5_press_free(press_ptr);

    if(db->mem_records[i] == NULL){
        SLOW5_ERROR("Error when converting the read %d to memory\n",i);
        exit(EXIT_FAILURE);
    }

}


/* partially free a data batch - only the read dependent allocations are freed */
static void slow5_free_db_tmp(slow5_batch_t* db) {
    int32_t i = 0;
    for (i = 0; i < db->n_rec; ++i) {
        free(db->mem_records[i]);
    }
}

/* completely free a data batch */
static void slow5_free_db(slow5_batch_t* db) {

    free(db->mem_records);
    free(db->mem_bytes);;

    free(db);
}


static inline int32_t steal_work(slow5_pt_arg_t* all_args, int32_t num_thread) {
	int32_t i, c_i = -1;
	int32_t k;
	for (i = 0; i < num_thread; ++i){
        slow5_pt_arg_t args = all_args[i];
        //fprintf(stderr,"endi : %d, starti : %d\n",args.endi,args.starti);
		if (args.endi-args.starti > SLOW5_STEAL_THRESH) {
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


static void* slow5_pthread_single(void* voidargs) {
    int32_t i;
    slow5_pt_arg_t* args = (slow5_pt_arg_t*)voidargs;
    slow5_batch_t* db = args->db;
    slow5_mt_t* core = args->core;

#ifndef SLOW5_WORK_STEAL
    for (i = args->starti; i < args->endi; i++) {
        args->func(core,db,i);
    }
#else
    slow5_pt_arg_t* all_args = (slow5_pt_arg_t*)(args->all_pthread_args);
    //adapted from kthread.c in minimap2
    for (;;) {
		i = __sync_fetch_and_add(&args->starti, 1);
		if (i >= args->endi) {
            break;
        }
		args->func(core,db,i);
	}
	while ((i = steal_work(all_args,core->num_thread)) >= 0){
		args->func(core,db,i);
    }
#endif

    //fprintf(stderr,"Thread %d done\n",(myargs->position)/THREADS);
    pthread_exit(0);
}

static void slow5_pthread_db(slow5_mt_t* core, slow5_batch_t* db, void (*func)(slow5_mt_t*,slow5_batch_t*,int)){
    //create threads
    pthread_t tids[core->num_thread];
    slow5_pt_arg_t pt_args[core->num_thread];
    int32_t t, ret;
    int32_t i = 0;
    int32_t num_thread = core->num_thread;
    int32_t step = (db->n_rec + num_thread - 1) / num_thread;
    //todo : check for higher num of threads than the data
    //current works but many threads are created despite

    SLOW5_LOG_DEBUG("Creating %d threads\n",num_thread);
    //set the data structures
    for (t = 0; t < num_thread; t++) {
        pt_args[t].core = core;
        pt_args[t].db = db;
        pt_args[t].starti = i;
        i += step;
        if (i > db->n_rec) {
            pt_args[t].endi = db->n_rec;
        } else {
            pt_args[t].endi = i;
        }
        pt_args[t].func=func;
    #ifdef SLOW5_WORK_STEAL
        pt_args[t].all_pthread_args =  (void *)pt_args;
    #endif
        //fprintf(stderr,"t%d : %d-%d\n",t,pt_args[t].starti,pt_args[t].endi);

    }

    //create threads
    for(t = 0; t < core->num_thread; t++){
        ret = pthread_create(&tids[t], NULL, slow5_pthread_single,
                                (void*)(&pt_args[t]));
        if(ret < 0){
            SLOW5_ERROR("Error creating thread %d\n",t);
            exit(EXIT_FAILURE);
        }
    }

    //pthread joining
    for (t = 0; t < core->num_thread; t++) {
        int ret = pthread_join(tids[t], NULL);
        if(ret < 0){
            SLOW5_ERROR("Error creating thread %d\n",t);
            exit(EXIT_FAILURE);
        }
    }
}

/* process all reads in the given batch db */
static void slow5_work_db(slow5_mt_t* core, slow5_batch_t* db, void (*func)(slow5_mt_t*,slow5_batch_t*,int)){

    if (core->num_thread == 1) {
        int32_t i=0;
        for (i = 0; i < db->n_rec; i++) {
            func(core,db,i);
        }

    }

    else {
        slow5_pthread_db(core,db,func);
    }
}

int slow5_get_batch(slow5_mt_t *core, slow5_batch_t *db, char **rid, int num_rid){

    if(num_rid>db->capacity_rec){
        SLOW5_ERROR("Requested %d is greater than the capacity %d",num_rid,db->capacity_rec);
        exit(EXIT_FAILURE);
    }

    db->rid = rid;
    db->n_rec = num_rid;
    slow5_work_db(core,db,slow5_work_per_single_read2);
    SLOW5_LOG_DEBUG("loaded and parsed %d recs\n",num_rid);
    slow5_free_db_tmp(db);

    return num_rid;
}


int slow5_get_next_batch(slow5_mt_t *core, slow5_batch_t *db, int num_reads){

    if(num_reads>db->capacity_rec){
        SLOW5_ERROR("Requested %d is greater than the capacity %d",num_reads,db->capacity_rec);
        exit(EXIT_FAILURE);
    }
    db->n_rec = num_reads;
    int n = slow5_load_db(core,db);
    db->n_rec = n;
    SLOW5_LOG_DEBUG("Loaded %d recs\n",n);
    slow5_work_db(core,db,slow5_work_per_single_read);
    SLOW5_LOG_DEBUG("Parsed %d recs\n",n);
    slow5_free_db_tmp(db);

    return n;
}

int slow5_encode_batch(slow5_mt_t *core, slow5_batch_t *db, int num_reads){
    db->n_rec = num_reads;
    slow5_work_db(core,db,slow5_work_per_single_read3);
    return db->n_rec;
}

int slow5_write_batch(slow5_mt_t *core, slow5_batch_t *db, int num_reads){

    db->n_rec = num_reads;
    slow5_work_db(core,db,slow5_work_per_single_read3);
    SLOW5_LOG_DEBUG("Processed %d recs\n",num_reads);

    int num_wr=slow5_write_db(core,db);
    SLOW5_LOG_DEBUG("Written %d recs\n",num_wr);

    slow5_free_db_tmp(db);

    return num_wr;
}

void slow5_free_batch(slow5_batch_t *db){

    slow5_rec_t **reads = db->slow5_rec;
    if(reads != NULL){
        for(int i=0;i<db->capacity_rec;i++){
            slow5_rec_free(reads[i]);
        }
    }

    free(reads);
    //*read = NULL;

    //slow5_free_db_tmp(db);
    slow5_free_db(db);
}


int slow5_get_next_batch_lazy(slow5_rec_t ***read, slow5_file_t *s5p, int batch_size, int num_threads){

    slow5_mt_t *core = slow5_init_mt(num_threads, s5p);
    slow5_batch_t* db = slow5_init_batch(batch_size);

    int ret = slow5_get_next_batch(core,db,batch_size);

    *read = db->slow5_rec;
    db->slow5_rec = NULL;

    slow5_free_batch(db);
    slow5_free_mt(core);

    return ret;

}
int slow5_get_batch_lazy(slow5_rec_t ***read, slow5_file_t *s5p, char **rid, int num_rid, int num_threads){

    slow5_mt_t *core = slow5_init_mt(num_threads, s5p);
    slow5_batch_t* db = slow5_init_batch(num_rid);

    int ret = slow5_get_batch(core,db,rid,num_rid);

    *read = db->slow5_rec;
    db->slow5_rec = NULL;

    slow5_free_batch(db);
    slow5_free_mt(core);

    return ret;

}
int slow5_write_batch_lazy(slow5_rec_t **read, slow5_file_t *s5p, int batch_size, int num_threads){
    slow5_mt_t *core = slow5_init_mt(num_threads, s5p);
    slow5_batch_t* db = slow5_init_batch(batch_size);


    free(db->slow5_rec);
    db->slow5_rec = read;
    int ret = slow5_write_batch(core,db,batch_size);
    db->slow5_rec = NULL;

    slow5_free_batch(db);
    slow5_free_mt(core);

    return ret;
}

void slow5_free_batch_lazy(slow5_rec_t ***read, int num_rec){

    slow5_rec_t **reads = *read;
    for(int i=0;i<num_rec;i++){
        slow5_rec_free(reads[i]);
    }

    free(reads);
    *read = NULL;
}

#else

#include <stdio.h>
#include <stdlib.h>
#include <slow5/slow5.h>
#include <slow5/slow5_mt.h>

int slow5_get_next_batch_lazy(slow5_rec_t ***read, slow5_file_t *s5p, int batch_size, int num_threads){
    fprintf(stderr,"slow5lib has not been compiled with lazy multithreading support\n");
    exit(EXIT_FAILURE);
}
int slow5_get_batch_lazy(slow5_rec_t ***read, slow5_file_t *s5p, char **rid, int num_rid, int num_threads){
    fprintf(stderr,"slow5lib has not been compiled with lazy multithreading support\n");
    exit(EXIT_FAILURE);
}

int slow5_write_batch_lazy(slow5_rec_t **read, slow5_file_t *s5p, int batch_size, int num_threads){
    fprintf(stderr,"slow5lib has not been compiled with lazy multithreading support\n");
    exit(EXIT_FAILURE);
}
void slow5_free_batch_lazy(slow5_rec_t ***read, int num_rec){
    fprintf(stderr,"slow5lib has not been compiled with lazy multithreading support\n");
    exit(EXIT_FAILURE);
}

slow5_mt_t *slow5_init_mt(int num_thread, slow5_file_t *s5p){
    fprintf(stderr,"slow5lib has not been compiled with lazy multithreading support\n");
    exit(EXIT_FAILURE);
}
slow5_batch_t* slow5_init_batch(int batch_capacity){
    fprintf(stderr,"slow5lib has not been compiled with lazy multithreading support\n");
    exit(EXIT_FAILURE);
}
int slow5_get_next_batch(slow5_mt_t *mt, slow5_batch_t *read_batch, int num_reads){
    fprintf(stderr,"slow5lib has not been compiled with lazy multithreading support\n");
    exit(EXIT_FAILURE);
}
int slow5_get_batch(slow5_mt_t *mt, slow5_batch_t *read_batch, char **rid, int num_rid){
    fprintf(stderr,"slow5lib has not been compiled with lazy multithreading support\n");
    exit(EXIT_FAILURE);
}
int slow5_encode_batch(slow5_mt_t *core, slow5_batch_t *db, int num_reads){
    fprintf(stderr,"slow5lib has not been compiled with lazy multithreading support\n");
    exit(EXIT_FAILURE);
}
int slow5_write_batch(slow5_mt_t *mt, slow5_batch_t *read_batch, int num_reads){
    fprintf(stderr,"slow5lib has not been compiled with lazy multithreading support\n");
    exit(EXIT_FAILURE);
}
void slow5_free_batch(slow5_batch_t *read_batch){
    fprintf(stderr,"slow5lib has not been compiled with lazy multithreading support\n");
    exit(EXIT_FAILURE);
}
void slow5_free_mt(slow5_mt_t *mt){
    fprintf(stderr,"slow5lib has not been compiled with lazy multithreading support\n");
    exit(EXIT_FAILURE);
}


#endif