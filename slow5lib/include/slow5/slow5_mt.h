/* @file slow5_mt.h
**
******************************************************************************/
#ifndef SLOW5_MT_H
#define SLOW5_MT_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/**************************************************************************************************
 ***  Easy Multi-thread API *************************************************************************
 **************************************************************************************************/

/*
This is a easy multi-thread API that can fetch a batch of slow5 records using multiple threads in parallel.
This API uses a fork-join thread model. It is not meant to be used by a programmer who has the expertise to write multi-threaded code and use the slow5 low-level API directly.
While fetching using this API would be faster than using a single thread, it would not be as efficient as using low-level API functions
directly from the programmer's own multi-threaded code. This API can be optimised using a thread pool, but that is not implemented yet.
*/

/* a batch of read data (dynamic data based on the reads) */
typedef struct {
    int32_t n_rec;
    int32_t capacity_rec;

    char **mem_records; //unused in get()
    size_t *mem_bytes;

    slow5_rec_t **slow5_rec;
    char **rid; //only used in get()

} slow5_batch_t;

/* mt data structure (mostly static data throughout the program lifetime) */
typedef struct {
    //slow5
    slow5_file_t *sf;
    int num_thread;
} slow5_mt_t;

/*
these functions will lazily exit on error (need to do proper error handling, but a bit too much work at the moment)
also these functions are not optimised for cases that are unlikely to be bottlenecks
that is they do superfluous mallocs and free and computations in cases which are unlikely to be bottlenecks
also each batch call will create and destruct threads rather than using a thread pool
*/

slow5_mt_t *slow5_init_mt(int num_thread, slow5_file_t *s5p);
slow5_batch_t* slow5_init_batch(int batch_capacity);
int slow5_get_next_batch(slow5_mt_t *mt, slow5_batch_t *read_batch, int num_reads);
int slow5_get_batch(slow5_mt_t *mt, slow5_batch_t *read_batch, char **rid, int num_rid);
int slow5_write_batch(slow5_mt_t *mt, slow5_batch_t *read_batch, int num_reads);
int slow5_encode_batch(slow5_mt_t *core, slow5_batch_t *db, int num_reads); //used by slow5curl
void slow5_free_batch(slow5_batch_t *read_batch);
void slow5_free_mt(slow5_mt_t *mt);

/*
These functions are for those who are even more lazier to init and free the slow5_mt_t and slow5_batch_t structures.
Mainly for use by the pyslow5 python module. Recommended to use above functions instead.
*/
int slow5_get_next_batch_lazy(slow5_rec_t ***read, slow5_file_t *s5p, int batch_size, int num_threads);
int slow5_get_batch_lazy(slow5_rec_t ***read, slow5_file_t *s5p, char **rid, int num_rid, int num_threads);
int slow5_write_batch_lazy(slow5_rec_t **read, slow5_file_t *s5p, int batch_size, int num_threads);
void slow5_free_batch_lazy(slow5_rec_t ***read, int num_rec);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
