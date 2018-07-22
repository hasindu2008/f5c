#ifndef F5CMISC_H
#define F5CMISC_H
 
#include <errno.h> 
#include <sys/time.h>
#include <sys/resource.h>
 
 
#define WARN "[%s::WARNING]\033[1;33m "
#define ERR "[%s::ERROR]\033[1;31m "
#define CEND	"\033[0m\n"
 
#define WARNING(arg, ...) fprintf(stderr,"[%s::WARNING]\033[1;33m " arg "\033[0m\n",__func__, __VA_ARGS__) 
#define ERROR(arg, ...) fprintf(stderr,"[%s::ERROR]\033[1;31m " arg "\033[0m\n",__func__, __VA_ARGS__) 
#define INFO(arg, ...) fprintf(stderr,"[%s::INFO]\033[1;34m " arg "\033[0m\n",__func__, __VA_ARGS__) 
#define SUCCESS(arg, ...) fprintf(stderr,"[%s::SUCCESS]\033[1;32m " arg "\033[0m\n",__func__, __VA_ARGS__) 
#define DEBUG(arg, ...) fprintf(stderr,"[%s::DEBUG]\033[1;35m Error occured at %s:%d. " arg "\033[0m\n",__func__,__FILE__, __LINE__-2, __VA_ARGS__) 

#define MALLOC_CHK(ret) ({\
    if (ret==NULL){ \
        fprintf(stderr,"[%s::ERROR]\033[1;31m Failed to allocate memory : %s.\033[0m\
						\n[%s::DEBUG]\033[1;35m Error occured at %s:%d.\033[0m\n\n" \
						,__func__,strerror(errno),__func__,__FILE__, __LINE__);\
        exit(EXIT_FAILURE);\
    }\
    }) 
 
/*********************** Some error checks *********************/
/*Die on error. Print the error and exit if the return value of the previous function NULL*/
#define NULL_CHK(ret) ({\
    if (ret==NULL){ \
        fprintf(stderr,"[%s::ERROR]\033[1;31m %s.\033[0m\
						\n[%s::DEBUG]\033[1;35m Error occured at %s:%d.\033[0m\n\n" \
						,__func__,strerror(errno),__func__,__FILE__, __LINE__);\
        exit(EXIT_FAILURE);\
    }\
    })

/*Die on error. Print the error and exit if the return value of the previous function is -1*/
#define NEG_CHECK(ret) ({\
    if (ret<0){ \
        fprintf(stderr,"[%s::ERROR]\033[1;31m %s.\033[0m\
						\n[%s::DEBUG]\033[1;35m Error occured at %s:%d.\033[0m\n\n" \
						,__func__,strerror(errno),__func__,__FILE__, __LINE__);\
        exit(EXIT_FAILURE);\
    }\
    })

event_table getevents(size_t nsample,float *rawptr);	
	
	
//taken from minimap2/misc	
static inline double realtime(void)
{
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

//taken from minimap2/misc	
static inline double cputime(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

#endif