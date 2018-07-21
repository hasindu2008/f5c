#ifndef F5CMISC_H
#define F5CMISC_H
 
#include <errno.h> 
#include <sys/time.h>
#include <sys/resource.h>
 
#define WARN "[WARNING]\033[1;33m "
#define ERR "[ERROR]\033[1;31m "
#define CEND	"\033[0m\n"
 

#define MALLOC_CHK(ret) ({\
    if (ret==NULL){ \
        fprintf(stderr,"[ERROR]\033[1;31m Failed to allocate memory.\033[0m\n\tError at File %s line number %d : %s\n",__FILE__, __LINE__,strerror(errno));\
        exit(EXIT_FAILURE);\
    }\
    }) 
 
/*********************** Some error checks *********************/
/*Die on error. Print the error and exit if the return value of the previous function NULL*/
#define errorCheckNULL(ret) ({\
    if (ret==NULL){ \
        fprintf(stderr,"[ERROR]\033[1;31m Error at File %s line number %d : %s\033[0m\n",__FILE__, __LINE__,strerror(errno));\
        exit(EXIT_FAILURE);\
    }\
    })

/*Die on error. Print the error and exit if the return value of the previous function is -1*/
#define errorCheck(ret) ({\
    if (ret<0){ \
        fprintf(stderr,"[ERROR]\033[1;31m Error at File %s line number %d : %s\033[0m\n",__FILE__, __LINE__,strerror(errno));\
        exit(EXIT_FAILURE);\
    }\
    })

event_table getevents(int nsample,float *rawptr);	
	
	
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