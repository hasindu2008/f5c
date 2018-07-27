#ifndef F5CMISC_H
#define F5CMISC_H

#include <errno.h>
#include <sys/resource.h>
#include <sys/time.h>

#include "error.h"

event_table getevents(size_t nsample, float* rawptr);

//taken from minimap2/misc
static inline double realtime(void) {
    struct timeval tp;
    struct timezone tzp;
    gettimeofday(&tp, &tzp);
    return tp.tv_sec + tp.tv_usec * 1e-6;
}

//taken from minimap2/misc
static inline double cputime(void) {
    struct rusage r;
    getrusage(RUSAGE_SELF, &r);
    return r.ru_utime.tv_sec + r.ru_stime.tv_sec +
           1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

#endif
