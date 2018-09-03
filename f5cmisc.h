#ifndef F5CMISC_H
#define F5CMISC_H

#include "error.h"
#include "f5c.h"
#include <errno.h>
#include <sys/resource.h>
#include <sys/time.h>

event_table getevents(size_t nsample, float* rawptr);
void read_model(model_t* model, const char* file);
void set_model(model_t* model);
void set_cpgmodel(model_t* model);
scalings_t estimate_scalings_using_mom(char* sequence, int32_t sequence_len, model_t* pore_model,
                                       event_table et);
int32_t align(AlignedPair* out_2, char* sequence, int32_t sequence_len,event_table events, model_t* models,
                   scalings_t scaling, float sample_rate);
int32_t postalign(event_alignment_t* alignment, char* sequence, int32_t n_kmers,AlignedPair* event_alignment,
                             int32_t n_events);
bool recalibrate_model(model_t* pore_model, event_table et,
                       scalings_t* scallings,
                       const event_alignment_t* alignment_output,
                       int32_t num_alignments, bool scale_var);
// taken from minimap2/misc
static inline double realtime(void) {
    struct timeval tp;
    struct timezone tzp;
    gettimeofday(&tp, &tzp);
    return tp.tv_sec + tp.tv_usec * 1e-6;
}

// taken from minimap2/misc
static inline double cputime(void) {
    struct rusage r;
    getrusage(RUSAGE_SELF, &r);
    return r.ru_utime.tv_sec + r.ru_stime.tv_sec +
           1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

#endif
