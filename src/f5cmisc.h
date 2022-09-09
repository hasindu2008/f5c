/* @file f5cmisc.h
**
** miscellaneous definitions and function prototypes to f5c framework
** @author: Hasindu Gamaarachchi (hasindu@unsw.edu.au)
** @@
******************************************************************************/

#ifndef F5CMISC_H
#define F5CMISC_H

#include <sys/resource.h>
#include <sys/time.h>
#include "f5c.h"
#include "error.h"

#define MIN_CALIBRATION_VAR 2.5
#define MAX_EVENT_TO_BP_RATIO 20
#define AVG_EVENTS_PER_KMER_MAX 15.0f //if average events per base of a read >AVG_EVENTS_PER_KMER_MAX do not process

//model types
#define MODEL_TYPE_NUCLEOTIDE 1
#define MODEL_TYPE_METH 2

//default model IDs
#define MODEL_ID_DNA_NUCLEOTIDE 1
#define MODEL_ID_DNA_CPG 2
#define MODEL_ID_RNA_NUCLEOTIDE 3


// Flags to modify the behaviour of the HMM
enum HMMAlignmentFlags
{
    HAF_ALLOW_PRE_CLIP = 1, // allow events to go unmatched before the aligning region
    HAF_ALLOW_POST_CLIP = 2 // allow events to go unmatched after the aligning region
};

/* performance related parameter profiles for various systems*/
int set_profile(char* profile, opt_t *opt);

/* I/O processes for fast5 reading*/
void init_iop(core_t* core,opt_t opt);
void free_iop(core_t* core,opt_t opt);

/* input */
ret_status_t load_db1(core_t* core, db_t* db);
ret_status_t load_db2(core_t* core, db_t* db);

/* models */
uint32_t read_model(model_t* model, const char* file, uint32_t type);
uint32_t set_model(model_t* model, uint32_t model_id);

/* events */
event_table getevents(size_t nsample, float* rawptr, int8_t rna);

/* alignment related */
scalings_t estimate_scalings_using_mom(char* sequence, int32_t sequence_len, model_t* pore_model, uint32_t kmer_size, event_table et);
int32_t align(AlignedPair* out_2, char* sequence, int32_t sequence_len,
              event_table events, model_t* models, uint32_t kmer_size, scalings_t scaling,
              float sample_rate);
int32_t postalign(event_alignment_t* alignment, index_pair_t* base_to_event_map, double* events_per_base,
                  char* sequence, int32_t n_kmers, AlignedPair* event_alignment, int32_t n_events, uint32_t kmer_size);
bool recalibrate_model(model_t* pore_model, uint32_t kmer_size, event_table et, scalings_t* scallings,
                       const event_alignment_t* alignment_output, int32_t num_alignments, bool scale_var, int32_t minNumEventsToRescale);

/* methylation call */
void calculate_methylation_for_read(std::map<int, ScoredSite>* site_score_map, char* ref, bam1_t* record,
                                    int32_t read_length, event_t* event, index_pair_t* base_to_event_map,
                                    scalings_t scaling, model_t* cpgmodel, uint32_t kmer_size, double events_per_base);

/* event alignment */
char *emit_event_alignment_tsv(uint32_t strand_idx,
                              const event_table* et, model_t* model, uint32_t kmer_size,  scalings_t scalings,
                              const std::vector<event_alignment_t>& alignments,
                              int8_t print_read_names, int8_t scale_events, int8_t write_samples, int8_t write_signal_index, int8_t collapse,
                              int64_t read_index, char* read_name, char *ref_name, float sample_rate, float *samples);
void emit_event_alignment_tsv_header(FILE* fp, int8_t print_read_names, int8_t write_samples, int8_t write_signal_index);
void emit_sam_header(samFile* fp, const bam_hdr_t* hdr);
void emit_event_alignment_sam(htsFile* fp, char* read_name, bam_hdr_t* base_hdr, bam1_t* base_record,
                              const std::vector<event_alignment_t>& alignments);
void realign_read(std::vector<event_alignment_t>* event_alignment_result, EventalignSummary *summary, FILE *summary_fp,char* ref,
                  const bam_hdr_t* hdr, const bam1_t* record, int32_t read_length,
                  size_t read_idx, int region_start, int region_end,
                  event_table* events, model_t* model, uint32_t kmer_size, index_pair_t* base_to_event_map,
                  scalings_t scaling,double events_per_base, float sample_rate);


/* hmm */
float profile_hmm_score(const char *m_seq,const char *m_rc_seq, event_t* event, scalings_t scaling,
                        model_t* cpgmodel, uint32_t kmer_size, uint32_t event_start_idx,
                        uint32_t event_stop_idx, uint8_t strand, int8_t event_stride,
                        uint8_t rc,double events_per_base,uint32_t hmm_flags);
float profile_hmm_score_r9(const char *m_seq,
                            const char *m_rc_seq,
                            event_t* event,
                            scalings_t scaling,
                            model_t* cpgmodel, uint32_t kmer_size,
                            uint32_t event_start_idx,
                            uint32_t event_stop_idx,
                            uint8_t strand,
                            int8_t event_stride,
                            uint8_t rc,
                            double events_per_base,
                            uint32_t hmm_flags);

//void pthread_db(core_t* core, db_t* db, void (*func)(core_t*,db_t*,int));

#ifdef HAVE_CUDA
    /* alignment on GPU */
    void align_cuda(core_t* core, db_t* db);
#endif

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

//taken from minimap2
static inline long peakrss(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
#ifdef __linux__
	return r.ru_maxrss * 1024;
#else
	return r.ru_maxrss;
#endif
}

// Prints to the provided buffer a nice number of bytes (KB, MB, GB, etc)
//from https://www.mbeckler.org/blog/?p=114
static inline void print_size(const char* name, uint64_t bytes)
{
    const char* suffixes[7];
    suffixes[0] = "B";
    suffixes[1] = "KB";
    suffixes[2] = "MB";
    suffixes[3] = "GB";
    suffixes[4] = "TB";
    suffixes[5] = "PB";
    suffixes[6] = "EB";
    uint64_t s = 0; // which suffix to use
    double count = bytes;
    while (count >= 1024 && s < 7)
    {
        s++;
        count /= 1024;
    }
    if (count - floor(count) == 0.0)
        fprintf(stderr, "[%s] %s : %d %s\n", __func__ , name, (int)count, suffixes[s]);
    else
        fprintf(stderr, "[%s] %s : %.1f %s\n", __func__, name, count, suffixes[s]);
}

//replace u with t in a string
static inline void replace_char(char *str, char u, char t){
    while(*str){
        if(*str == u){
            *str = t;
        }
        str++;
    }
}

static inline int64_t mm_parse_num(const char* str) //taken from minimap2
{
    double x;
    char* p;
    x = strtod(str, &p);
    if (*p == 'G' || *p == 'g')
        x *= 1e9;
    else if (*p == 'M' || *p == 'm')
        x *= 1e6;
    else if (*p == 'K' || *p == 'k')
        x *= 1e3;
    return (int64_t)(x + .499);
}


#endif
