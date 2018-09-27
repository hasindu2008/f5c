#include "f5c.h"
#include <assert.h>
#include "f5cmisc.cuh"

#ifdef ALIGN_KERNEL_SLICED

//#define DEBUG_ESTIMATED_SCALING 1
//#define DEBUG_RECALIB_SCALING 1
//#define DEBUG_ADAPTIVE 1


//todo : can make more efficient using bit encoding
//todo : is inlining correct?
__forceinline__ __device__  uint32_t get_rank(char base) {
    if (base == 'A') { //todo: do we neeed simple alpha?
        return 0;
    } else if (base == 'C') {
        return 1;
    } else if (base == 'G') {
        return 2;
    } else if (base == 'T') {
        return 3;
    } else {
        //WARNING("A None ACGT base found : %c", base); //todo : fix this in gpu code
        return 0;
    }
}

// return the lexicographic rank of the kmer amongst all strings of
// length k for this alphabet
__forceinline__ __device__  uint32_t get_kmer_rank(const char* str, uint32_t k) {
    //uint32_t p = 1;
    uint32_t r = 0;

    // from last base to first
    for (uint32_t i = 0; i < k; ++i) {
        //r += rank(str[k - i - 1]) * p;
        //p *= size();
        r += get_rank(str[k - i - 1]) << (i << 1);
    }
    return r;
}

//copy a kmer from a reference
__forceinline__ __device__ void kmer_cpy(char* dest, char* src, uint32_t k) {
    uint32_t i = 0;
    for (i = 0; i < k; i++) {
        dest[i] = src[i];
    }
    dest[i] = '\0';
}

#define log_inv_sqrt_2pi  -0.918938f // Natural logarithm

__forceinline__ __device__ float
log_normal_pdf(float x, float gp_mean, float gp_stdv, float gp_log_stdv) {
    /*INCOMPLETE*/
    //float log_inv_sqrt_2pi = -0.918938f; // Natural logarithm
    float a = (x - gp_mean) / gp_stdv;
    return log_inv_sqrt_2pi - gp_log_stdv + (-0.5f * a * a);
    // return 1;
}

__forceinline__ __device__ float
log_probability_match_r9(scalings_t scaling, model_t* models, event_t* event,
                         int event_idx, uint32_t kmer_rank) {
    // event level mean, scaled with the drift value
    //strand = 0;
 #ifdef DEBUG_ADAPTIVE   
    assert(kmer_rank < 4096);
 #endif   
    //float level = read.get_drift_scaled_level(event_idx, strand);

    //float time =
    //    (events.event[event_idx].start - events.event[0].start) / sample_rate;
    float unscaledLevel = event[event_idx].mean;
    float scaledLevel = unscaledLevel;
    //float scaledLevel = unscaledLevel - time * scaling.shift;

    //fprintf(stderr, "level %f\n",scaledLevel);
    //GaussianParameters gp = read.get_scaled_gaussian_from_pore_model_state(pore_model, strand, kmer_rank);
    float gp_mean =
        scaling.scale * models[kmer_rank].level_mean + scaling.shift;
    float gp_stdv = models[kmer_rank].level_stdv * 1; //scaling.var = 1;
    // float gp_stdv = 0;
    // float gp_log_stdv = models[kmer_rank].level_log_stdv + scaling.log_var;
    // if(models[kmer_rank].level_stdv <0.01 ){
    // 	fprintf(stderr,"very small std dev %f\n",models[kmer_rank].level_stdv);
    // }
    float gp_log_stdv =
        log(models[kmer_rank].level_stdv + 0); // scaling.log_var = log(1)=0;

    float lp = log_normal_pdf(scaledLevel, gp_mean, gp_stdv, gp_log_stdv);
    return lp;
}

#define event_kmer_to_band(ei, ki) (ei + 1) + (ki + 1)
#define band_event_to_offset(bi, ei) band_lower_left[bi].event_idx - (ei)
#define band_kmer_to_offset(bi, ki) (ki) - band_lower_left[bi].kmer_idx
#define is_offset_valid(offset) (offset) >= 0 && (offset) < bandwidth
#define event_at_offset(bi, offset) band_lower_left[(bi)].event_idx - (offset)
#define kmer_at_offset(bi, offset) band_lower_left[(bi)].kmer_idx + (offset)

#define move_down(curr_band)                                                   \
    { curr_band.event_idx + 1, curr_band.kmer_idx }
#define move_right(curr_band)                                                  \
    { curr_band.event_idx, curr_band.kmer_idx + 1 }

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

#define BAND_ARRAY(r, c) ( bands[((r)*(ALN_BANDWIDTH)+(c))] )
#define TRACE_ARRAY(r, c) ( trace[((r)*(ALN_BANDWIDTH)+(c))] )

#define FROM_D  0
#define FROM_U  1
#define FROM_L  2

#define min_average_log_emission  -5.0
#define max_gap_threshold  50
#define bandwidth  ALN_BANDWIDTH
#define half_bandwidth  ALN_BANDWIDTH/2

#define epsilon 1e-10

__global__ void adaptive_align_inner(float *bands,uint8_t *trace,EventKmerPair* band_lower_left,
    int min_offset, int max_offset, int band_idx, int32_t* kmer_ranks,
    scalings_t scaling, model_t* models, event_t* events,double p_stay){
    
    //double epsilon = 1e-10;
    double lp_skip = log(epsilon);
    double lp_stay = log(p_stay);
    double lp_step = log(1.0 - exp(lp_skip) - exp(lp_stay));
    //double lp_trim = log(0.01);

    int tid=blockIdx.x*blockDim.x+threadIdx.x;
    //for (int offset = min_offset; offset < max_offset; ++offset) {
    if(tid < (max_offset-min_offset)) {
        int offset=min_offset+tid;

        int event_idx = event_at_offset(band_idx, offset);
        int kmer_idx = kmer_at_offset(band_idx, offset);

        int32_t kmer_rank = kmer_ranks[kmer_idx];

        int offset_up = band_event_to_offset(band_idx - 1, event_idx - 1);
        int offset_left = band_kmer_to_offset(band_idx - 1, kmer_idx - 1);
        int offset_diag = band_kmer_to_offset(band_idx - 2, kmer_idx - 1);

#ifdef DEBUG_ADAPTIVE
        // verify loop conditions
        assert(kmer_idx >= 0 && kmer_idx < n_kmers);
        assert(event_idx >= 0 && event_idx < n_events);
        assert(offset_diag ==
               band_event_to_offset(band_idx - 2, event_idx - 1));
        assert(offset_up - offset_left == 1);
        assert(offset >= 0 && offset < bandwidth);
#endif //DEBUG_ADAPTIVE

        float up = is_offset_valid(offset_up)
                       ? BAND_ARRAY(band_idx - 1,offset_up)
                       : -INFINITY;
        float left = is_offset_valid(offset_left)
                         ? BAND_ARRAY(band_idx - 1,offset_left)
                         : -INFINITY;
        float diag = is_offset_valid(offset_diag)
                         ? BAND_ARRAY(band_idx - 2,offset_diag)
                         : -INFINITY;

        float lp_emission = log_probability_match_r9(
            scaling, models, events, event_idx, kmer_rank);
        //fprintf(stderr, "lp emiision : %f , event idx %d, kmer rank %d\n", lp_emission,event_idx,kmer_rank);
        float score_d = diag + lp_step + lp_emission;
        float score_u = up + lp_stay + lp_emission;
        float score_l = left + lp_skip;

        float max_score = score_d;
        uint8_t from = FROM_D;

        max_score = score_u > max_score ? score_u : max_score;
        from = max_score == score_u ? FROM_U : from;
        max_score = score_l > max_score ? score_l : max_score;
        from = max_score == score_l ? FROM_L : from;

#ifdef DEBUG_ADAPTIVE
        fprintf(stderr,
                "[adafill] offset-up: %d offset-diag: %d offset-left: %d\n",
                offset_up, offset_diag, offset_left);
        fprintf(stderr, "[adafill] up: %.2lf diag: %.2lf left: %.2lf\n", up,
                diag, left);
        fprintf(stderr,
                "[adafill] bi: %d o: %d e: %d k: %d s: %.2lf f: %d emit: "
                "%.2lf\n",
                band_idx, offset, event_idx, kmer_idx, max_score, from,
                lp_emission);
#endif //DEBUG_ADAPTIVE
        BAND_ARRAY(band_idx,offset) = max_score;
        TRACE_ARRAY(band_idx,offset) = from;
        //fills += 1;
    }

}

__forceinline__ __device__ void align_single_pre(AlignedPair* out_2, char* sequence,
                                int32_t sequence_len, event_t* events,
                                int32_t n_event, model_t* models,
                                scalings_t scaling,int32_t* kmer_ranks,float *bands,uint8_t *trace, EventKmerPair* band_lower_left) {
    //fprintf(stderr, "%s\n", sequence);
    //fprintf(stderr, "Scaling %f %f", scaling.scale, scaling.shift);

    //size_t strand_idx = 0;
    //size_t k = 6;

    // size_t n_events = events[strand_idx].n;
    int32_t n_events = n_event;
    int32_t n_kmers = sequence_len - KMER_SIZE + 1;
    //fprintf(stderr,"n_kmers : %d\n",n_kmers);
    // backtrack markers
    //const uint8_t FROM_D = 0;
    //const uint8_t FROM_U = 1;
    //const uint8_t FROM_L = 2;

    // qc
    //double min_average_log_emission = -5.0;
    //int max_gap_threshold = 50;

    // banding
    //int bandwidth = ALN_BANDWIDTH;
    //half_bandwidth = bandwidth / 2;

    // transition penalties
    float events_per_kmer = (float)n_events / n_kmers;
    float p_stay = 1 - (1 / (events_per_kmer + 1));

    // setting a tiny skip penalty helps keep the true alignment within the adaptive band
    // this was empirically determined
    //double epsilon = 1e-10;
    double lp_skip = log(epsilon);
    double lp_stay = log(p_stay);
    double lp_step = log(1.0 - exp(lp_skip) - exp(lp_stay));
    double lp_trim = log(0.01);

    // dp matrix
    int32_t n_rows = n_events + 1;
    int32_t n_cols = n_kmers + 1;
    int32_t n_bands = n_rows + n_cols;

    // Initialize

    // Precompute k-mer ranks to avoid doing this in the inner loop
    //size_t* kmer_ranks = (size_t*)malloc(sizeof(size_t) * n_kmers);
    //MALLOC_CHK(kmer_ranks); //todo : fix these to error check

    for (int32_t i = 0; i < n_kmers; ++i) {
        //>>>>>>>>> New replacement begin
        char* substring = &sequence[i];
        kmer_ranks[i] = get_kmer_rank(substring, KMER_SIZE);
        //<<<<<<<<< New replacement over
    }

    //float** bands = (float**)malloc(sizeof(float*) * n_bands);
    //MALLOC_CHK(bands);
    //uint8_t** trace = (uint8_t**)malloc(sizeof(uint8_t*) * n_bands);
    //MALLOC_CHK(trace);

    for (int32_t i = 0; i < n_bands; i++) {
        //bands[i] = (float*)malloc(sizeof(float) * bandwidth);
        //MALLOC_CHK(bands[i]);
        //trace[i] = (uint8_t*)malloc(sizeof(uint8_t) * bandwidth);
        //MALLOC_CHK(trace[i]);

        for (int j = 0; j < bandwidth; j++) {
            BAND_ARRAY(i,j) = -INFINITY;
            TRACE_ARRAY(i,j) = 0;
        }
    }

    // Keep track of the event/kmer index for the lower left corner of the band
    // these indices are updated at every iteration to perform the adaptive banding
    // Only the first two bands have their coordinates initialized, the rest are computed adaptively

    // struct EventKmerPair {
    //     int event_idx;
    //     int kmer_idx;
    // };
    //>>>>>>>>>>>>>>>>>New Replacement Begin
    //struct EventKmerPair* band_lower_left =
    //    (struct EventKmerPair*)malloc(sizeof(struct EventKmerPair) * n_bands);
    //MALLOC_CHK(band_lower_left);
    //std::vector<EventKmerPair> band_lower_left(n_bands);
    //<<<<<<<<<<<<<<<<<New Replacement over

    // initialize range of first two bands
    band_lower_left[0].event_idx = half_bandwidth - 1;
    band_lower_left[0].kmer_idx = -1 - half_bandwidth;
    band_lower_left[1] = move_down(band_lower_left[0]);

    int start_cell_offset = band_kmer_to_offset(0, -1);
    assert(is_offset_valid(start_cell_offset));
    assert(band_event_to_offset(0, -1) == start_cell_offset);
    BAND_ARRAY(0,start_cell_offset) = 0.0f;

    // band 1: first event is trimmed
    int first_trim_offset = band_event_to_offset(1, 0);
    assert(kmer_at_offset(1, first_trim_offset) == -1);
    assert(is_offset_valid(first_trim_offset));
    BAND_ARRAY(1,first_trim_offset) = lp_trim;
    TRACE_ARRAY(1,first_trim_offset) = FROM_U;

    //int fills = 0;
#ifdef DEBUG_ADAPTIVE
    fprintf(stderr, "[trim] bi: %d o: %d e: %d k: %d s: %.2lf\n", 1,
            first_trim_offset, 0, -1, BAND_ARRAY(1,first_trim_offset);
#endif

}

__forceinline__ __device__ void align_single(AlignedPair* out_2, char* sequence,
                                int32_t sequence_len, event_t* events,
                                int32_t n_event, model_t* models,
                                scalings_t scaling,int32_t* kmer_ranks,float *bands,uint8_t *trace, EventKmerPair* band_lower_left) {

    //fprintf(stderr, "%s\n", sequence);
    //fprintf(stderr, "Scaling %f %f", scaling.scale, scaling.shift);

    //size_t strand_idx = 0;
    //size_t k = 6;

    // size_t n_events = events[strand_idx].n;
    int32_t n_events = n_event;
    int32_t n_kmers = sequence_len - KMER_SIZE + 1;
    //fprintf(stderr,"n_kmers : %d\n",n_kmers);
    // backtrack markers
    //const uint8_t FROM_D = 0;
    //const uint8_t FROM_U = 1;
    //const uint8_t FROM_L = 2;

    // qc
    //double min_average_log_emission = -5.0;
    //int max_gap_threshold = 50;

    // banding
    //int bandwidth = ALN_BANDWIDTH;
    //half_bandwidth = bandwidth / 2;

    // transition penalties
    float events_per_kmer = (float)n_events / n_kmers;
    float p_stay = 1 - (1 / (events_per_kmer + 1));

    // setting a tiny skip penalty helps keep the true alignment within the adaptive band
    // this was empirically determined
    //double epsilon = 1e-10;
    double lp_skip = log(epsilon);
    double lp_stay = log(p_stay);
    double lp_step = log(1.0 - exp(lp_skip) - exp(lp_stay));
    double lp_trim = log(0.01);

    // dp matrix
    int32_t n_rows = n_events + 1;
    int32_t n_cols = n_kmers + 1;
    int32_t n_bands = n_rows + n_cols;                                    

    // fill in remaining bands
    for (int32_t band_idx = 2; band_idx < n_bands; ++band_idx) {
        // Determine placement of this band according to Suzuki's adaptive algorithm
        // When both ll and ur are out-of-band (ob) we alternate movements
        // otherwise we decide based on scores
        float ll = BAND_ARRAY((band_idx - 1), 0);
        float ur = BAND_ARRAY((band_idx - 1),(bandwidth - 1));
        bool ll_ob = ll == -INFINITY;
        bool ur_ob = ur == -INFINITY;

        bool right = false;
        if (ll_ob && ur_ob) {
            right = band_idx % 2 == 1;
        } else {
            right = ll < ur; // Suzuki's rule
        }

        if (right) {
            band_lower_left[band_idx] =
                move_right(band_lower_left[band_idx - 1]);
        } else {
            band_lower_left[band_idx] =
                move_down(band_lower_left[band_idx - 1]);
        }
        // If the trim state is within the band, fill it in here
        int trim_offset = band_kmer_to_offset(band_idx, -1);
        if (is_offset_valid(trim_offset)) {
            int32_t event_idx = event_at_offset(band_idx, trim_offset);
            if (event_idx >= 0 && event_idx < n_events) {
                BAND_ARRAY(band_idx,trim_offset) = lp_trim * (event_idx + 1);
                TRACE_ARRAY(band_idx,trim_offset) = FROM_U;
            } else {
                BAND_ARRAY(band_idx,trim_offset) = -INFINITY;
            }
        }

        // Get the offsets for the first and last event and kmer
        // We restrict the inner loop to only these values
        int kmer_min_offset = band_kmer_to_offset(band_idx, 0);
        int kmer_max_offset = band_kmer_to_offset(band_idx, n_kmers);
        int event_min_offset = band_event_to_offset(band_idx, n_events - 1);
        int event_max_offset = band_event_to_offset(band_idx, -1);

        int min_offset = MAX(kmer_min_offset, event_min_offset);
        min_offset = MAX(min_offset, 0);

        int max_offset = MIN(kmer_max_offset, event_max_offset);
        max_offset = MIN(max_offset, bandwidth);

#ifdef DYNAMIC_PARALLELISM
        if(max_offset-min_offset > DYNAMIC_THRESH){
            int32_t BLOCK_LEN = DYNAMIC_BLOCK_LEN;
            int grid=((max_offset-min_offset + BLOCK_LEN - 1) / BLOCK_LEN);
            int block=(BLOCK_LEN);
            adaptive_align_inner<<<grid, block>>>(bands,trace,band_lower_left,min_offset, max_offset, band_idx, kmer_ranks,scaling, models, events,p_stay);        
            cudaDeviceSynchronize();
            cudaError_t code = cudaGetLastError();
            if (code != cudaSuccess) {
                printf("grid %d, block %d\n",grid,block);
                printf("Cuda error: %s \n in file : %s line number : %lu\n",
                        cudaGetErrorString(code), __FILE__, __LINE__);
                return ; //todo : generalise
            }
        }
        else{
#endif  //dynamic paralellism 

        for (int offset = min_offset; offset < max_offset; ++offset) {
            int event_idx = event_at_offset(band_idx, offset);
            int kmer_idx = kmer_at_offset(band_idx, offset);

            int32_t kmer_rank = kmer_ranks[kmer_idx];

            int offset_up = band_event_to_offset(band_idx - 1, event_idx - 1);
            int offset_left = band_kmer_to_offset(band_idx - 1, kmer_idx - 1);
            int offset_diag = band_kmer_to_offset(band_idx - 2, kmer_idx - 1);

#ifdef DEBUG_ADAPTIVE
            // verify loop conditions
            assert(kmer_idx >= 0 && kmer_idx < n_kmers);
            assert(event_idx >= 0 && event_idx < n_events);
            assert(offset_diag ==
                   band_event_to_offset(band_idx - 2, event_idx - 1));
            assert(offset_up - offset_left == 1);
            assert(offset >= 0 && offset < bandwidth);
#endif //DEBUG_ADAPTIVE

            float up = is_offset_valid(offset_up)
                           ? BAND_ARRAY(band_idx - 1,offset_up)
                           : -INFINITY;
            float left = is_offset_valid(offset_left)
                             ? BAND_ARRAY(band_idx - 1,offset_left)
                             : -INFINITY;
            float diag = is_offset_valid(offset_diag)
                             ? BAND_ARRAY(band_idx - 2,offset_diag)
                             : -INFINITY;

            float lp_emission = log_probability_match_r9(
                scaling, models, events, event_idx, kmer_rank);
            //fprintf(stderr, "lp emiision : %f , event idx %d, kmer rank %d\n", lp_emission,event_idx,kmer_rank);
            float score_d = diag + lp_step + lp_emission;
            float score_u = up + lp_stay + lp_emission;
            float score_l = left + lp_skip;

            float max_score = score_d;
            uint8_t from = FROM_D;

            max_score = score_u > max_score ? score_u : max_score;
            from = max_score == score_u ? FROM_U : from;
            max_score = score_l > max_score ? score_l : max_score;
            from = max_score == score_l ? FROM_L : from;

#ifdef DEBUG_ADAPTIVE
            fprintf(stderr,
                    "[adafill] offset-up: %d offset-diag: %d offset-left: %d\n",
                    offset_up, offset_diag, offset_left);
            fprintf(stderr, "[adafill] up: %.2lf diag: %.2lf left: %.2lf\n", up,
                    diag, left);
            fprintf(stderr,
                    "[adafill] bi: %d o: %d e: %d k: %d s: %.2lf f: %d emit: "
                    "%.2lf\n",
                    band_idx, offset, event_idx, kmer_idx, max_score, from,
                    lp_emission);
#endif //DEBUG_ADAPTIVE
            BAND_ARRAY(band_idx,offset) = max_score;
            TRACE_ARRAY(band_idx,offset) = from;
            //fills += 1;
        }

#ifdef DYNAMIC_PARALLELISM
        }
#endif  //dynamic paralellism     

    }

}

__forceinline__ __device__ int32_t align_single_post(AlignedPair* out_2, char* sequence,
    int32_t sequence_len, event_t* events,
    int32_t n_event, model_t* models,
    scalings_t scaling,int32_t* kmer_ranks,float *bands,uint8_t *trace, EventKmerPair* band_lower_left) {

    //fprintf(stderr, "%s\n", sequence);
    //fprintf(stderr, "Scaling %f %f", scaling.scale, scaling.shift);

    //size_t strand_idx = 0;
    //size_t k = 6;

    // size_t n_events = events[strand_idx].n;
    int32_t n_events = n_event;
    int32_t n_kmers = sequence_len - KMER_SIZE + 1;
    //fprintf(stderr,"n_kmers : %d\n",n_kmers);
    // backtrack markers
    //const uint8_t FROM_D = 0;
    //const uint8_t FROM_U = 1;
    //const uint8_t FROM_L = 2;

    // qc
    //double min_average_log_emission = -5.0;
    //int max_gap_threshold = 50;

    // banding
    //int bandwidth = ALN_BANDWIDTH;
    //half_bandwidth = bandwidth / 2;

    // transition penalties
    float events_per_kmer = (float)n_events / n_kmers;
    float p_stay = 1 - (1 / (events_per_kmer + 1));

    // setting a tiny skip penalty helps keep the true alignment within the adaptive band
    // this was empirically determined
    //double epsilon = 1e-10;
    double lp_skip = log(epsilon);
    double lp_stay = log(p_stay);
    double lp_step = log(1.0 - exp(lp_skip) - exp(lp_stay));
    double lp_trim = log(0.01);

    // dp matrix
    int32_t n_rows = n_events + 1;
    int32_t n_cols = n_kmers + 1;
    int32_t n_bands = n_rows + n_cols;
    //
    // Backtrack to compute alignment
    //
    double sum_emission = 0;
    double n_aligned_events = 0;

    //>>>>>>>>>>>>>> New replacement begin
    // std::vector<AlignedPair> out;

    int outIndex = 0;
    //<<<<<<<<<<<<<<<<New Replacement over

    float max_score = -INFINITY;
    int curr_event_idx = 0;
    int curr_kmer_idx = n_kmers - 1;

    // Find best score between an event and the last k-mer. after trimming the remaining evnets
    for (int32_t event_idx = 0; event_idx < n_events; ++event_idx) {
        int band_idx = event_kmer_to_band(event_idx, curr_kmer_idx);

        //>>>>>>>New  replacement begin
        /*assert(band_idx < bands.size());*/

        assert(band_idx < n_bands);

        //<<<<<<<<New Replacement over
        int offset = band_event_to_offset(band_idx, event_idx);
        if (is_offset_valid(offset)) {
            float s =
                BAND_ARRAY(band_idx,offset) + (n_events - event_idx) * lp_trim;
            if (s > max_score) {
                max_score = s;
                curr_event_idx = event_idx;
            }
        }
    }

#ifdef DEBUG_ADAPTIVE
    fprintf(stderr, "[adaback] ei: %d ki: %d s: %.2f\n", curr_event_idx,
            curr_kmer_idx, max_score);
#endif

    int curr_gap = 0;
    int max_gap = 0;
    while (curr_kmer_idx >= 0 && curr_event_idx >= 0) {
        // emit alignment
        //>>>>>>>New Repalcement begin
        assert(outIndex < n_events * 2);
        out_2[outIndex].ref_pos = curr_kmer_idx;
        out_2[outIndex].read_pos = curr_event_idx;
        outIndex++;
        // out.push_back({curr_kmer_idx, curr_event_idx});
        //<<<<<<<<<New Replacement over

#ifdef DEBUG_ADAPTIVE
        fprintf(stderr, "[adaback] ei: %d ki: %d\n", curr_event_idx,
                curr_kmer_idx);
#endif
        // qc stats
        //>>>>>>>>>>>>>>New Replacement begin
        char* substring = &sequence[curr_kmer_idx];
        int32_t kmer_rank = get_kmer_rank(substring, KMER_SIZE);
        //<<<<<<<<<<<<<New Replacement over
        float tempLogProb = log_probability_match_r9(
            scaling, models, events, curr_event_idx, kmer_rank);

        sum_emission += tempLogProb;
        //fprintf(stderr, "lp_emission %f \n", tempLogProb);
        //fprintf(stderr,"lp_emission %f, sum_emission %f, n_aligned_events %d\n",tempLogProb,sum_emission,outIndex);

        n_aligned_events += 1;

        int band_idx = event_kmer_to_band(curr_event_idx, curr_kmer_idx);
        int offset = band_event_to_offset(band_idx, curr_event_idx);
        assert(band_kmer_to_offset(band_idx, curr_kmer_idx) == offset);

        uint8_t from = TRACE_ARRAY(band_idx,offset);
        if (from == FROM_D) {
            curr_kmer_idx -= 1;
            curr_event_idx -= 1;
            curr_gap = 0;
        } else if (from == FROM_U) {
            curr_event_idx -= 1;
            curr_gap = 0;
        } else {
            curr_kmer_idx -= 1;
            curr_gap += 1;
            max_gap = MAX(curr_gap, max_gap);
        }
    }

    //>>>>>>>>New replacement begin
    // std::reverse(out.begin(), out.end());
    int c;
    int end = outIndex - 1;
    for (c = 0; c < outIndex / 2; c++) {
        int ref_pos_temp = out_2[c].ref_pos;
        int read_pos_temp = out_2[c].read_pos;
        out_2[c].ref_pos = out_2[end].ref_pos;
        out_2[c].read_pos = out_2[end].read_pos;
        out_2[end].ref_pos = ref_pos_temp;
        out_2[end].read_pos = read_pos_temp;
        end--;
    }

    // if(outIndex>1){
    //   AlignedPair temp={out_2[0].ref_pos,out[0].read_pos};
    //   int i;
    //   for(i=0;i<outIndex-1;i++){
    //     out_2[i]={out_2[outIndex-1-i].ref_pos,out[outIndex-1-i].read_pos};
    //   }
    //   out[outIndex-1]={temp.ref_pos,temp.read_pos};
    // }
    //<<<<<<<<<New replacement over

    // QC results
    double avg_log_emission = sum_emission / n_aligned_events;
    //fprintf(stderr,"sum_emission %f, n_aligned_events %f, avg_log_emission %f\n",sum_emission,n_aligned_events,avg_log_emission);
    //>>>>>>>>>>>>>New replacement begin
    bool spanned = out_2[0].ref_pos == 0 &&
                   out_2[outIndex - 1].ref_pos == int(n_kmers - 1);
    // bool spanned = out.front().ref_pos == 0 && out.back().ref_pos == n_kmers - 1;
    //<<<<<<<<<<<<<New replacement over
    //bool failed = false;
    if (avg_log_emission < min_average_log_emission || !spanned ||
        max_gap > max_gap_threshold) {
        //failed = true;
        //>>>>>>>>>>>>>New replacement begin
        outIndex = 0;
        // out.clear();
        //free(out_2);
        //out_2 = NULL;
        //<<<<<<<<<<<<<New replacement over
    }

    // free(kmer_ranks);
    // for (size_t i = 0; i < n_bands; i++) {
    //     free(bands[i]);
    //     free(trace[i]);
    // }
    // free(bands);
    // free(trace);
    // free(band_lower_left);
    //fprintf(stderr, "ada\t%s\t%s\t%.2lf\t%zu\t%.2lf\t%d\t%d\t%d\n", read.read_name.substr(0, 6).c_str(), failed ? "FAILED" : "OK", events_per_kmer, sequence.size(), avg_log_emission, curr_event_idx, max_gap, fills);
    //outSize=outIndex;
    //if(outIndex>500000)fprintf(stderr, "Max outSize %d\n", outIndex);
    return outIndex;
}


__global__ void align_kernel_pre(AlignedPair* event_align_pairs,
    int32_t* n_event_align_pairs, char* read,
    int32_t* read_len, int32_t* read_ptr,
    event_t* event_table, int32_t* n_events,
    int32_t* event_ptr, model_t* model,
    scalings_t* scalings, int32_t n_bam_rec,int32_t* kmer_ranks,float *bands,uint8_t *trace, EventKmerPair* band_lower_left) {

    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < n_bam_rec) {
        AlignedPair* out_2 = &event_align_pairs[event_ptr[i] * 2];
        char* sequence = &read[read_ptr[i]];
        int32_t sequence_len = read_len[i];
        event_t* events = &event_table[event_ptr[i]];
        int32_t n_event = n_events[i];
        scalings_t scaling = scalings[i];
        int32_t* kmer_rank = &kmer_ranks[read_ptr[i]];
        float *band = &bands[(read_ptr[i]+event_ptr[i])*ALN_BANDWIDTH];
        uint8_t *trace1 = &trace[(read_ptr[i]+event_ptr[i])*ALN_BANDWIDTH];
        EventKmerPair* band_lower_left1 = &band_lower_left[read_ptr[i]+event_ptr[i]];;

        align_single_pre(out_2, sequence, sequence_len,
                            events, n_event, model, scaling,kmer_rank,band,trace1,band_lower_left1);
    }
}


__global__ void 
//__launch_bounds__(MY_KERNEL_MAX_THREADS, MY_KERNEL_MIN_BLOCKS)
align_kernel_core(AlignedPair* event_align_pairs,
    int32_t* n_event_align_pairs, char* read,
    int32_t* read_len, int32_t* read_ptr,
    event_t* event_table, int32_t* n_events,
    int32_t* event_ptr, model_t* model,
    scalings_t* scalings, int32_t n_bam_rec,int32_t* kmer_ranks,float *bands,uint8_t *trace, EventKmerPair* band_lower_left) {

#ifndef WARP_HACK        
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < n_bam_rec) {
#else
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    int i = tid/32;    
    if (i < n_bam_rec && tid%32==0) {
#endif        
        AlignedPair* out_2 = &event_align_pairs[event_ptr[i] * 2];
        char* sequence = &read[read_ptr[i]];
        int32_t sequence_len = read_len[i];
        event_t* events = &event_table[event_ptr[i]];
        int32_t n_event = n_events[i];
        scalings_t scaling = scalings[i];
        int32_t* kmer_rank = &kmer_ranks[read_ptr[i]];
        float *band = &bands[(read_ptr[i]+event_ptr[i])*ALN_BANDWIDTH];
        uint8_t *trace1 = &trace[(read_ptr[i]+event_ptr[i])*ALN_BANDWIDTH];
        EventKmerPair* band_lower_left1 = &band_lower_left[read_ptr[i]+event_ptr[i]];;

        align_single(out_2, sequence, sequence_len,
                            events, n_event, model, scaling,kmer_rank,band,trace1,band_lower_left1);
    }
}


__global__ void align_kernel_post(AlignedPair* event_align_pairs,
    int32_t* n_event_align_pairs, char* read,
    int32_t* read_len, int32_t* read_ptr,
    event_t* event_table, int32_t* n_events,
    int32_t* event_ptr, model_t* model,
    scalings_t* scalings, int32_t n_bam_rec,int32_t* kmer_ranks,float *bands,uint8_t *trace, EventKmerPair* band_lower_left) {

    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < n_bam_rec) {
        AlignedPair* out_2 = &event_align_pairs[event_ptr[i] * 2];
        char* sequence = &read[read_ptr[i]];
        int32_t sequence_len = read_len[i];
        event_t* events = &event_table[event_ptr[i]];
        int32_t n_event = n_events[i];
        scalings_t scaling = scalings[i];
        int32_t* kmer_rank = &kmer_ranks[read_ptr[i]];
        float *band = &bands[(read_ptr[i]+event_ptr[i])*ALN_BANDWIDTH];
        uint8_t *trace1 = &trace[(read_ptr[i]+event_ptr[i])*ALN_BANDWIDTH];
        EventKmerPair* band_lower_left1 = &band_lower_left[read_ptr[i]+event_ptr[i]];;

        n_event_align_pairs[i] = align_single_post(out_2, sequence, sequence_len,
                            events, n_event, model, scaling,kmer_rank,band,trace1,band_lower_left1);
    }
}



__global__ void 
//__launch_bounds__(MY_KERNEL_MAX_THREADS, MY_KERNEL_MIN_BLOCKS)
align_kernel_core_2d(AlignedPair* event_align_pairs,
    int32_t* n_event_align_pairs, char* read,
    int32_t* read_len, int32_t* read_ptr,
    event_t* event_table, int32_t* n_events1,
    int32_t* event_ptr, model_t* models,
    scalings_t* scalings, int32_t n_bam_rec,int32_t* kmer_rank,float *band,uint8_t *traces, EventKmerPair* band_lower_lefts) {
   
    int i = blockDim.y * blockIdx.y + threadIdx.y;
    int tid=blockIdx.x*blockDim.x+threadIdx.x;

    if (i < n_bam_rec) {   

        AlignedPair* out_2 = &event_align_pairs[event_ptr[i] * 2];
        char* sequence = &read[read_ptr[i]];
        int32_t sequence_len = read_len[i];
        event_t* events = &event_table[event_ptr[i]];
        int32_t n_event = n_events1[i];
        scalings_t scaling = scalings[i];
        int32_t* kmer_ranks = &kmer_rank[read_ptr[i]];
        float *bands = &band[(read_ptr[i]+event_ptr[i])*ALN_BANDWIDTH];
        uint8_t *trace = &traces[(read_ptr[i]+event_ptr[i])*ALN_BANDWIDTH];
        EventKmerPair* band_lower_left = &band_lower_lefts[read_ptr[i]+event_ptr[i]];;

        // size_t n_events = events[strand_idx].n;
        int32_t n_events = n_event;
        int32_t n_kmers = sequence_len - KMER_SIZE + 1;
        //fprintf(stderr,"n_kmers : %d\n",n_kmers);

        // transition penalties
        float events_per_kmer = (float)n_events / n_kmers;
        float p_stay = 1 - (1 / (events_per_kmer + 1));

        // setting a tiny skip penalty helps keep the true alignment within the adaptive band
        // this was empirically determined
        //double epsilon = 1e-10;
        double lp_skip = log(epsilon);
        double lp_stay = log(p_stay);
        double lp_step = log(1.0 - exp(lp_skip) - exp(lp_stay));
        double lp_trim = log(0.01);

        // dp matrix
        int32_t n_rows = n_events + 1;
        int32_t n_cols = n_kmers + 1;
        int32_t n_bands = n_rows + n_cols;                                    

        // fill in remaining bands
        for (int32_t band_idx = 2; band_idx < n_bands; ++band_idx) {
            // Determine placement of this band according to Suzuki's adaptive algorithm
            // When both ll and ur are out-of-band (ob) we alternate movements
            // otherwise we decide based on scores
            float ll = BAND_ARRAY((band_idx - 1), 0);
            float ur = BAND_ARRAY((band_idx - 1),(bandwidth - 1));
            bool ll_ob = ll == -INFINITY;
            bool ur_ob = ur == -INFINITY;

            bool right = false;
            if (ll_ob && ur_ob) {
                right = band_idx % 2 == 1;
            } else {
                right = ll < ur; // Suzuki's rule
            }

            if (right) {
                band_lower_left[band_idx] =
                    move_right(band_lower_left[band_idx - 1]);
            } else {
                band_lower_left[band_idx] =
                    move_down(band_lower_left[band_idx - 1]);
            }
            // If the trim state is within the band, fill it in here
            int trim_offset = band_kmer_to_offset(band_idx, -1);
            if (is_offset_valid(trim_offset)) {
                int32_t event_idx = event_at_offset(band_idx, trim_offset);
                if (event_idx >= 0 && event_idx < n_events) {
                    BAND_ARRAY(band_idx,trim_offset) = lp_trim * (event_idx + 1);
                    TRACE_ARRAY(band_idx,trim_offset) = FROM_U;
                } else {
                    BAND_ARRAY(band_idx,trim_offset) = -INFINITY;
                }
            }

            // Get the offsets for the first and last event and kmer
            // We restrict the inner loop to only these values
            int kmer_min_offset = band_kmer_to_offset(band_idx, 0);
            int kmer_max_offset = band_kmer_to_offset(band_idx, n_kmers);
            int event_min_offset = band_event_to_offset(band_idx, n_events - 1);
            int event_max_offset = band_event_to_offset(band_idx, -1);

            int min_offset = MAX(kmer_min_offset, event_min_offset);
            min_offset = MAX(min_offset, 0);

            int max_offset = MIN(kmer_max_offset, event_max_offset);
            max_offset = MIN(max_offset, bandwidth);

            __syncthreads();    
   
            if(tid < (max_offset-min_offset)) {
                int offset=min_offset+tid;

                int event_idx = event_at_offset(band_idx, offset);
                int kmer_idx = kmer_at_offset(band_idx, offset);

                int32_t kmer_rank = kmer_ranks[kmer_idx];

                int offset_up = band_event_to_offset(band_idx - 1, event_idx - 1);
                int offset_left = band_kmer_to_offset(band_idx - 1, kmer_idx - 1);
                int offset_diag = band_kmer_to_offset(band_idx - 2, kmer_idx - 1);

    #ifdef DEBUG_ADAPTIVE
                // verify loop conditions
                assert(kmer_idx >= 0 && kmer_idx < n_kmers);
                assert(event_idx >= 0 && event_idx < n_events);
                assert(offset_diag ==
                       band_event_to_offset(band_idx - 2, event_idx - 1));
                assert(offset_up - offset_left == 1);
                assert(offset >= 0 && offset < bandwidth);
    #endif //DEBUG_ADAPTIVE

                float up = is_offset_valid(offset_up)
                               ? BAND_ARRAY(band_idx - 1,offset_up)
                               : -INFINITY;
                float left = is_offset_valid(offset_left)
                                 ? BAND_ARRAY(band_idx - 1,offset_left)
                                 : -INFINITY;
                float diag = is_offset_valid(offset_diag)
                                 ? BAND_ARRAY(band_idx - 2,offset_diag)
                                 : -INFINITY;

                float lp_emission = log_probability_match_r9(
                    scaling, models, events, event_idx, kmer_rank);
                //fprintf(stderr, "lp emiision : %f , event idx %d, kmer rank %d\n", lp_emission,event_idx,kmer_rank);
                float score_d = diag + lp_step + lp_emission;
                float score_u = up + lp_stay + lp_emission;
                float score_l = left + lp_skip;

                float max_score = score_d;
                uint8_t from = FROM_D;

                max_score = score_u > max_score ? score_u : max_score;
                from = max_score == score_u ? FROM_U : from;
                max_score = score_l > max_score ? score_l : max_score;
                from = max_score == score_l ? FROM_L : from;

    #ifdef DEBUG_ADAPTIVE
                fprintf(stderr,
                        "[adafill] offset-up: %d offset-diag: %d offset-left: %d\n",
                        offset_up, offset_diag, offset_left);
                fprintf(stderr, "[adafill] up: %.2lf diag: %.2lf left: %.2lf\n", up,
                        diag, left);
                fprintf(stderr,
                        "[adafill] bi: %d o: %d e: %d k: %d s: %.2lf f: %d emit: "
                        "%.2lf\n",
                        band_idx, offset, event_idx, kmer_idx, max_score, from,
                        lp_emission);
    #endif //DEBUG_ADAPTIVE
                BAND_ARRAY(band_idx,offset) = max_score;
                TRACE_ARRAY(band_idx,offset) = from;
                //fills += 1;
            }

            __syncthreads();  
   
        }

    }

}



#endif