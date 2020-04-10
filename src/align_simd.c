#include "f5c.h"
#include <emmintrin.h>
#include <assert.h>
#include <sys/time.h>
//#define DEBUG_ESTIMATED_SCALING 1
//#define DEBUG_RECALIB_SCALING 1
//#define DEBUG_ADAPTIVE 1

//Code was adapted from Nanopolish : nanopolish_raw_loader.cpp

//todo : can make more efficient using bit encoding
static inline uint32_t get_rank(char base) {
    if (base == 'A') { //todo: do we neeed simple alpha?
        return 0;
    } else if (base == 'C') {
        return 1;
    } else if (base == 'G') {
        return 2;
    } else if (base == 'T') {
        return 3;
    } else {
        WARNING("A None ACGT base found : %c", base);
        return 0;
    }
}

// return the lexicographic rank of the kmer amongst all strings of
// length k for this alphabet
static inline uint32_t get_kmer_rank(const char* str, uint32_t k) {
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
static inline void kmer_cpy(char* dest, char* src, uint32_t k) {
    uint32_t i = 0;
    for (i = 0; i < k; i++) {
        dest[i] = src[i];
    }
    dest[i] = '\0';
}

static inline float log_normal_pdf(float x, float gp_mean, float gp_stdv,
                                   float gp_log_stdv) {
    /*INCOMPLETE*/
    float log_inv_sqrt_2pi = -0.918938f; // Natural logarithm
    float a = (x - gp_mean) / gp_stdv;
    return log_inv_sqrt_2pi - gp_log_stdv + (-0.5f * a * a);
    // return 1;
}

static inline float log_probability_match_r9(scalings_t scaling,
                                             model_t* models,
                                             event_table events, int32_t event_idx,
                                             int32_t kmer_rank, uint8_t strand,
                                             float sample_rate) {
    // event level mean, scaled with the drift value
    strand = 0;
    assert(kmer_rank < 4096);
    //float level = read.get_drift_scaled_level(event_idx, strand);

    //float time =
    //    (events.event[event_idx].start - events.event[0].start) / sample_rate;
    float unscaledLevel = events.event[event_idx].mean;
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
    #ifdef CACHED_LOG
        float gp_log_stdv = models[kmer_rank].level_log_stdv;
    #else
        float gp_log_stdv =
        log(models[kmer_rank].level_stdv); // scaling.log_var = log(1)=0;
    #endif

    float lp = log_normal_pdf(scaledLevel, gp_mean, gp_stdv, gp_log_stdv);
    return lp;
}

//SSE2 helper function to compare two floating point vectors. Returns 1 in positions where the value is the same and 0 if not.
__m128i compare_from_vec(__m128 vec1, __m128 vec2){
    return _mm_add_epi32((__m128i)_mm_cmpgt_ps(vec1,vec2),_mm_set1_epi32(1));
}

//Helper function to print a single int vector
void print_int_vec(__m128i vec,const char *vec_name){
    int32_t * arr = (int32_t *)malloc(4 * sizeof(int32_t));;
    arr[0] = _mm_cvtsi128_si32(vec);
    arr[1] = _mm_cvtsi128_si32(_mm_shuffle_epi32(vec,1));
    arr[2] = _mm_cvtsi128_si32(_mm_shuffle_epi32(vec,2));
    arr[3] = _mm_cvtsi128_si32 (_mm_shuffle_epi32(vec,3));

    fprintf(stderr, "%s: (%d,%d,%d,%d)\n",vec_name,arr[0],arr[1],arr[2],arr[3]);

    free(arr);
}



//Helper function to print a single band vector and a from vector
void print_float_vec(__m128 vec, const char *vec_name){
    float * arr = (float *)malloc(4 * sizeof(float));

    _mm_store_ps(arr,vec);

    fprintf(stderr, "%s: (%.2f,%.2f,%.2f,%.2f)\n",vec_name,arr[3],arr[2],arr[1],arr[0]);

    free(arr);
}

#define event_kmer_to_band(ei, ki) (ei + 1) + (ki + 1)
#define band_event_to_offset(bi, ei) band_lower_left[bi].event_idx - (ei)
#define band_kmer_to_offset(bi, ki) (ki) - band_lower_left[bi].kmer_idx
#define is_offset_valid(offset) (offset) >= 0 && (offset) < bandwidth
#define event_at_offset(bi, offset) band_lower_left[(bi)].event_idx - (offset)
#define kmer_at_offset(bi, offset) band_lower_left[(bi)].kmer_idx + (offset)

#define move_down(curr_band) \
    { curr_band.event_idx + 1, curr_band.kmer_idx }
#define move_right(curr_band) \
    { curr_band.event_idx, curr_band.kmer_idx + 1 }

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))


#ifdef ALIGN_2D_ARRAY
    #define BAND_ARRAY(r, c) ( bands[(r)][(c)] )
    #define TRACE_ARRAY(r, c) ( trace[(r)][(c)] )
#else
    #define BAND_ARRAY(r, c) ( bands[((r)*(ALN_BANDWIDTH)+(c))] )
    #define TRACE_ARRAY(r, c) ( trace[((r)*(ALN_BANDWIDTH)+(c))] )
#endif


int32_t align_simd(AlignedPair* out_2, char* sequence, int32_t sequence_len,
              event_table events, model_t* models, scalings_t scaling,
              float sample_rate) {
    //fprintf(stderr, "%s\n", sequence);
    //fprintf(stderr, "Scaling %f %f", scaling.scale, scaling.shift);

    //INFO("%s", "Running in simd mode");

    int32_t strand_idx = 0;
    int32_t k = 6;

    // int32_t n_events = events[strand_idx].n;
    int32_t n_events = events.n;
    int32_t n_kmers = sequence_len - k + 1;
    //fprintf(stderr,"n_kmers : %d\n",n_kmers);
    // backtrack markers
    const int32_t FROM_D = 0;
    const int32_t FROM_U = 1;
    const int32_t FROM_L = 2;

    // qc
    float min_average_log_emission = -5.0;
    int max_gap_threshold = 50;

    // banding
    int32_t bandwidth = ALN_BANDWIDTH;
    int32_t half_bandwidth = ALN_BANDWIDTH / 2;

    //SSE2 number of vectors in a band
    // int32_t bandwidth_vec = sse2_convert_size(bandwidth);

    // transition penalties
    float events_per_kmer = (float)n_events / n_kmers;
    float p_stay = 1 - (1 / (events_per_kmer + 1));

    // setting a tiny skip penalty helps keep the true alignment within the adaptive band
    // this was empirically determined
    float epsilon = 1e-10;
    float lp_skip = log(epsilon);
    float lp_stay = log(p_stay);
    float lp_step = log(1.0 - exp(lp_skip) - exp(lp_stay));
    float lp_trim = log(0.01);

    //SSE2 constant vectors
    __m128 lp_skip_vec = _mm_set1_ps(lp_skip);
    __m128 lp_stay_vec = _mm_set1_ps(lp_stay);
    __m128 lp_step_vec = _mm_set1_ps(lp_step);

    // dp matrix
    int32_t n_rows = n_events + 1;
    int32_t n_cols = n_kmers + 1;
    int32_t n_bands = n_rows + n_cols;

    // Initialize

    // Precompute k-mer ranks to avoid doing this in the inner loop
    int32_t* kmer_ranks = (int32_t*)malloc(sizeof(int32_t) * n_kmers);
    MALLOC_CHK(kmer_ranks);

    for (int32_t i = 0; i < n_kmers; ++i) {
        //>>>   >>>>>> New replacement begin
        char* substring = &sequence[i];
        kmer_ranks[i] = get_kmer_rank(substring, k);
        //<<<<<<<<< New replacement over
    }

#ifdef ALIGN_2D_ARRAY
    float** bands = (float**)malloc(sizeof(float*) * n_bands);
    MALLOC_CHK(bands);
    int32_t** trace = (int32_t**)malloc(sizeof(int32_t*) * n_bands);
    MALLOC_CHK(trace);
#else
    float* bands = (float*)malloc(sizeof(float) * n_bands * bandwidth);
    MALLOC_CHK(bands);
    int32_t* trace = (int32_t*)malloc(sizeof(int32_t) * n_bands * bandwidth);
    MALLOC_CHK(trace);

#endif

    //Initialize default values
    for (int32_t i = 0; i < n_bands; i++) {
    #ifdef ALIGN_2D_ARRAY
        bands[i] = (float*)malloc(sizeof(float) * bandwidth);
        MALLOC_CHK(bands[i]);
        trace[i] = (int32_t*)malloc(sizeof(int32_t) * bandwidth);
        MALLOC_CHK(trace[i]);
    #endif

        for (int32_t j = 0; j < bandwidth; j++) {
            BAND_ARRAY(i,j) = -INFINITY;
            TRACE_ARRAY(i,j) = 0;
        }
    }

    // Keep track of the event/kmer index for the lower left corner of the band
    // these indices are updated at every iteration to perform the adaptive banding
    // Only the first two  have their coordinates initialized, the rest are computed adaptively

    struct EventKmerPair {
        int32_t event_idx;
        int32_t kmer_idx;
    };
    //>>>>>>>>>>>>>>>>>New Replacement Begin
    struct EventKmerPair* band_lower_left =
        (struct EventKmerPair*)malloc(sizeof(struct EventKmerPair) * n_bands);
    MALLOC_CHK(band_lower_left);
    //std::vector<EventKmerPair> band_lower_left(n_);
    //<<<<<<<<<<<<<<<<<New Replacement over

    // initialize range of first two
    band_lower_left[0].event_idx = half_bandwidth - 1;
    band_lower_left[0].kmer_idx = -1 - half_bandwidth;
    band_lower_left[1] = move_down(band_lower_left[0]);

    int32_t start_cell_offset = band_kmer_to_offset(0, -1);
    // assert(is_offset_valid(start_cell_offset));
    // assert(band_event_to_offset(0, -1) == start_cell_offset);
    BAND_ARRAY(0,start_cell_offset) = 0.0f;

    // band 1: first event is trimmed
    int32_t first_trim_offset = band_event_to_offset(1, 0);
    // assert(kmer_at_offset(1, first_trim_offset) == -1);
    // assert(is_offset_valid(first_trim_offset));
    BAND_ARRAY(1,first_trim_offset) = lp_trim;
    TRACE_ARRAY(1,first_trim_offset) = FROM_U;

    //Declare/initialise variables needed for the loops
    int32_t event_idx,kmer_idx,kmer_rank,offset_up,offset_left,offset_diag;
    float lp_emission;
    int32_t fills = 0;
    float *up_arr = (float *)malloc(sizeof(float) * 4);
    MALLOC_CHK(up_arr);
    float *left_arr = (float *)malloc(sizeof(float) * 4);
    MALLOC_CHK(left_arr);
    float *diag_arr = (float *)malloc(sizeof(float) * 4);
    MALLOC_CHK(diag_arr);
    float *lp_emission_arr = (float *)malloc(sizeof(float) * 4);
    MALLOC_CHK(lp_emission_arr);
    int32_t *from_arr = (int32_t *)malloc(sizeof(int32_t) * 4);
    MALLOC_CHK(from_arr);
    float * band_arr = (float *)malloc(sizeof(float) * 4);
    MALLOC_CHK(band_arr);

#ifdef DEBUG_ADAPTIVE
    fprintf(stderr, "[trim] bi: %d o: %d e: %d k: %d s: %.2lf\n", 1,
            first_trim_offset, 0, -1, bands[1][first_trim_offset]);
#endif

    // fill in remaining bands
    for (int32_t band_idx = 2; band_idx < n_bands; ++band_idx) {
        // Determine placement of this band according to Suzuki's adaptive algorithm
        // When both ll and ur are out-of-band (ob) we alternate movements
        // otherwise we decide based on scores
        float ll = BAND_ARRAY(band_idx - 1,0);
        float ur = BAND_ARRAY(band_idx - 1,bandwidth - 1);

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
        int32_t trim_offset = band_kmer_to_offset(band_idx, -1);
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
        int32_t kmer_min_offset = band_kmer_to_offset(band_idx, 0);
        int32_t kmer_max_offset = band_kmer_to_offset(band_idx, n_kmers);
        int32_t event_min_offset = band_event_to_offset(band_idx, n_events - 1);
        int32_t event_max_offset = band_event_to_offset(band_idx, -1);

        int32_t min_offset = MAX(kmer_min_offset, event_min_offset);
        min_offset = MAX(min_offset, 0);

        int32_t max_offset = MIN(kmer_max_offset, event_max_offset);
        max_offset = MIN(max_offset, bandwidth);

        //Inner loop: Parallelised with SSE2. Jump up 4 every time
        for (int32_t offset = min_offset; offset < max_offset; offset += 4) {
            if(offset + 4 >= max_offset){
                //If we don't have >= 4 cells left to fill, compute individually using a sequential loop
                for(int32_t seq_offset = offset; seq_offset < max_offset; seq_offset++){
                    event_idx = event_at_offset(band_idx, seq_offset);
                    kmer_idx = kmer_at_offset(band_idx, seq_offset);
                    kmer_rank = kmer_ranks[kmer_idx];
                    //simd_debug
                    // fprintf(stderr, "event idx %d, kmer_idx: %d, kmer rank %d, band_idx: %d, seq_offset: %d, bllevent: %d, bllkmer: %d\n",
                    // event_idx,kmer_idx,kmer_rank,band_idx,seq_offset,band_lower_left[band_idx].event_idx,band_lower_left[band_idx].kmer_idx);

                    //Offset of the up, left, and diagonal positions
                    offset_up = band_event_to_offset(band_idx - 1, event_idx - 1);
                    offset_left = band_kmer_to_offset(band_idx - 1, kmer_idx - 1);
                    offset_diag = band_kmer_to_offset(band_idx - 2, kmer_idx - 1);

#ifdef DEBUG_ADAPTIVE
                    // verify loop conditions
                    assert(kmer_idx >= 0 && kmer_idx < n_kmers);
                    assert(event_idx >= 0 && event_idx < n_events);
                    assert(offset_diag == band_event_to_offset(band_idx - 2, event_idx - 1));
                    assert(offset_up - offset_left == 1);
                    assert(seq_offset >= 0 && seq_offset < bandwidth);
#endif

                    float up = is_offset_valid(offset_up)
                            ? BAND_ARRAY(band_idx - 1,offset_up)
                            : -INFINITY;

                    float left = is_offset_valid(offset_left)
                                ? BAND_ARRAY(band_idx - 1,offset_left)
                                : -INFINITY;

                    float diag = is_offset_valid(offset_diag)
                                ? BAND_ARRAY(band_idx - 2,offset_diag)
                                : -INFINITY;

                    lp_emission = log_probability_match_r9(scaling, models, events, event_idx,
                                            kmer_rank, strand_idx, sample_rate);

                    //Compute score and from values for single entry
                    float score_d_s = diag + lp_step + lp_emission;
                    float score_u_s = up + lp_stay + lp_emission;
                    float score_l_s = left + lp_skip;

                    //A single max_score/from calculation
                    float max_score_single = score_d_s;
                    int32_t from_single = FROM_D;

                    max_score_single = score_u_s > max_score_single ? score_u_s : max_score_single;
                    from_single = max_score_single == score_u_s ? FROM_U : from_single;
                    max_score_single = score_l_s > max_score_single ? score_l_s : max_score_single;
                    from_single = max_score_single == score_l_s ? FROM_L : from_single;

                    //simd_debug
                    // fprintf(stderr, "offset: %d, score: %f, from : %d\n",seq_offset,max_score_single,from_single);

                    //Store in arrays
                    BAND_ARRAY(band_idx,seq_offset) = max_score_single;
                    TRACE_ARRAY(band_idx,seq_offset) = from_single;
                    fills += 1;
                }
            }else{
                //Compute using SIMD
                //Need to load values sequentially because the __m128i vectors don't overlap
                //Load 4 values corresponding to the left, up and diagonal bands into arrays
                for (int32_t vec_pos = 0; vec_pos < 4 ; ++vec_pos) {
                    //Index of the first element of the vector we are looking at
                    event_idx = event_at_offset(band_idx, offset + vec_pos);
                    kmer_idx = kmer_at_offset(band_idx, offset + vec_pos);
                    //simd_debug
                    // fprintf(stderr, "event idx %d, kmer_idx: %d, kmer rank %d, band_idx: %d, vec_pos: %d, offset: %d, bllevent: %d, bllkmer: %d\n", event_idx,kmer_idx,kmer_rank,
                    // band_idx,vec_pos,offset,band_lower_left[band_idx].event_idx,band_lower_left[band_idx].kmer_idx);
                    kmer_rank = kmer_ranks[kmer_idx];

                    //Offset of the up, left, and diagonal positions
                    offset_up = band_event_to_offset(band_idx - 1, event_idx - 1);
                    offset_left = band_kmer_to_offset(band_idx - 1, kmer_idx - 1);
                    offset_diag = band_kmer_to_offset(band_idx - 2, kmer_idx - 1);

#ifdef DEBUG_ADAPTIVE
                    // verify loop conditions
                    assert(kmer_idx >= 0 && kmer_idx < n_kmers);
                    assert(event_idx >= 0 && event_idx < n_events);
                    assert(offset_diag == band_event_to_offset(band_idx - 2, event_idx - 1));
                    assert(offset_up - offset_left == 1);
                    assert(offset >= 0 && offset < bandwidth);
#endif

                    up_arr[vec_pos] = is_offset_valid(offset_up)
                            ? BAND_ARRAY(band_idx - 1,offset_up)
                            : -INFINITY;

                    left_arr[vec_pos] = is_offset_valid(offset_left)
                                ? BAND_ARRAY(band_idx - 1,offset_left)
                                : -INFINITY;

                    diag_arr[vec_pos] = is_offset_valid(offset_diag)
                                ? BAND_ARRAY(band_idx - 2,offset_diag)
                                : -INFINITY;

                    lp_emission_arr[vec_pos] = log_probability_match_r9(scaling, models, events, event_idx,
                                            kmer_rank, strand_idx, sample_rate);
                }

                //convert data from the arrays to __m128
                __m128 up_vec = _mm_set_ps(up_arr[0],up_arr[1],up_arr[2],up_arr[3]);
                __m128 left_vec = _mm_set_ps(left_arr[0],left_arr[1],left_arr[2],left_arr[3]);
                __m128 diag_vec = _mm_set_ps(diag_arr[0],diag_arr[1],diag_arr[2],diag_arr[3]);
                __m128 lp_emission_vec = _mm_set_ps(lp_emission_arr[0],lp_emission_arr[1],lp_emission_arr[2],lp_emission_arr[3]);

                __m128 score_d = _mm_add_ps(diag_vec,_mm_add_ps(lp_step_vec,lp_emission_vec));
                __m128 score_u = _mm_add_ps(up_vec,_mm_add_ps(lp_stay_vec,lp_emission_vec));
                __m128 score_l = _mm_add_ps(left_vec,lp_skip_vec);

                // __m128 max_score = score_d;

                __m128 max_score = _mm_max_ps(score_l,_mm_max_ps(score_u,score_d));

                //These vectors have a 1 where the max_score corresponds to the direction, and 0 where it doesn't
                __m128i compare_up = compare_from_vec(max_score,score_u);
                __m128i compare_left = compare_from_vec(max_score,score_l);

                //FROM_D=0, FROM_U=1, FROM_L=2, so only need to add compare_up to 2 * compare_left
                __m128i from = _mm_add_epi32(compare_up,_mm_add_epi32(compare_left,compare_left));

                //simd_debug
                // print_int_vec(compare_up,"compare_up");
                // print_int_vec(compare_left,"compare_left");
                // print_int_vec(from,"from");
                // print_float_vec(score_u,"score_u");
                // print_float_vec(score_l,"score_l");
                // print_float_vec(max_score,"max_score");

                //Store results in BAND and TRACE array
                float *band_position = &(BAND_ARRAY(band_idx,offset));
                int32_t *trace_position = &(TRACE_ARRAY(band_idx,offset));

                //reverse array to store scores in correct order. need to use temp array because storer needs aligned boundary
                _mm_storer_ps(band_arr,max_score);
                band_position[0] = band_arr[3];
                band_position[1] = band_arr[2];
                band_position[2] = band_arr[1];
                band_position[3] = band_arr[0];

                trace_position[3] = _mm_cvtsi128_si32(from);
                trace_position[2] = _mm_cvtsi128_si32(_mm_shuffle_epi32(from,1));
                trace_position[1] = _mm_cvtsi128_si32(_mm_shuffle_epi32(from,2));
                trace_position[0] = _mm_cvtsi128_si32 (_mm_shuffle_epi32(from,3));

                //simd_debug
                //fprintf(stderr,"bands: %f %f %f %f\n",band_position[0],band_position[1],band_position[2],band_position[3]);
                // fprintf(stderr,"\nband array: %f %f %f %f\n",BAND_ARRAY(band_idx,offset),BAND_ARRAY(band_idx,offset+1),BAND_ARRAY(band_idx,offset+2),BAND_ARRAY(band_idx,offset+3));

                fills += 4;
            }
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
#endif
        }
    }
    // Backtrack to compute alignment
    //
    double sum_emission = 0;
    double n_aligned_events = 0;

    //>>>>>>>>>>>>>> New replacement begin
    // std::vector<AlignedPair> out;

    int32_t outIndex = 0;
    //<<<<<<<<<<<<<<<<New Replacement over

    float max_score = -INFINITY;
    int32_t curr_event_idx = 0;
    int32_t curr_kmer_idx = n_kmers - 1;

    // Find best score between an event and the last k-mer. after trimming the remaining evnets
    for (int32_t event_idx = 0; event_idx < n_events; ++event_idx) {
        int32_t band_idx = event_kmer_to_band(event_idx, curr_kmer_idx);

        //>>>>>>>New  replacement begin
        /*assert(band_idx < bands.size());*/
        assert((int32_t)band_idx < n_bands);
        //<<<<<<<<New Replacement over
        int32_t offset = band_event_to_offset(band_idx, event_idx);

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

    int32_t curr_gap = 0;
    int32_t max_gap = 0;
    while (curr_kmer_idx >= 0 && curr_event_idx >= 0) {
        // emit alignment
        //>>>>>>>New Repalcement begin
        assert(outIndex < (int32_t)(n_events * 2));
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
        int32_t kmer_rank = get_kmer_rank(substring, k);
        //<<<<<<<<<<<<<New Replacement over
        float tempLogProb = log_probability_match_r9(
            scaling, models, events, curr_event_idx, kmer_rank, 0, sample_rate);


        sum_emission += tempLogProb;
        // fprintf(stderr, "lp_emission %f \n", tempLogProb);

        n_aligned_events += 1;

        int32_t band_idx = event_kmer_to_band(curr_event_idx, curr_kmer_idx);
        int32_t offset = band_event_to_offset(band_idx, curr_event_idx);
        assert(band_kmer_to_offset(band_idx, curr_kmer_idx) == offset);

        int32_t from = TRACE_ARRAY(band_idx,offset);
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
    int32_t c;
    int32_t end = outIndex - 1;
    for (c = 0; c < outIndex / 2; c++) {
        int32_t ref_pos_temp = out_2[c].ref_pos;
        int32_t read_pos_temp = out_2[c].read_pos;
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

    //>>>>>>>>>>>>>New replacement begin
    bool spanned = out_2[0].ref_pos == 0 &&
                   out_2[outIndex - 1].ref_pos == int32_t(n_kmers - 1);
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

    free(kmer_ranks);
#ifdef ALIGN_2D_ARRAY
    for (int32_t i = 0; i < n_bands; i++) {
        free(bands[i]);
        free(trace[i]);
    }
#endif

    //free regular mallocs
    free(bands);
    free(trace);
    free(band_lower_left);

    //free sse2 mallocs
    free(up_arr);
    free(left_arr);
    free(diag_arr);
    free(lp_emission_arr);
    free(from_arr);
    free(band_arr);

    //fprintf(stderr, "ada\t%s\t%s\t%.2lf\t%zu\t%.2lf\t%d\t%d\t%d\n", read.read_name.substr(0, 6).c_str(), failed ? "FAILED" : "OK", events_per_kmer, sequence.size(), avg_log_emission, curr_event_idx, max_gap, fills);
    //outSize=outIndex;
    //if(outIndex>500000)fprintf(stderr, "Max outSize %d\n", outIndex);
    return outIndex;
}
