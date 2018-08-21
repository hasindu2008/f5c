#include "f5c.h"
#include <vector>
#include <algorithm>
#include <cassert>


#define DEBUG_PRINT_STATS

#define event_kmer_to_band(ei, ki) (ei + 1) + (ki + 1)
#define band_event_to_offset(bi, ei) band_lower_left[bi].event_idx - (ei)
#define band_kmer_to_offset(bi, ki) (ki) - band_lower_left[bi].kmer_idx
#define is_offset_valid(offset) (offset) >= 0 && (offset) < bandwidth
#define event_at_offset(bi, offset) band_lower_left[(bi)].event_idx - (offset)
#define kmer_at_offset(bi, offset) band_lower_left[(bi)].kmer_idx + (offset)
#define move_down(curr_band) { curr_band.event_idx + 1, curr_band.kmer_idx }
#define move_right(curr_band) { curr_band.event_idx, curr_band.kmer_idx + 1 }

float log_normal_pdf(float x, float gp_mean, float gp_stdv, float gp_log_stdv)
{
    /*INCOMPLETE*/
    float log_inv_sqrt_2pi=-0.39908993417f;// <<<<<<<<<<< Are we working on log base 10?
    float a = (x - mean) / stdv;
    return log_inv_sqrt_2pi - gp_log_stdv + (-0.5f * a * a);
    return 1;
}


float log_probability_match_r9(scalings_t scaling,
                                      model_t* models,
                                      event_table events,
                                      int event_idx,
                                      uint32_t kmer_rank,
                                      uint8_t strand)
{
    // event level mean, scaled with the drift value
    strand=0;

    //float level = read.get_drift_scaled_level(event_idx, strand);

    float time=events.event[event_idx].start- events.event[0].start;
    float unscaledLevel=events.event[event_idx].mean;
    float scaledLevel=unscaledLevel - time* scaling.shift;

    //GaussianParameters gp = read.get_scaled_gaussian_from_pore_model_state(pore_model, strand, kmer_rank);
    float gp_mean = scaling.scale * models[kmer_rank].level_mean + scaling.shift;
    //float gp_stdv = models[kmer_rank].level_stdv * scaling.var;
    float gp_stdv = 0;
    // float gp_log_stdv = models[kmer_rank].level_log_stdv + scaling.log_var;
    float gp_log_stdv = models[kmer_rank].level_stdv + scaling.log_var;

    float lp = log_normal_pdf(scaledLevel, gp_mean,gp_stdv,gp_log_stdv);
    return lp;
}

//todo : can make more efficient using bit encoding
static inline uint32_t rank(char base) {
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
static inline uint32_t kmer_rank_function(const char* str, uint32_t k) {
    //uint32_t p = 1;
    uint32_t r = 0;

    // from last base to first
    for (uint32_t i = 0; i < k; ++i) {
        //r += rank(str[k - i - 1]) * p;
        //p *= size();
        r += rank(str[k - i - 1]) << 2;
    }
    return r;
}

scalings_t estimate_scalings_using_mom(char* sequence, model_t* pore_model,
                                       event_table et) {
    scalings_t out;
    size_t n_kmers =
        strlen(sequence) - KMER_SIZE + 1; //todo :strlen can be pre-calculated

    //const Alphabet* alphabet = pore_model.pmalphabet;

    // Calculate summary statistics over the events and
    // the model implied by the read
    double event_level_sum = 0.0f; //do we need double?
    for (size_t i = 0; i < et.n; ++i) {
        event_level_sum += et.event[i].mean;
    }

    double kmer_level_sum = 0.0f;
    double kmer_level_sq_sum = 0.0f;
    for (size_t i = 0; i < n_kmers; ++i) {
        size_t kr = kmer_rank_function(&sequence[i], KMER_SIZE);
        double l = pore_model[kr].level_mean;
        kmer_level_sum += l;
        kmer_level_sq_sum += l * l;
    }

    double shift = event_level_sum / et.n - kmer_level_sum / n_kmers;

    // estimate scale
    double event_level_sq_sum = 0.0f;
    for (size_t i = 0; i < et.n; ++i) {
        event_level_sq_sum +=
            (et.event[i].mean - shift) * (et.event[i].mean - shift);
    }

    double scale = (event_level_sq_sum / et.n) / (kmer_level_sq_sum / n_kmers);

    //out.set4(shift, scale, 0.0, 1.0);
    out.shift = (float)shift;
    out.scale = (float)scale;

#ifdef DEBUG_PRINT_STATS
    fprintf(stderr, "event mean: %.2lf kmer mean: %.2lf shift: %.2lf\n",
            event_level_sum / et.n, kmer_level_sum / n_kmers, out.shift);
    fprintf(stderr, "event sq-mean: %.2lf kmer sq-mean: %.2lf scale: %.2lf\n",
            event_level_sq_sum / et.n, kmer_level_sq_sum / n_kmers, out.scale);
    //fprintf(stderr, "truth shift: %.2lf scale: %.2lf\n", pore_model.shift, pore_model.scale);
#endif
    return out;
}
AlignedPair* align(char* sequence,event_table events,model_t* models, scalings_t scaling){
// std::vector<AlignedPair> align(char* sequence,event_table events,model_t* models, scalings_t scaling){
    /*
     * This function should be tested now :-)
     */

    size_t strand_idx=0;
    size_t k = 6;

  // size_t n_events = events[strand_idx].n;
 size_t n_events = events.n;
  size_t n_kmers=strlen(sequence)-k+1;

  // backtrack markers
  const uint8_t FROM_D = 0;
  const uint8_t FROM_U = 1;
  const uint8_t FROM_L = 2;

  // qc
  double min_average_log_emission = -5.0;
  int max_gap_threshold = 50;

  // banding
  int bandwidth = 100;
  int half_bandwidth = bandwidth / 2;

  // transition penalties
  double events_per_kmer = (double)n_events / n_kmers;
  double p_stay = 1 - (1 / (events_per_kmer + 1));

  // setting a tiny skip penalty helps keep the true alignment within the adaptive band
  // this was empirically determined
  double epsilon = 1e-10;
  double lp_skip = log(epsilon);
  double lp_stay = log(p_stay);
  double lp_step = log(1.0 - exp(lp_skip) - exp(lp_stay));
  double lp_trim = log(0.01);

  // dp matrix
  size_t n_rows = n_events + 1;
  size_t n_cols = n_kmers + 1;
  size_t n_bands = n_rows + n_cols;

  // Initialize

  // Precompute k-mer ranks to avoid doing this in the inner loop
  size_t* kmer_ranks=(size_t*)malloc(sizeof(size_t)*n_kmers);

  for(size_t i = 0; i < n_kmers; ++i) {
      //>>>>>>>>> New replacement begin
      char* substring=&sequence[i];
      kmer_ranks[i] = kmer_rank_function(substring, k);
      //<<<<<<<<< New replacement over
  }

  float** bands=(float**)malloc(sizeof(float*) * n_bands);
  uint8_t** trace=(uint8_t**)malloc(sizeof(uint8_t*) * n_bands);

  for(size_t i=0;i<n_bands;i++){
    bands[i]=(float*)malloc(sizeof(float)*bandwidth);
    trace[i]=(uint8_t*)malloc(sizeof(uint8_t)*bandwidth);

    for(size_t j=0;j<bandwidth;j++){
      bands[i][j]=-INFINITY;
      trace[i][j]=0;
    }
  }


  // Keep track of the event/kmer index for the lower left corner of the band
  // these indices are updated at every iteration to perform the adaptive banding
  // Only the first two bands have their coordinates initialized, the rest are computed adaptively




  struct EventKmerPair
  {
      int event_idx;
      int kmer_idx;
  };
  //>>>>>>>>>>>>>>>>>New Replacement Begin
  struct EventKmerPair* band_lower_left=(struct EventKmerPair*) malloc(sizeof(struct EventKmerPair) * n_bands);
  //std::vector<EventKmerPair> band_lower_left(n_bands);
  //<<<<<<<<<<<<<<<<<New Replacement over


  // initialize range of first two bands
  band_lower_left[0].event_idx = half_bandwidth - 1;
  band_lower_left[0].kmer_idx = -1 - half_bandwidth;
  band_lower_left[1] = move_down(band_lower_left[0]);


  int start_cell_offset = band_kmer_to_offset(0, -1);
  assert(is_offset_valid(start_cell_offset));
  assert(band_event_to_offset(0, -1) == start_cell_offset);
  bands[0][start_cell_offset] = 0.0f;

  // band 1: first event is trimmed
  int first_trim_offset = band_event_to_offset(1, 0);
  assert(kmer_at_offset(1, first_trim_offset) == -1);
  assert(is_offset_valid(first_trim_offset));
  bands[1][first_trim_offset] = lp_trim;
  trace[1][first_trim_offset] = FROM_U;



      int fills = 0;
  #ifdef DEBUG_ADAPTIVE
      fprintf(stderr, "[trim] bi: %d o: %d e: %d k: %d s: %.2lf\n", 1, first_trim_offset, 0, -1, bands[1][first_trim_offset]);
  #endif

      // fill in remaining bands
    for(int band_idx = 2; band_idx < n_bands; ++band_idx) {
        // Determine placement of this band according to Suzuki's adaptive algorithm
        // When both ll and ur are out-of-band (ob) we alternate movements
        // otherwise we decide based on scores
        float ll = bands[band_idx - 1][0];
        float ur = bands[band_idx - 1][bandwidth - 1];
        bool ll_ob = ll == -INFINITY;
        bool ur_ob = ur == -INFINITY;

        bool right = false;
        if(ll_ob && ur_ob) {
            right = band_idx % 2 == 1;
        } else {
            right = ll < ur; // Suzuki's rule
        }

        if(right) {
            band_lower_left[band_idx] = move_right(band_lower_left[band_idx - 1]);
        } else {
            band_lower_left[band_idx] = move_down(band_lower_left[band_idx - 1]);
        }
        // If the trim state is within the band, fill it in here
        int trim_offset = band_kmer_to_offset(band_idx, -1);
        if(is_offset_valid(trim_offset)) {
            int event_idx = event_at_offset(band_idx, trim_offset);
            if(event_idx >= 0 && event_idx < n_events) {
                bands[band_idx][trim_offset] = lp_trim * (event_idx + 1);
                trace[band_idx][trim_offset] = FROM_U;
            } else {
                bands[band_idx][trim_offset] = -INFINITY;
            }
        }

        // Get the offsets for the first and last event and kmer
        // We restrict the inner loop to only these values
        int kmer_min_offset = band_kmer_to_offset(band_idx, 0);
        int kmer_max_offset = band_kmer_to_offset(band_idx, n_kmers);
        int event_min_offset = band_event_to_offset(band_idx, n_events - 1);
        int event_max_offset = band_event_to_offset(band_idx, -1);

        int min_offset = std::max(kmer_min_offset, event_min_offset);
        min_offset = std::max(min_offset, 0);

        int max_offset = std::min(kmer_max_offset, event_max_offset);
        max_offset = std::min(max_offset, bandwidth);

        for(int offset = min_offset; offset < max_offset; ++offset) {
            int event_idx = event_at_offset(band_idx, offset);
            int kmer_idx = kmer_at_offset(band_idx, offset);

            size_t kmer_rank = kmer_ranks[kmer_idx];

            int offset_up   = band_event_to_offset(band_idx - 1, event_idx - 1);
            int offset_left = band_kmer_to_offset(band_idx - 1, kmer_idx - 1);
            int offset_diag = band_kmer_to_offset(band_idx - 2, kmer_idx - 1);

#ifdef DEBUG_ADAPTIVE
            // verify loop conditions
            assert(kmer_idx >= 0 && kmer_idx < n_kmers);
            assert(event_idx >= 0 && event_idx < n_events);
            assert(offset_diag == band_event_to_offset(band_idx - 2, event_idx - 1));
            assert(offset_up - offset_left == 1);
            assert(offset >= 0 && offset < bandwidth);
#endif

            float up   = is_offset_valid(offset_up)   ? bands[band_idx - 1][offset_up]   : -INFINITY;
            float left = is_offset_valid(offset_left) ? bands[band_idx - 1][offset_left] : -INFINITY;
            float diag = is_offset_valid(offset_diag) ? bands[band_idx - 2][offset_diag] : -INFINITY;

            float lp_emission = log_probability_match_r9(scaling, models, events, event_idx, kmer_rank,strand_idx);


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
            fprintf(stderr, "[adafill] offset-up: %d offset-diag: %d offset-left: %d\n", offset_up, offset_diag, offset_left);
            fprintf(stderr, "[adafill] up: %.2lf diag: %.2lf left: %.2lf\n", up, diag, left);
            fprintf(stderr, "[adafill] bi: %d o: %d e: %d k: %d s: %.2lf f: %d emit: %.2lf\n", band_idx, offset, event_idx, kmer_idx, max_score, from, lp_emission);
#endif
            bands[band_idx][offset] = max_score;
            trace[band_idx][offset] = from;
            fills += 1;
        }
    }

    //
    // Backtrack to compute alignment
    //
    double sum_emission = 0;
    double n_aligned_events = 0;

    //>>>>>>>>>>>>>> New replacement begin
    std::vector<AlignedPair> out;
    int capacity = 1000000;
    AlignedPair* out_2=(AlignedPair*)malloc(sizeof(AlignedPair)*capacity);
    int outIndex=0;
    //<<<<<<<<<<<<<<<<New Replacement over

    float max_score = -INFINITY;
    int curr_event_idx = 0;
    int curr_kmer_idx = n_kmers -1;

    // Find best score between an event and the last k-mer. after trimming the remaining evnets
    for(int event_idx = 0; event_idx < n_events; ++event_idx) {
        int band_idx = event_kmer_to_band(event_idx, curr_kmer_idx);

        //>>>>>>>New  replacement begin
        /*assert(band_idx < bands.size());*/
        assert(band_idx < n_bands);
        //<<<<<<<<New Replacement over
        int offset = band_event_to_offset(band_idx, event_idx);
        if(is_offset_valid(offset)) {
            float s = bands[band_idx][offset] + (n_events - event_idx) * lp_trim;
            if(s > max_score) {
                max_score = s;
                curr_event_idx = event_idx;
            }
        }
    }

#ifdef DEBUG_ADAPTIVE
    fprintf(stderr, "[adaback] ei: %d ki: %d s: %.2f\n", curr_event_idx, curr_kmer_idx, max_score);
#endif

    int curr_gap = 0;
    int max_gap = 0;
    while(curr_kmer_idx >= 0 && curr_event_idx >= 0) {

        // emit alignment
        //>>>>>>>New Repalcement begin
        out_2[outIndex].ref_pos = curr_kmer_idx;
        out_2[outIndex].read_pos = curr_event_idx;
        outIndex++;
        out.push_back({curr_kmer_idx, curr_event_idx});
        //<<<<<<<<<New Replacement over

#ifdef DEBUG_ADAPTIVE
        fprintf(stderr, "[adaback] ei: %d ki: %d\n", curr_event_idx, curr_kmer_idx);
#endif
        // qc stats
        //>>>>>>>>>>>>>>New Replacement begin
        char* substring=&sequence[curr_kmer_idx];
        size_t kmer_rank = kmer_rank_function(sequence, k);
        //<<<<<<<<<<<<<New Replacement over
        sum_emission += log_probability_match_r9(scaling, models, events, kmer_rank, curr_event_idx,0);
        n_aligned_events += 1;

        int band_idx = event_kmer_to_band(curr_event_idx, curr_kmer_idx);
        int offset = band_event_to_offset(band_idx, curr_event_idx);
        assert(band_kmer_to_offset(band_idx, curr_kmer_idx) == offset);

        uint8_t from = trace[band_idx][offset];
        if(from == FROM_D) {
            curr_kmer_idx -= 1;
            curr_event_idx -= 1;
            curr_gap = 0;
        } else if(from == FROM_U) {
            curr_event_idx -= 1;
            curr_gap = 0;
        } else {
            curr_kmer_idx -= 1;
            curr_gap += 1;
            max_gap = std::max(curr_gap, max_gap);
        }
    }

    //>>>>>>>>New replacement begin
    std::reverse(out.begin(), out.end());
    int c;
    int end = outIndex-1;
    for (c = 0; c < outIndex/2; c++) {
      int ref_pos_temp  = out_2[c].ref_pos;
      int read_pos_temp = out[c].read_pos;
      out_2[c].ref_pos  = out_2[end].ref_pos;
      out_2[c].read_pos = out[end].read_pos;
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
    bool spanned = out_2[0].ref_pos == 0 && out_2[outIndex-1].ref_pos == n_kmers - 1;
    // bool spanned = out.front().ref_pos == 0 && out.back().ref_pos == n_kmers - 1;
    //<<<<<<<<<<<<<New replacement over
    bool failed = false;
    if(avg_log_emission < min_average_log_emission || !spanned || max_gap > max_gap_threshold) {
        failed = true;
        //>>>>>>>>>>>>>New replacement begin
        //outIndex=0;
        out.clear();
        free(out_2);
        //<<<<<<<<<<<<<New replacement over
    }


    free(kmer_ranks);
    for(size_t i=0;i<n_bands;i++){
      free(bands[i]);
      free(trace[i]);
    }
    free(bands);
    free(trace);
    free(band_lower_left);
    //fprintf(stderr, "ada\t%s\t%s\t%.2lf\t%zu\t%.2lf\t%d\t%d\t%d\n", read.read_name.substr(0, 6).c_str(), failed ? "FAILED" : "OK", events_per_kmer, sequence.size(), avg_log_emission, curr_event_idx, max_gap, fills);
    //outSize=outIndex;
    //if(outIndex>500000)fprintf(stderr, "Max outSize %d\n", outIndex);
    return out_2;







}


// float log_probability_match_r9(scalings_t scaling,
//                                       model_t* models,
//                                       event_table events,
//                                       int event_idx,
//                                       uint32_t kmer_rank,
//                                       uint8_t strand)
// {
//     // event level mean, scaled with the drift value
//     strand=0;

//     //float level = read.get_drift_scaled_level(event_idx, strand);

//     float time=events[event_idx].start- events[0].start;
//     float unscaledLevel=events[event_idx].mean;
//     float scaledLevel=unscaledLevel - time* scaling->shift;

//     //GaussianParameters gp = read.get_scaled_gaussian_from_pore_model_state(pore_model, strand, kmer_rank);
//     float gp_mean = scaling.scale * models[kmer_rank].level_mean + scaling.shift;
//     float gp_stdv = model[kmer_rank].level_stdv * scaling.var;
//     float gp_log_stdv = model[kmer_rank].level_log_stdv + scaling.log_var;


//     float lp = log_normal_pdf(scaledLevel, gp_mean,gp_stdv,gp_log_stdv);
//     return lp;
// }
