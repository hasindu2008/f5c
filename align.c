#include "f5c.h"
#include <assert.h>
//#define DEBUG_PRINT_STATS 1

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
        size_t kr = get_kmer_rank(&sequence[i], KMER_SIZE);
        double l = pore_model[kr].level_mean;
        //fprintf(stderr,"Kmer : %c%c%c%c%c%c, kmer_rank : %d , kmer_mean : %f \n",sequence[i],sequence[i+1],sequence[i+2],sequence[i+3],sequence[i+4],sequence[i+5],kr,l);
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

AlignedPair* align(char* sequence, event_table events, model_t* models,
                   scalings_t scaling) {
    return (NULL);
}

event_alignment_t* postalign(char* sequence, AlignedPair* event_alignment,
                             int32_t n_events) {
    if (n_events > 0) {
        /* transform alignment into the base-to-event map*/

        int32_t n_kmers = strlen(sequence) - KMER_SIZE +
                          1; //todo :strlen can be pre-calculated

        // create base-to-event map
        index_pair_t* base_to_event_map =
            (index_pair_t*)(malloc(sizeof(index_pair_t) * n_kmers));
        MALLOC_CHK(base_to_event_map);

        //initialisesing (todo : check if really required)
        int32_t i = 0;
        for (i = 0; i < n_kmers; i++) {
            base_to_event_map[i].start = -1;
            base_to_event_map[i].stop = -1;
        }

        int32_t max_event = 0;
        int32_t min_event = INT32_MAX;

        int32_t prev_event_idx = -1;

        for (i = 0; i < n_events; ++i) {
            int32_t k_idx = event_alignment[i].ref_pos;
            int32_t event_idx = event_alignment[i].read_pos;
            index_pair_t elem = base_to_event_map[k_idx];
            if (event_idx != prev_event_idx) {
                if (elem.start == -1) {
                    elem.start = event_idx;
                }
                elem.stop = event_idx;
            }

            max_event = max_event > event_idx ? max_event : event_idx;
            min_event = min_event < event_idx ? min_event : event_idx;
            prev_event_idx = event_idx;
        }

        //events_per_base[strand_idx] = (double)(max_event - min_event) / n_kmers;

        /*prepare data structures for the final calibration*/

        event_alignment_t* alignment =
            (event_alignment_t*)malloc(sizeof(event_alignment_t) * n_kmers);
        MALLOC_CHK(alignment);
        int32_t alignment_index = 0;

        int32_t prev_kmer_rank = -1;

        int32_t ki;
        for (ki = 0; ki < n_kmers; ++ki) {
            index_pair_t event_range_for_kmer = base_to_event_map[ki];

            // skip kmers without events
            if (event_range_for_kmer.start == -1)
                continue;

            // skip k-mers that cannot be shifted to a valid position
            // if(ki + shift_offset < 0 || ki + shift_offset >= n_kmers) {
            //     continue;
            // }

            for (int32_t event_idx = event_range_for_kmer.start;
                 event_idx <= event_range_for_kmer.stop; event_idx++) {
                assert(event_idx < n_events);

                // since we use the 1D read seqence here we never have to reverse complement

                int32_t kmer_rank = get_kmer_rank(&sequence[ki], KMER_SIZE);

                event_alignment_t ea;
                // ref data
                //ea.ref_name = "read";
                ea.read_idx = -1; // not needed
                kmer_cpy(ea.ref_kmer, &sequence[ki], KMER_SIZE);
                ea.ref_position = ki;
                //ea.strand_idx = strand_idx;
                ea.event_idx = event_idx;
                ea.rc = false;
                kmer_cpy(ea.model_kmer, &sequence[ki], KMER_SIZE);
                ea.hmm_state = prev_kmer_rank != kmer_rank ? 'M' : 'E';
                if (alignment_index > n_events) {
                    ERROR("We have run out of space in event_alignment_t* "
                          "alignment. Assumption fialed. Current size %d",
                          n_events);
                    exit(EXIT_FAILURE);
                }
                alignment[alignment_index] = ea;
                alignment_index++;
                prev_kmer_rank = kmer_rank;
            }
        }

        free(base_to_event_map);
        return alignment;
    }

    else {
        return NULL;
    }
}

// recalculate shift, scale, drift, scale_sd from an alignment and the read
// returns true if the recalibration was performed
// in either case, sets residual to the L1 norm of the residual
bool recalibrate_model(model_t* pore_model, event_table et,
                       scalings_t* scallings,
                       const event_alignment_t* alignment_output,
                       int32_t num_alignments, bool scale_var) {
    //std::vector<double> raw_events, times, level_means, level_stdvs;

    //count the number of M states. Then can do on the fly without mallocs.

    //std::cout << "Previous pore model parameters: " << sr.pore_model[strand_idx].shift << ", "
    //                                                << sr.pore_model[strand_idx].scale << ", "
    //                                                << sr.pore_model[strand_idx].drift << ", "
    //                                                << sr.pore_model[strand_idx].var   << std::endl;

    // extract necessary vectors from the read and the pore model; note do not want scaled values
    int32_t num_M_state = 0;
    for (int32_t ei = 0; ei < num_alignments; ++ei) {
        event_alignment_t ea = alignment_output[ei];
        if (ea.hmm_state == 'M') {
            num_M_state++;
            //
            //fprintf(stdout, "recalibrate ei: %zu level: %.2lf kmer: %s model: %.2lf\n",
            //        ei, sr.get_uncorrected_level(ea.event_idx, strand_idx), model_kmer.c_str(),
            //        sr.pore_model[strand_idx].states[rank].level_mean);
            //
        }
    }

    const int32_t minNumEventsToRescale = 200;
    bool recalibrated = false;
    if (num_M_state >= minNumEventsToRescale) {
        // Assemble linear system corresponding to weighted least squares problem
        // Can just directly call a weighted least squares solver, but there's enough
        // structure in our problem it's a little faster just to build the normal eqn
        // matrices ourselves

        double A00 = 0, A01 = 0, A10 = 0, A11 = 0;
        double b0 = 0, b1 = 0;
        double x0 = 0, x1 = 0;

        for (int32_t ei = 0; ei < num_alignments; ++ei) {
            event_alignment_t ea = alignment_output[ei];
            if (ea.hmm_state == 'M') {
                //std::string model_kmer = ea.rc ? pore_model.pmalphabet->reverse_complement(ea.ref_kmer) : ea.ref_kmer;
                uint32_t rank = get_kmer_rank(ea.ref_kmer, KMER_SIZE);

                double raw_event = et.event[ea.event_idx].mean;
                double level_mean = pore_model[rank].level_mean;
                double level_stdv = pore_model[rank].level_stdv;

                double inv_var = 1. / (level_stdv * level_stdv);
                double mu = level_mean;
                double e = raw_event;

                A00 += inv_var;
                A01 += mu * inv_var;
                A11 += mu * mu * inv_var;

                b0 += e * inv_var;
                b1 += mu * e * inv_var;
            }
        }

        A10 = A01;

        // perform the linear solve
        //Eigen::VectorXd x = A.fullPivLu().solve(b);
        double div = A00 * A11 - A01 * A10;
        x0 = -(A01 * b1 - A11 * b0) / div;
        x1 = (A00 * b1 - A10 * b0) / div;

        double shift = x0;
        double scale = x1;
        //double drift = 0.;
        double var = 1.0;

        if (scale_var) {
            var = 0.;
            for (int32_t ei = 0; ei < num_alignments; ++ei) {
                event_alignment_t ea = alignment_output[ei];
                if (ea.hmm_state == 'M') {
                    uint32_t rank = get_kmer_rank(ea.ref_kmer, KMER_SIZE);
                    double raw_event = et.event[ea.event_idx].mean;
                    double level_mean = pore_model[rank].level_mean;
                    double level_stdv = pore_model[rank].level_stdv;
                    double yi = (raw_event - shift - scale * level_mean);
                    var += yi * yi / (level_stdv * level_stdv);
                }
            }
            var /= num_M_state;
            var = sqrt(var);
        }

        scallings->shift = shift;
        scallings->scale = scale;
        //scallings->drift=drift;
        scallings->var = var;

        recalibrated = true;
    }

    return recalibrated;
}
