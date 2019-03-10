#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>

#include "f5c.h"
#include "f5cmisc.h"

struct EventAlignmentParameters
{
    EventAlignmentParameters()
    {
        //sr = NULL;
        fai = NULL;
        hdr = NULL;
        record = NULL;
        //strand_idx = NUM_STRANDS;
        
        alphabet = "";
        read_idx = -1;
        region_start = -1;
        region_end = -1;
    }


    // Mandatory
    //SquiggleRead* sr;
    const faidx_t* fai;
    const bam_hdr_t* hdr;
    const bam1_t* record;
    //size_t strand_idx;
    
    // optional
    std::string alphabet;
    int read_idx;
    int region_start;
    int region_end;
};



std::vector<EventAlignment> align_read_to_ref(const EventAlignmentParameters& params)
{
    // Sanity check input parameters
    //assert(params.sr != NULL);
    assert(params.fai != NULL);
    assert(params.hdr != NULL);
    assert(params.record != NULL);
    //assert(params.strand_idx < NUM_STRANDS);
    assert( (params.region_start == -1 && params.region_end == -1) || (params.region_start <= params.region_end));
    //const PoreModel* pore_model = params.get_model();

    std::vector<EventAlignment> alignment_output;

    // Extract the reference subsequence for the entire alignment
    // int fetched_len = 0;
    // int ref_offset = params.record->core.pos;
    // std::string ref_name(params.hdr->target_name[params.record->core.tid]);
    // std::string ref_seq = get_reference_region_ts(params.fai, ref_name.c_str(), ref_offset, 
    //                                               bam_endpos(params.record), &fetched_len);

    // convert to upper case
    std::transform(ref_seq.begin(), ref_seq.end(), ref_seq.begin(), ::toupper);
    
    // k from read pore model
    const uint32_t k = KMER_SIZE;

    // If the reference sequence contains ambiguity codes
    // switch them to the lexicographically lowest base
    ref_seq =disambiguate(ref_seq);
    std::string rc_ref_seq = reverse_complement(ref_seq);

    // Skip unmapped
    // if((params.record->core.flag & BAM_FUNMAP) != 0) {
    //     return alignment_output;
    // }

    // Get the read-to-reference aligned segments
    std::vector<AlignedSegment> aligned_segments = get_aligned_segments(params.record);
    for(size_t segment_idx = 0; segment_idx < aligned_segments.size(); ++segment_idx) {

        AlignedSegment& aligned_pairs = aligned_segments[segment_idx];

        if(params.region_start != -1 && params.region_end != -1) {
            trim_aligned_pairs_to_ref_region(aligned_pairs, params.region_start, params.region_end);
        }

        // Trim the aligned pairs to be within the range of the maximum kmer index
        int max_kmer_idx = params.sr->read_sequence.size() - k;
        trim_aligned_pairs_to_kmer(aligned_pairs, max_kmer_idx);

        if(aligned_pairs.empty()) {
            return alignment_output;
        }

        bool do_base_rc = bam_is_rev(params.record);
        bool rc_flags[2] = { do_base_rc, !do_base_rc }; // indexed by strand
        const int align_stride = 100; // approximately how many reference bases to align to at once
        const int output_stride = 50; // approximately how many event alignments to output at once

        // get the event range of the read to re-align
        int read_kidx_start = aligned_pairs.front().read_pos;
        int read_kidx_end = aligned_pairs.back().read_pos;
        
        if(do_base_rc) {
            read_kidx_start = params.sr->flip_k_strand(read_kidx_start, k);
            read_kidx_end = params.sr->flip_k_strand(read_kidx_end, k);
        }
        
        assert(read_kidx_start >= 0);
        assert(read_kidx_end >= 0);

        int first_event = params.sr->get_closest_event_to(read_kidx_start, params.strand_idx);
        int last_event = params.sr->get_closest_event_to(read_kidx_end, params.strand_idx);
        bool forward = first_event < last_event;

        int curr_start_event = first_event;
        int curr_start_ref = aligned_pairs.front().ref_pos;
        int curr_pair_idx = 0;

        while( (forward && curr_start_event < last_event) ||
               (!forward && curr_start_event > last_event)) {

            // Get the index of the aligned pair approximately align_stride away
            int end_pair_idx = get_end_pair(aligned_pairs, curr_start_ref + align_stride, curr_pair_idx);
        
            int curr_end_ref = aligned_pairs[end_pair_idx].ref_pos;
            int curr_end_read = aligned_pairs[end_pair_idx].read_pos;

            if(do_base_rc) {
                curr_end_read = params.sr->flip_k_strand(curr_end_read, k);
            }
            assert(curr_end_read >= 0);

            int s = curr_start_ref - ref_offset;
            int l = curr_end_ref - curr_start_ref + 1;

            std::string fwd_subseq = ref_seq.substr(s, l);
            std::string rc_subseq = rc_ref_seq.substr(ref_seq.length() - s - l, l);
            assert(fwd_subseq.length() == rc_subseq.length());

            HMMInputSequence hmm_sequence(fwd_subseq, rc_subseq, pore_model->pmalphabet);
            
            // Require a minimum amount of sequence to align to
            if(hmm_sequence.length() < 2 * k) {
                break;
            }

            // Set up HMM input
            HMMInputData input;
            input.read = params.sr;
            input.pore_model = pore_model;
            assert(input.pore_model != NULL);

            input.event_start_idx = curr_start_event;
            input.event_stop_idx = params.sr->get_closest_event_to(curr_end_read, params.strand_idx);
            //printf("[SEGMENT_START] read: %s event start: %zu event end: %zu\n", params.sr->read_name.c_str(), input.event_start_idx, input.event_stop_idx);

            // A limitation of the segment-by-segment alignment is that we can't jump
            // over very large deletions wrt to the reference. The effect of this
            // is that we can get segments that have very few alignable events. We
            // just stop processing them for now
            if(abs((int)input.event_start_idx - (int)input.event_stop_idx) < 2)
                break;

            input.strand = params.strand_idx;
            input.event_stride = input.event_start_idx < input.event_stop_idx ? 1 : -1;
            input.rc = rc_flags[params.strand_idx];

            std::vector<HMMAlignmentState> event_alignment = profile_hmm_align(hmm_sequence, input);
            
            // Output alignment
            size_t num_output = 0;
            size_t event_align_idx = 0;

            // If we aligned to the last event, output everything and stop
            bool last_section = end_pair_idx == (int)aligned_pairs.size() - 1;

            /*
            // Don't allow the segment to end on an E state or else we get alignment
            // artifacts at the segment boundary
            if(!last_section) {
                size_t last_match_index = event_alignment.size() - 1;
                while(event_alignment[last_match_index].state != 'M') {
                    last_match_index -= 1;
                }

                event_alignment.resize(last_match_index + 1);
                if(event_alignment.empty()) {
                    break;
                }
                assert(event_alignment.back().state == 'M');
            }
            */

            int last_event_output = 0;
            int last_ref_kmer_output = 0;

            for(; event_align_idx < event_alignment.size() && 
                  (num_output < output_stride || last_section); event_align_idx++) {

                HMMAlignmentState& as = event_alignment[event_align_idx];
                if(as.state != 'K' && (int)as.event_idx != curr_start_event) {

                    EventAlignment ea;
                    
                    // ref
                    ea.ref_name = ref_name;
                    ea.ref_position = curr_start_ref + as.kmer_idx;
                    ea.ref_kmer = ref_seq.substr(ea.ref_position - ref_offset, k);

                    // event
                    ea.read_idx = params.read_idx;
                    ea.strand_idx = params.strand_idx;
                    ea.event_idx = as.event_idx;
                    ea.rc = input.rc;

                    // hmm
                    ea.hmm_state = as.state;

                    if(ea.hmm_state != 'B') {
                        ea.model_kmer = hmm_sequence.get_kmer(as.kmer_idx, k, input.rc);
                    } else {
                        ea.model_kmer = std::string(k, 'N');
                    }

                    // store
                    alignment_output.push_back(ea);

                    // update
                    last_event_output = as.event_idx;
                    last_ref_kmer_output = curr_start_ref + as.kmer_idx;

                    num_output += 1;
                }
            }

            // Advance the pair iterator to the ref base
            curr_start_event = last_event_output;
            curr_start_ref = last_ref_kmer_output;
            //printf("[SEGMENT_END] read: %s last event output: %zu ref pos: %zu (%s)\n", params.sr->read_name.c_str(), last_event_output, last_ref_kmer_output, ref_seq.substr(last_ref_kmer_output - ref_offset, k).c_str());
            curr_pair_idx = get_end_pair(aligned_pairs, curr_start_ref, curr_pair_idx);

#if EVENTALIGN_TRAIN
            // update training data for read
            params.sr->parameters[params.strand_idx].add_training_from_alignment(hmm_sequence, input, event_alignment);
            global_training[params.strand_idx].add_training_from_alignment(hmm_sequence, input, event_alignment);
#endif

            if(num_output == 0) {
                break;
            }
        } // for realignmentsegment
    } // for bam aligned segment

    return alignment_output;
}


// Realign the read in event space
void realign_read(char* ref,
                  const faidx_t* fai,
                  const bam_hdr_t* hdr,
                  const bam1_t* record,
                  size_t read_idx,
                  int region_start,
                  int region_end, event_table* events)
{
    // Load a squiggle read for the mapped read
    std::string read_name = bam_get_qname(record);

    // load read
    //SquiggleRead sr(read_name, read_db, opt::write_samples ? SRF_LOAD_RAW_SAMPLES : 0);

    if(opt::verbose > 1) {
        fprintf(stderr, "Realigning %s [%zu]\n",
                read_name.c_str(), sr.events.n);
    }

    //for(int strand_idx = 0; strand_idx < 2; ++strand_idx) {

        // Do not align this strand if it was not sequenced
        // if(!sr.has_events_for_strand(strand_idx)) {
        //     continue;
        // }

        EventAlignmentParameters params;
        //params.sr = &sr;
        params.fai = fai;
        params.hdr = hdr;
        params.record = record;
       // params.strand_idx = strand_idx;
        
        params.read_idx = read_idx;
        params.region_start = region_start;
        params.region_end = region_end;

        std::vector<EventAlignment> alignment = align_read_to_ref(params);

        // EventalignSummary summary;
        // //if(writer.summary_fp != NULL) {
        //     summary = summarize_alignment(sr, strand_idx, params, alignment);
        // //}

        // if(opt::output_sam) {
        //     emit_event_alignment_sam(writer.sam_fp, sr, hdr, record, alignment);
        // } else {
        //     emit_event_alignment_tsv(writer.tsv_fp, sr, strand_idx, params, alignment);
        // }

        // if(writer.summary_fp != NULL && summary.num_events > 0) {
        //     assert(params.alphabet == "");
        //     const PoreModel* pore_model = params.get_model();
        //     SquiggleScalings& scalings = sr.scalings[strand_idx];
        //     fprintf(writer.summary_fp, "%zu\t%s\t%s\t", read_idx, read_name.c_str(), sr.fast5_path.c_str());
        //     fprintf(writer.summary_fp, "%s\t%s\t", pore_model->name.c_str(), strand_idx == 0 ? "template" : "complement");
        //     fprintf(writer.summary_fp, "%d\t%d\t%d\t%d\t", summary.num_events, summary.num_steps, summary.num_skips, summary.num_stays);
        //     fprintf(writer.summary_fp, "%.2lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\n", summary.sum_duration, scalings.shift, scalings.scale, scalings.drift, scalings.var);
        // }
        
    //}
}