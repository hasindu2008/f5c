
Conversation opened. 42 messages. All messages read.

Skip to content
Using Gmail with screen readers
2 of 6,197
Porting adaptive_banded_simple_event_align [was: Re: Biocomputing/algorithm/coding work]
Inbox
	x
Roshan G. Ragel
	
	Jul 20, 2018, 8:01 AM
Hi Gihan and Hiruna, The tool we are talking about is one called nanopolish. IIIR, it is used in a process called methylation in the genomic pipeline. If you ar
39
Hasindu Gamaarachchi
	
	Mar 25, 2019, 11:31 AM
Thanks, I will check as soon as I get some free time.
Hasindu Gamaarachchi
	
Attachments8:48 AM (1 hour ago)
	
to me
Hi Hiruna,

I did a few changes so that your eventalign.c is compatible with the current repository. Could you send a pull request to the https://github.com/hasindu2008/f5c/tree/eventalign
so that you get the necessary credit for the work you did.
fork a f5c
checkout to eventalign branch
put this attached eventalign.c  
commit
and send a pull request to the eventalign branch (not the master)
Attachments area
	
Sure, I will do that.
Done.
No, I can't.
	
	

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include "matrix.h"
#include <algorithm>

#include "f5c.h"
#include "f5cmisc.h"


#define METHYLATED_SYMBOL 'M'


typedef std::vector<AlignedPair> AlignedSegment;



enum ProfileStateR9
{
    PSR9_KMER_SKIP = 0,
    PSR9_BAD_EVENT,
    PSR9_MATCH,
    PSR9_NUM_STATES = 3,
    PSR9_PRE_SOFT // intentionally after PS_NUM_STATES
};

enum HMMMovementType
{
    HMT_FROM_SAME_M = 0,
    HMT_FROM_PREV_M,
    HMT_FROM_SAME_B,
    HMT_FROM_PREV_B,
    HMT_FROM_PREV_K,
    HMT_FROM_SOFT,
    HMT_NUM_MOVEMENT_TYPES
};
typedef struct { float x[HMT_NUM_MOVEMENT_TYPES]; } HMMUpdateScores;

// from https://github.com/jts/nanopolish/blob/3180474dc1a79c076a21e70bb4d19d2705b39b56/src/common/nanopolish_common.h
// A representation of an event->kmer alignment
struct HMMAlignmentState
{
    uint32_t event_idx;
    uint32_t kmer_idx;
    double l_posterior;
    double l_fm;
    double log_transition_probability;
    char state;
};

static inline void profile_hmm_forward_initialize_r9(FloatMatrix& fm)
{
    // initialize forward calculation
    for(uint32_t si = 0; si < fm.n_cols; si++) {
        set(fm, 0, si, -INFINITY);
    }

    for(uint32_t ri = 0; ri < fm.n_rows; ri++) {
        set(fm, ri, PSR9_KMER_SKIP, -INFINITY);
        set(fm, ri, PSR9_BAD_EVENT, -INFINITY);
        set(fm, ri, PSR9_MATCH, -INFINITY);
    }
}

// Convert an enumerated state into a symbol
static inline char ps2char(ProfileStateR9 ps) { return "KBMNS"[ps]; }


static inline std::vector<HMMAlignmentState> profile_hmm_align(std::string fwd_subseq,std::string rc_subseq,
    uint32_t k, //KMER_SIZE
    uint32_t e_start, //curr_start_event;
    uint32_t e_end, //input_event_stop_idx
    int8_t event_stride) // event_stride )
{
    std::vector<HMMAlignmentState> alignment;
    uint32_t n_kmers = fwd_subseq.length() - k + 1;
    uint32_t n_states = PSR9_NUM_STATES * (n_kmers + 2); // + 2 for explicit terminal states
    uint32_t n_events = 0;
    if(e_end > e_start)
        n_events = e_end - e_start + 1;
    else
        n_events = e_start - e_end + 1;
    assert(n_events >= 2);
    uint32_t n_rows = n_events + 1;
    // Allocate a matrix to hold the HMM result
    FloatMatrix vm;
    allocate_matrix(vm, n_rows, n_states);
    UInt8Matrix bm;
    allocate_matrix(bm, n_rows, n_states);

    profile_hmm_forward_initialize_r9(vm);
   //? profile_hmm_fill_generic_r9(sequence, data, e_start, flags, output);
     // Traverse the backtrack matrix to compute the results
    int traversal_stride = event_stride;

#if HMM_REVERSE_FIX
    // Hack to support the fixed HMM
    // TODO: clean up
    traversal_stride = 1;
    if(event_stride == -1) {
        e_start = event_stop_idx;
    }
#endif

    // start from the last event matched to the last kmer
    uint32_t row = n_rows - 1;
    uint32_t col = PSR9_NUM_STATES * n_kmers + PSR9_MATCH;

    while(row > 0) {
        
        uint32_t event_idx = e_start + (row - 1) * traversal_stride;
        uint32_t block = col / PSR9_NUM_STATES;
        uint32_t kmer_idx = block - 1;
        ProfileStateR9 curr_ps = (ProfileStateR9) (col % PSR9_NUM_STATES);

#if DEBUG_BACKTRACK
        printf("backtrace %zu %zu coord: (%zu, %zu, %zu) state: %d\n", event_idx, kmer_idx, row, col, block, curr_ps);
#endif
        assert(block > 0);
        assert(get(vm, row, col) != -INFINITY);

        HMMAlignmentState as;
        as.event_idx = event_idx;
        as.kmer_idx = kmer_idx;
        as.l_posterior = -INFINITY; // not computed
        as.l_fm = get(vm, row, col);
        as.log_transition_probability = -INFINITY; // not computed
        as.state = ps2char(curr_ps);
        alignment.push_back(as);

        // Update the event (row) and k-mer using the backtrack matrix
        HMMMovementType movement = (HMMMovementType)get(bm, row, col);
        if(movement == HMT_FROM_SOFT) {
            break;
        }

         // update kmer_idx and state
        ProfileStateR9 next_ps;
        switch(movement) {
            case HMT_FROM_SAME_M: 
                next_ps = PSR9_MATCH;
                break;
            case HMT_FROM_PREV_M: 
                kmer_idx -= 1;
                next_ps = PSR9_MATCH;
                break;
            case HMT_FROM_SAME_B:
                next_ps = PSR9_BAD_EVENT;
                break;
            case HMT_FROM_PREV_B:
                kmer_idx -= 1;
                next_ps = PSR9_BAD_EVENT;
                break;
            case HMT_FROM_PREV_K:
                kmer_idx -= 1;
                next_ps = PSR9_KMER_SKIP;
                break;
            case HMT_FROM_SOFT:
                assert(false);
                break;
        }

         // update row (event) idx only if this isn't a kmer skip, which is silent
        if(curr_ps != PSR9_KMER_SKIP) {
            row -= 1;
        }

        col = PSR9_NUM_STATES * (kmer_idx + 1) + next_ps;
    }

#if HMM_REVERSE_FIX
// change the strand of the kmer indices if we aligned to the reverse strand
if(event_stride == -1) {
    for(size_t ai = 0; ai < alignment.size(); ++ai) {
        size_t k_idx = alignment[ai].kmer_idx;
        alignment[ai].kmer_idx = fwd_subseq.length() - k_idx - k;
    }
} else {
    std::reverse(alignment.begin(), alignment.end());
}
#else
    //
    std::reverse(alignment.begin(), alignment.end());
#endif

 //
    free_matrix(vm);
    free_matrix(bm);

    return alignment;

}




//from nanopolish eventalign.c
// Returns the index into the aligned_pairs vector that has the highest ref_pos
// that is not greater than ref_pos_max. It starts the search at pair_idx
static int get_end_pair(const std::vector<AlignedPair>& aligned_pairs, int ref_pos_max, int pair_idx)
{
    while(pair_idx < (int)aligned_pairs.size()) {
        if(aligned_pairs[pair_idx].ref_pos > ref_pos_max)
            return pair_idx - 1;
        pair_idx += 1;
    }

    return aligned_pairs.size() - 1;
}
//from nanopolish eventalign.c
// Modify the aligned_pairs vector to ensure there are no alignments
// outside of the given reference coordinates
static void trim_aligned_pairs_to_ref_region(std::vector<AlignedPair>& aligned_pairs, int ref_start, int ref_end)
{
    std::vector<AlignedPair> trimmed;
    for(size_t i = 0; i < aligned_pairs.size(); ++i) {
        if(aligned_pairs[i].ref_pos >= ref_start && 
           aligned_pairs[i].ref_pos <= ref_end) {
            trimmed.push_back(aligned_pairs[i]);
        }
    }
    
    aligned_pairs.swap(trimmed);
}
//from nanopolish eventalign.c
// Modify the aligned_pairs vector to ensure the highest read position
// does not exceed max_kmer
static void trim_aligned_pairs_to_kmer(std::vector<AlignedPair>& aligned_pairs, int max_kmer_idx)
{
    int idx = aligned_pairs.size() - 1;
    while(idx >= 0 && aligned_pairs[idx].read_pos > max_kmer_idx)
        idx -= 1;

    if(idx < 0)
        aligned_pairs.clear(); // no valid data
    else
        aligned_pairs.resize(idx + 1);
}


//from https://github.com/hasindu2008/f5c/blob/64d9d51cd773e46e1a72dbf2e2bf333be24df7c4/meth.c
//helper for get_closest_event_to
static inline int get_next_event(int start, int stop, int stride,index_pair_t* base_to_event_map) 
{
    while(start != stop) {

        int ei = base_to_event_map[start].start;
        if(ei != -1)
            return ei;
        start += stride;
    }
    return -1;
}

//from https://github.com/hasindu2008/f5c/blob/64d9d51cd773e46e1a72dbf2e2bf333be24df7c4/meth.c
static inline int get_closest_event_to(int k_idx, index_pair_t* base_to_event_map, int base_to_event_map_size)
{
    int stop_before = std::max(0, k_idx - 1000);
    int stop_after = std::min(k_idx + 1000, base_to_event_map_size - 1);

    int event_before = get_next_event(k_idx, stop_before, -1, base_to_event_map);
    int event_after = get_next_event(k_idx, stop_after, 1, base_to_event_map);

    // TODO: better selection of "best" event to return
    if(event_before == -1)
        return event_after;
    return event_before;
}

//from https://github.com/hasindu2008/f5c/blob/64d9d51cd773e46e1a72dbf2e2bf333be24df7c4/meth.c
//Returns true if c is a valid ambiguity code
static inline bool isAmbiguous(char c) {
    switch (c) {
        case 'M':
        case 'R':
        case 'W':
        case 'S':
        case 'Y':
        case 'K':
        case 'V':
        case 'H':
        case 'D':
        case 'B':
        case 'N':
            return true;
        default:
            return false;
    }
}
//from https://github.com/hasindu2008/f5c/blob/64d9d51cd773e46e1a72dbf2e2bf333be24df7c4/meth.c
//IUPAC alphabet
static inline bool isUnambiguous(char c) {
    switch (c) {
        case 'A':
        case 'C':
        case 'G':
        case 'T':
            return true;
        default:
            return false;
    }
}

//from https://github.com/hasindu2008/f5c/blob/64d9d51cd773e46e1a72dbf2e2bf333be24df7c4/meth.c
//Returns true if c is a valid symbol in this alphabet
static inline bool isValid(char c) { return isUnambiguous(c) || isAmbiguous(c); }
static const char* complement_dna = "TGCA";
static const uint8_t rank_dna[256] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

//from https://github.com/hasindu2008/f5c/blob/64d9d51cd773e46e1a72dbf2e2bf333be24df7c4/meth.c
//reverse-complement a DNA string
static inline std::string reverse_complement(const std::string& str) {
    std::string out(str.length(), 'A');
    size_t i = 0;             // input
    int j = str.length() - 1; // output
    while (i < str.length()) {
        // complement a single base
        assert(str[i] != METHYLATED_SYMBOL);
        out[j--] = complement_dna[rank_dna[(int)str[i++]]];
    }
    return out;
}
//from https://github.com/hasindu2008/f5c/blob/64d9d51cd773e46e1a72dbf2e2bf333be24df7c4/meth.c
static inline std::string getPossibleSymbols(char c) {
    switch (c) {
        case 'A':
            return "A";
        case 'C':
            return "C";
        case 'G':
            return "G";
        case 'T':
            return "T";
        case 'M':
            return "AC";
        case 'R':
            return "AG";
        case 'W':
            return "AT";
        case 'S':
            return "CG";
        case 'Y':
            return "CT";
        case 'K':
            return "GT";
        case 'V':
            return "ACG";
        case 'H':
            return "ACT";
        case 'D':
            return "AGT";
        case 'B':
            return "CGT";
        case 'N':
            return "ACGT";
        default:
            return "";
    }
}
//from https://github.com/hasindu2008/f5c/blob/64d9d51cd773e46e1a72dbf2e2bf333be24df7c4/meth.c
static inline std::string disambiguate(const std::string& str) {
    // create output and convert lower case to upper case
    std::string out(str);
    std::transform(out.begin(), out.end(), out.begin(), ::toupper);

    size_t i = 0;
    while (i < out.length()) {
        size_t stride = 1;
        //bool is_recognition_site = false;

        assert(isValid(out[i]));
        out[i] = getPossibleSymbols(out[i])[0];
        stride = 1;


        i += stride;
    }
    return out;
}

// from https://github.com/hasindu2008/f5c/blob/64d9d51cd773e46e1a72dbf2e2bf333be24df7c4/meth.c
static std::vector<AlignedSegment> get_aligned_segments_two_params(const bam1_t* record, int read_stride)
{
    std::vector<AlignedSegment> out;
    // Initialize first segment
    out.push_back(AlignedSegment());

    // This code is derived from bam_fillmd1_core
    //uint8_t *ref = NULL;
    //uint8_t *seq = bam_get_seq(record);
    uint32_t *cigar = bam_get_cigar(record);
    const bam1_core_t *c = &record->core;

    // read pos is an index into the original sequence that is present in the FASTQ
    // on the strand matching the reference
    int read_pos = 0;

    // query pos is an index in the query string that is recorded in the bam
    // we record this as a sanity check
    //int query_pos = 0;
    
    int ref_pos = c->pos;

    for (uint32_t ci = 0; ci < c->n_cigar; ++ci) {
        
        int cigar_len = cigar[ci] >> 4;
        int cigar_op = cigar[ci] & 0xf;

        // Set the amount that the ref/read positions should be incremented
        // based on the cigar operation
        int read_inc = 0;
        int ref_inc = 0;
 
        // Process match between the read and the reference
        bool is_aligned = false;
        if(cigar_op == BAM_CMATCH || cigar_op == BAM_CEQUAL || cigar_op == BAM_CDIFF) {
            is_aligned = true;
            read_inc = read_stride;
            ref_inc = 1;
        } else if(cigar_op == BAM_CDEL) {
            ref_inc = 1;   
        } else if(cigar_op == BAM_CREF_SKIP) {
            // end the current segment and start a new one
            out.push_back(AlignedSegment());
            ref_inc = 1;
        } else if(cigar_op == BAM_CINS) {
            read_inc = read_stride;
        } else if(cigar_op == BAM_CSOFT_CLIP) {
            read_inc = 1; // special case, do not use read_stride
        } else if(cigar_op == BAM_CHARD_CLIP) {
            read_inc = 0;
        } else {
            printf("Cigar: %d\n", cigar_op);
            assert(false && "Unhandled cigar operation");
        }

        // Iterate over the pairs of aligned bases
        for(int j = 0; j < cigar_len; ++j) {
            if(is_aligned) {
                out.back().push_back({ref_pos, read_pos});
            }

            // increment
            read_pos += read_inc;
            ref_pos += ref_inc;
        }
    }
    return out;
}


// from https://github.com/hasindu2008/f5c/blob/64d9d51cd773e46e1a72dbf2e2bf333be24df7c4/meth.c
static inline int32_t flip_k_strand(int32_t read_length,int32_t k_idx, uint32_t k)
{
    return read_length - k_idx - k;
}

// from https://github.com/hasindu2008/f5c/blob/64d9d51cd773e46e1a72dbf2e2bf333be24df7c4/meth.c
struct EventAlignment
{
    // ref data
    // std::string ref_name;
    std::string ref_kmer;
    int ref_position;

    // event data
    size_t read_idx;
    // int strand_idx;
    int event_idx;
    bool rc;

    // hmm data
    std::string model_kmer;
    char hmm_state;
};








//-------------
struct EventAlignmentParameters
{
    EventAlignmentParameters()
    {
        //sr = NULL;
        et = NULL;
        fai = NULL;
        hdr = NULL;
        record = NULL;
        model = NULL;
        //strand_idx = NUM_STRANDS;
        
        alphabet = "";
        read_idx = -1;
        region_start = -1;
        region_end = -1;
    }


    // Mandatory
    //SquiggleRead* sr;
    //necessary parameters inside sr
    index_pair_t* base_to_event_map;
    int32_t read_length;


    std::string ref_name;

    ///
    const event_table* et;
    const faidx_t* fai;
    const bam_hdr_t* hdr;
    const bam1_t* record;
    model_t* model;
    //size_t strand_idx;
    
    // optional
    std::string alphabet;
    int read_idx;
    int region_start;
    int region_end;
};



 std::vector<event_alignment_t> align_read_to_ref(const EventAlignmentParameters& params, char *ref)
{


    // Sanity check input parameters
    //assert(params.sr != NULL);
    assert(params.et != NULL);
    assert(params.fai != NULL);
    assert(params.hdr != NULL);
    assert(params.record != NULL);
    assert(params.model != NULL);
    //assert(params.strand_idx < NUM_STRANDS);
    assert( (params.region_start == -1 && params.region_end == -1) || (params.region_start <= params.region_end));
    //const PoreModel* pore_model = params.get_model(); // --hasindu this is model now

    std::vector<event_alignment_t> alignment_output;
    

    // Extract the reference subsequence for the entire alignment
    // int fetched_len = 0;
    int ref_offset = params.record->core.pos;
    // std::string ref_name(params.hdr->target_name[params.record->core.tid]);
    // std::string ref_seq = get_reference_region_ts(params.fai, ref_name.c_str(), ref_offset, 
    //                                               bam_endpos(params.record), &fetched_len);
    // hasindu - a hack to get the reference sequence
    std::string ref_seq = ref;

    // convert to upper case
    std::transform(ref_seq.begin(), ref_seq.end(), ref_seq.begin(), ::toupper);
    
    // k from read pore model
    const uint32_t k = KMER_SIZE;

    // If the reference sequence contains ambiguity codes
    // switch them to the lexicographically lowest base
    ref_seq =disambiguate(ref_seq);
    std::string rc_ref_seq = reverse_complement(ref_seq);

    // --hasindu : this is already done outside of this function
    // Skip unmapped
    // if((params.record->core.flag & BAM_FUNMAP) != 0) {
    //     return alignment_output;
    // }

    // Get the read-to-reference aligned segments
    std::vector<AlignedSegment> aligned_segments = get_aligned_segments_two_params(params.record,1);
    for(size_t segment_idx = 0; segment_idx < aligned_segments.size(); ++segment_idx) {

        AlignedSegment& aligned_pairs = aligned_segments[segment_idx];

        if(params.region_start != -1 && params.region_end != -1) {
            trim_aligned_pairs_to_ref_region(aligned_pairs, params.region_start, params.region_end);
        }

        // Trim the aligned pairs to be within the range of the maximum kmer index
        int max_kmer_idx = params.read_length - k;
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
            read_kidx_start = flip_k_strand(params.read_length,read_kidx_start, k);
            read_kidx_end = flip_k_strand(params.read_length,read_kidx_end, k);
        }
        
        assert(read_kidx_start >= 0);
        assert(read_kidx_end >= 0);

        int first_event = get_closest_event_to(read_kidx_start, params.base_to_event_map, params.read_length-KMER_SIZE + 1);
        int last_event = get_closest_event_to(read_kidx_end,params.base_to_event_map, params.read_length-KMER_SIZE + 1);
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
                curr_end_read = flip_k_strand(params.read_length,curr_end_read, k);
            }
            assert(curr_end_read >= 0);

            int s = curr_start_ref - ref_offset;
            int l = curr_end_ref - curr_start_ref + 1;

            std::string fwd_subseq = ref_seq.substr(s, l);
            std::string rc_subseq = rc_ref_seq.substr(ref_seq.length() - s - l, l);
            assert(fwd_subseq.length() == rc_subseq.length());

     //?       HMMInputSequence hmm_sequence(fwd_subseq, rc_subseq, pore_model->pmalphabet);
            
            // Require a minimum amount of sequence to align to
            if(fwd_subseq.length() < 2 * k) {
                break;
            }


            //printf("[SEGMENT_START] read: %s event start: %zu event end: %zu\n", params.sr->read_name.c_str(), input.event_start_idx, input.event_stop_idx);

            // A limitation of the segment-by-segment alignment is that we can't jump
            // over very large deletions wrt to the reference. The effect of this
            // is that we can get segments that have very few alignable events. We
            // just stop processing them for now
            int input_event_stop_idx = get_closest_event_to(curr_end_read, params.base_to_event_map, params.read_length-KMER_SIZE + 1);
            if(abs((int)curr_start_event - input_event_stop_idx < 2))
                break;



            int8_t event_stride = curr_start_event < input_event_stop_idx ? 1 : -1;  
            std::vector<HMMAlignmentState> event_alignment = profile_hmm_align(fwd_subseq, rc_subseq, k,curr_start_event,input_event_stop_idx,event_stride);
            
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

                    event_alignment_t ea;
                    
                    
                    ea.ref_position = curr_start_ref + as.kmer_idx;
                    std::string ref__kmer = ref_seq.substr(ea.ref_position - ref_offset, k);
                    strcpy(ea.ref_kmer, ref__kmer.c_str()); 

                    // event
                    ea.read_idx = params.read_idx;
                    // ea.strand_idx = params.strand_idx;
                    ea.event_idx = as.event_idx;
                    ea.rc = rc_flags[0];

                    // hmm
                    ea.hmm_state = as.state;

                    if(ea.hmm_state != 'B') {
                        // hiruna
                        // since we are using one strand, the rc flag should be rc_flags[0] and removed calling get_kmer and replaced it with code inside
                        // from https://github.com/jts/nanopolish/blob/b9dc627e73816a415e4b96b14e8a8bf53622a6c9/src/hmm/nanopolish_hmm_input_sequence.h
                        // ea.model_kmer = ! rc_flags[0] ? fwd_subseq.substr(as.kmer_idx, k) : rc_subseq.substr(rc_subseq.length() - as.kmer_idx - k, k);
                        if(rc_flags[0]){
                            std::string kmer_one = rc_subseq.substr(rc_subseq.length() - as.kmer_idx - k, k);
                            strcpy(ea.ref_kmer, kmer_one.c_str());
                        }
                        else{
                            std::string kmer_one = fwd_subseq.substr(as.kmer_idx, k);
                            strcpy(ea.ref_kmer, kmer_one.c_str());
                        }
                    } else {
                        strcpy(ea.ref_kmer, std::string(k, 'N').c_str());
                        // ea.model_kmer = std::string(k, 'N');
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

//hiruna
// right now we don't use this to train
// #if EVENTALIGN_TRAIN
//             // update training data for read
//             params.sr->parameters[params.strand_idx].add_training_from_alignment(hmm_sequence, input, event_alignment);
//             global_training[params.strand_idx].add_training_from_alignment(hmm_sequence, input, event_alignment);
// #endif

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
                  int region_end, event_table* events, model_t* model)
{
    // Load a squiggle read for the mapped read
    std::string read_name = bam_get_qname(record);

    // load read -- hasindu : we have to get rid of this sr
    //SquiggleRead sr(read_name, read_db, opt::write_samples ? SRF_LOAD_RAW_SAMPLES : 0);

    //hiruna 
    //commented this
    // if(opt::verbose > 1) {
    //     fprintf(stderr, "Realigning %s [%zu]\n",
    //             read_name.c_str(), events->n);
    // }

    //for(int strand_idx = 0; strand_idx < 2; ++strand_idx) { -- hasindu : only 1 stand in our case

        // Do not align this strand if it was not sequenced
        // if(!sr.has_events_for_strand(strand_idx)) {
        //     continue;
        // }

        EventAlignmentParameters params;
        //params.sr = &sr; -- hasindu : we have to get rid of this sr
        params.et = events; // hasindu : what is required inside the sr is to be added like this
        params.model = model; // hasindu : what is required inside the sr is to be added like this
        params.fai = fai;
        params.hdr = hdr;
        params.record = record;
       // params.strand_idx = strand_idx; -- hasindu : only 1 stand in our case
        
        params.read_idx = read_idx;
        params.region_start = region_start;
        params.region_end = region_end;

        std::vector<event_alignment_t> alignment = align_read_to_ref(params,ref);

        //-- hasindu : output writing will be done outside of this function
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

eventalign.c
Displaying eventalign.c.
