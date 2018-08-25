#include "f5c.h"
#include <assert.h>
#include <vector>

// reverse-complement a string
// when the string contains methylated bases, the methylation
// symbol transfered to the output strand in the appropriate position
virtual std::string reverse_complement(const std::string& str) const {
    std::string out(str.length(), 'A');
    size_t i = 0;             // input
    int j = str.length() - 1; // output
    while (i < str.length()) {
        int recognition_index = -1;
        RecognitionMatch match;

        // Does this location (partially) match a methylated recognition site?
        for (size_t j = 0; j < num_recognition_sites(); ++j) {
            match = match_to_site(str, i, get_recognition_site_methylated(j),
                                  recognition_length());
            if (match.length > 0 && match.covers_methylated_site) {
                recognition_index = j;
                break;
            }
        }

        // If this subsequence matched a methylated recognition site,
        // copy the complement of the site to the output
        if (recognition_index != -1) {
            for (size_t k = match.offset; k < match.offset + match.length;
                 ++k) {
                out[j--] = get_recognition_site_methylated_complement(
                    recognition_index)[k];
                i += 1;
            }
        } else {
            // complement a single base
            assert(str[i] != METHYLATED_SYMBOL);
            out[j--] = complement(str[i++]);
        }
    }
    return out;
}

virtual std::string disambiguate(const std::string& str) const {
    // create output and convert lower case to upper case
    std::string out(str);
    std::transform(out.begin(), out.end(), out.begin(), ::toupper);

    size_t i = 0;
    while (i < out.length()) {
        size_t stride = 1;
        bool is_recognition_site = false;

        // Does this location (partially) match a methylated recognition site?
        // not to worry about in DNA alphabet
        // for(size_t j = 0; j < num_recognition_sites(); ++j) {

        //     RecognitionMatch match = match_to_site(out, i, get_recognition_site_methylated(j), recognition_length());
        //     if(match.length > 0) {
        //         stride = match.length; // skip to end of match
        //         is_recognition_site = true;
        //         break;
        //     }
        // }

        // disambiguate if not a recognition site
        //if(!is_recognition_site) {
        assert(IUPAC::isValid(out[i]));
        out[i] = IUPAC::getPossibleSymbols(out[i])[0];
        stride = 1;
        //}

        i += stride;
    }
    return out;
}

// Test CpG sites in this read for methylation
void calculate_methylation_for_read(char* ref, const OutputHandles& handles,
                                    const ReadDB& read_db, const faidx_t* fai,
                                    const bam_hdr_t* hdr, const bam1_t* record,
                                    size_t read_idx, int region_start,
                                    int region_end) {
    // Load a squiggle read for the mapped read
    // std::string read_name = bam_get_qname(record);
    // SquiggleRead sr(read_name, read_db);

    // An output map from reference positions to scored CpG sites
    std::map<int, ScoredSite> site_score_map;

    //todo : do this
    // if(!sr.has_events_for_strand(strand_idx)) {
    //     continue;
    // }

    //size_t k = sr.get_model_k(strand_idx);
    uint32_t k = KMER_SIZE;

    // // check if there is a cpg model for this strand
    // if(!PoreModelSet::has_model(sr.get_model_kit_name(strand_idx),
    //                             "cpg",
    //                             sr.get_model_strand_name(strand_idx),
    //                             k))
    // {
    //     continue;
    // }

    // Build the event-to-reference map for this read from the bam record
    //SequenceAlignmentRecord seq_align_record(record);
    //EventAlignmentRecord event_align_record(&sr, strand_idx, seq_align_record);

    std::vector<double> site_scores;
    std::vector<int> site_starts;
    std::vector<int> site_ends;
    std::vector<int> site_count;

    // std::string contig = hdr->target_name[record->core.tid];
    // int ref_start_pos = record->core.pos;
    // int ref_end_pos =  bam_endpos(record);

    // // Extract the reference sequence for this region
    // int fetched_len = 0;
    // assert(ref_end_pos >= ref_start_pos);
    // std::string ref_seq = get_reference_region_ts(fai, contig.c_str(), ref_start_pos,
    //                                               ref_end_pos, &fetched_len);

    // Remove non-ACGT bases from this reference segment
    std::string ref_seq = ref; //todo : this is a hack due to pure laziness
    ref_seq = disambiguate(ref_seq); //todo : do this

    // Scan the sequence for CpGs
    std::vector<int> cpg_sites;
    assert(ref_seq.size() != 0);
    for (size_t i = 0; i < ref_seq.size() - 1; ++i) {
        if (ref_seq[i] == 'C' && ref_seq[i + 1] == 'G') {
            cpg_sites.push_back(i);
        }
    }

    // Batch the CpGs together into groups that are separated by some minimum distance
    int min_separation = 10;
    std::vector<std::pair<int, int>> groups;

    size_t curr_idx = 0;
    while (curr_idx < cpg_sites.size()) {
        // Find the endpoint of this group of sites
        size_t end_idx = curr_idx + 1;
        while (end_idx < cpg_sites.size()) {
            if (cpg_sites[end_idx] - cpg_sites[end_idx - 1] > min_separation)
                break;
            end_idx += 1;
        }
        groups.push_back(std::make_pair(curr_idx, end_idx));
        curr_idx = end_idx;
    }

    for (size_t group_idx = 0; group_idx < groups.size(); ++group_idx) {
        size_t start_idx = groups[group_idx].first;
        size_t end_idx = groups[group_idx].second;

        // the coordinates on the reference substring for this group of sites
        int sub_start_pos = cpg_sites[start_idx] - min_separation;
        int sub_end_pos = cpg_sites[end_idx - 1] + min_separation;
        int span = cpg_sites[end_idx - 1] - cpg_sites[start_idx];

        // skip if too close to the start of the read alignment or
        // if the reference range is too large to efficiently call
        if (sub_start_pos <= min_separation || span > 200) {
            continue;
        }

        //todo : do this
        std::string subseq =
            ref_seq.substr(sub_start_pos, sub_end_pos - sub_start_pos + 1);
        std::string rc_subseq = reverse_complement(subseq);

        int calling_start = sub_start_pos + ref_start_pos;
        int calling_end = sub_end_pos + ref_start_pos;

        // using the reference-to-event map, look up the event indices for this segment
        int e1, e2;
        bool bounded = AlignmentDB::_find_by_ref_bounds(
            event_align_record.aligned_events, calling_start, calling_end, e1,
            e2);

        double ratio = fabs(e2 - e1) / (calling_start - calling_end);

        // Only process this region if the the read is aligned within the boundaries
        // and the span between the start/end is not unusually short
        if (!bounded || abs(e2 - e1) <= 10 || ratio > MAX_EVENT_TO_BP_RATIO) {
            continue;
        }

        uint32_t hmm_flags = HAF_ALLOW_PRE_CLIP | HAF_ALLOW_POST_CLIP;

        // Set up event data
        HMMInputData data;
        data.read = &sr;
        data.pore_model = sr.get_model(strand_idx, "cpg");
        data.strand = strand_idx;
        data.rc = event_align_record.rc;
        data.event_start_idx = e1;
        data.event_stop_idx = e2;
        data.event_stride =
            data.event_start_idx <= data.event_stop_idx ? 1 : -1;

        // Calculate the likelihood of the unmethylated sequence
        HMMInputSequence unmethylated(subseq, rc_subseq, mtest_alphabet);
        double unmethylated_score =
            profile_hmm_score(unmethylated, data, hmm_flags);

        // Methylate all CpGs in the sequence and score again
        std::string mcpg_subseq = mtest_alphabet->methylate(subseq);
        std::string rc_mcpg_subseq =
            mtest_alphabet->reverse_complement(mcpg_subseq);

        // Calculate the likelihood of the methylated sequence
        HMMInputSequence methylated(mcpg_subseq, rc_mcpg_subseq,
                                    mtest_alphabet);
        double methylated_score =
            profile_hmm_score(methylated, data, hmm_flags);

        // Aggregate score
        int start_position = cpg_sites[start_idx] + ref_start_pos;
        auto iter = site_score_map.find(start_position);
        if (iter == site_score_map.end()) {
            // insert new score into the map
            ScoredSite ss;
            ss.chromosome = contig;
            ss.start_position = start_position;
            ss.end_position = cpg_sites[end_idx - 1] + ref_start_pos;
            ss.n_cpg = end_idx - start_idx;

            // extract the CpG site(s) with a k-mers worth of surrounding context
            size_t site_output_start = cpg_sites[start_idx] - k + 1;
            size_t site_output_end = cpg_sites[end_idx - 1] + k;
            ss.sequence = ref_seq.substr(site_output_start,
                                         site_output_end - site_output_start);

            // insert into the map
            iter =
                site_score_map.insert(std::make_pair(start_position, ss)).first;
        }

        // set strand-specific score
        // upon output below the strand scores will be summed
        iter->second.ll_unmethylated[strand_idx] = unmethylated_score;
        iter->second.ll_methylated[strand_idx] = methylated_score;
        iter->second.strands_scored += 1;
    } // for group
}