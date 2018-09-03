#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "f5c.h"
#include "f5cmisc.h"

core_t* init_core(const char* bamfilename, const char* fastafile,
                  const char* fastqfile, opt_t opt) {
    core_t* core = (core_t*)malloc(sizeof(core_t));
    MALLOC_CHK(core);

    // load bam file
    core->m_bam_fh = sam_open(bamfilename, "r");
    NULL_CHK(core->m_bam_fh);

    // load bam index file
    core->m_bam_idx = sam_index_load(core->m_bam_fh, bamfilename);
    NULL_CHK(core->m_bam_idx);

    // read the bam header
    core->m_hdr = sam_hdr_read(core->m_bam_fh);
    NULL_CHK(core->m_hdr);

    core->itr = sam_itr_queryi(core->m_bam_idx, HTS_IDX_START, 0, 0);
    NULL_CHK(core->itr);

    // reference file
    core->fai = fai_load(fastafile);
    NULL_CHK(core->fai);

    // readbb
    core->readbb = new ReadDB;
    core->readbb->load(fastqfile);

    //model
    core->model = (model_t*)malloc(
        sizeof(model_t) * NUM_KMER); //4096 is 4^6 which os hardcoded now
    MALLOC_CHK(core->model);
    core->cpgmodel = (model_t*)malloc(
        sizeof(model_t) * NUM_KMER); //4096 is 4^6 which os hardcoded now
    MALLOC_CHK(core->cpgmodel);

    //load the model from files
    if (opt.model_file) {
        read_model(core->model, opt.model_file);
    } else {
        set_model(core->model);
    }

    //todo : load the cpg model from file
    set_cpgmodel(core->cpgmodel);

    core->opt = opt;
    return core;
}

void free_core(core_t* core) {
    free(core->model);
    free(core->cpgmodel);
    delete core->readbb;
    fai_destroy(core->fai);
    sam_itr_destroy(core->itr);
    bam_hdr_destroy(core->m_hdr);
    hts_idx_destroy(core->m_bam_idx);
    sam_close(core->m_bam_fh);
    free(core);
}

db_t* init_db(core_t* core) {
    db_t* db = (db_t*)(malloc(sizeof(db_t)));
    MALLOC_CHK(db);

    db->capacity_bam_rec = core->opt.batch_size;
    db->n_bam_rec = 0;

    db->bam_rec = (bam1_t**)(malloc(sizeof(bam1_t*) * db->capacity_bam_rec));
    MALLOC_CHK(db->bam_rec);

    int32_t i = 0;
    for (i = 0; i < db->capacity_bam_rec; ++i) {
        db->bam_rec[i] = bam_init1();
        NULL_CHK(db->bam_rec[i]);
    }

    db->fasta_cache = (char**)(malloc(sizeof(char*) * db->capacity_bam_rec));
    MALLOC_CHK(db->fasta_cache);
    db->read = (char**)(malloc(sizeof(char*) * db->capacity_bam_rec));

    db->f5 = (fast5_t**)malloc(sizeof(fast5_t*) * db->capacity_bam_rec);
    MALLOC_CHK(db->f5);

    db->et = (event_table*)malloc(sizeof(event_table) * db->capacity_bam_rec);
    MALLOC_CHK(db->et);

    return db;
}

int32_t load_db(core_t* core, db_t* db) {
    // get bams
    bam1_t* record;
    int32_t result = 0;
    db->n_bam_rec = 0;
    int32_t i = 0;
    while (db->n_bam_rec < db->capacity_bam_rec) {
        record = db->bam_rec[db->n_bam_rec];
        result = sam_itr_next(core->m_bam_fh, core->itr, record);

        if (result < 0) {
            break;
        } else {
            if ((record->core.flag & BAM_FUNMAP) == 0 &&
                record->core.qual >=
                    core->opt
                        .min_mapq) { // remove secondraies? //need to use the user parameter
                // printf("%s\t%d\n",bam_get_qname(db->bam_rec[db->n_bam_rec]),result);
                db->n_bam_rec++;
            }
        }
    }
    // fprintf(stderr,"%s:: %d queries read\n",__func__,db->n_bam_rec);

    // get ref sequences (can make efficient by taking the the start and end of
    // the sorted bam)
    for (i = 0; i < db->n_bam_rec; i++) {
        bam1_t* record = db->bam_rec[i];
        char* ref_name = core->m_hdr->target_name[record->core.tid];
        // printf("refname : %s\n",ref_name);
        int32_t ref_start_pos = record->core.pos;
        int32_t ref_end_pos = bam_endpos(record);
        assert(ref_end_pos >= ref_start_pos);

        // Extract the reference sequence for this region
        int32_t fetched_len = 0;
        char* refseq =
            faidx_fetch_seq(core->fai, ref_name, ref_start_pos, ref_end_pos,
                            &fetched_len); // error handle?
        db->fasta_cache[i] = refseq;
        // printf("seq : %s\n",db->fasta_cache[i]);

        // get the fast5

        // Get the read type from the fast5 file
        std::string qname = bam_get_qname(db->bam_rec[i]);
        char* fast5_path =
            (char*)malloc(core->readbb->get_signal_path(qname).size() +
                          10); // is +10 needed? do errorcheck
        strcpy(fast5_path, core->readbb->get_signal_path(qname).c_str());

        hid_t hdf5_file = fast5_open(fast5_path);
        if (hdf5_file >= 0) {
            db->f5[i] = (fast5_t*)calloc(1, sizeof(fast5_t));
            MALLOC_CHK(db->f5[i]);
            fast5_read(hdf5_file, db->f5[i]); // todo : errorhandle
            fast5_close(hdf5_file);
        } else {
            if (core->opt.flag & F5C_SKIP_UNREADABLE) {
                WARNING("Fast5 file is unreadable and will be skipped: %s",
                        fast5_path);
            } else {
                ERROR("Fast5 file could not be opened: %s", fast5_path);
                exit(EXIT_FAILURE);
            }
        }

        if (core->opt.flag & F5C_PRINT_RAW) {
            printf(">%s\tPATH:%s\tLN:%llu\n", qname.c_str(), fast5_path,
                   db->f5[i]->nsample);
            uint32_t j = 0;
            for (j = 0; j < db->f5[i]->nsample; j++) {
                printf("%d\t", (int)db->f5[i]->rawptr[j]);
            }
            printf("\n");
        }

        free(fast5_path);

        //get the read in ascci
        db->read[i] =
            (char*)malloc(core->readbb->get_read_sequence(qname).size() +
                          1); // is +1 needed? do errorcheck
        strcpy(db->read[i], core->readbb->get_read_sequence(qname).c_str());
    }
    // fprintf(stderr,"%s:: %d fast5 read\n",__func__,db->n_bam_rec);

    return db->n_bam_rec;
}

void process_db(core_t* core, db_t* db) {
    int32_t i;
    for (i = 0; i < db->n_bam_rec; i++) {
        float* rawptr = db->f5[i]->rawptr;
        float range = db->f5[i]->range;
        float digitisation = db->f5[i]->digitisation;
        float offset = db->f5[i]->offset;
        int32_t nsample = db->f5[i]->nsample;

        // convert to pA
        float raw_unit = range / digitisation;
        for (int32_t j = 0; j < nsample; j++) {
            rawptr[j] = (rawptr[j] + offset) * raw_unit;
        }
        db->et[i] = getevents(db->f5[i]->nsample, rawptr);

        //have to test if the computed events are correct
        //get the scalings
        scalings_t scalings =
            estimate_scalings_using_mom(db->read[i], core->model, db->et[i]);

        int n_events = 0;
        AlignedPair* event_alignment =
            align(db->read[i], db->et[i], core->model, scalings,
                  db->f5[i]->sample_rate, &n_events);
        assert(n_events > 0);
        //int32_t n_events=100;
        if (event_alignment != NULL) {
            event_alignment_t* alignment_output =
                postalign(db->read[i], event_alignment, n_events);
            recalibrate_model(core->model, db->et[i], &scalings,
                              alignment_output, n_events, 1);
            free(event_alignment);
        }
    }

    return;
}

void output_db(core_t* core, db_t* db) {
    if (core->opt.flag & F5C_PRINT_EVENTS) {
        int32_t i = 0;
        for (i = 0; i < db->n_bam_rec; i++) {
            printf(">%s\tLN:%d\tSTART:%d\tEND:%d\n",
                   bam_get_qname(db->bam_rec[i]), (int)db->et[i].n,
                   (int)db->et[i].start, (int)db->et[i].end);
            uint32_t j = 0;
            for (j = 0; j < db->et[i].n; j++) {
                printf("{%d,%f,%f,%f,%d,%d}\t", (int)db->et[i].event[j].start,
                       db->et[i].event[j].length, db->et[i].event[j].mean,
                       db->et[i].event[j].stdv, (int)db->et[i].event[j].pos,
                       (int)db->et[i].event[j].state);
            }
            printf("\n");
        }
    }
}

void free_db_tmp(db_t* db) {
    int32_t i = 0;
    for (i = 0; i < db->n_bam_rec; ++i) {
        bam_destroy1(db->bam_rec[i]);
        db->bam_rec[i] = bam_init1();
        free(db->fasta_cache[i]);
        free(db->read[i]);
        free(db->f5[i]->rawptr);
        free(db->f5[i]);
        free(db->et[i].event);
    }
}

void free_db(db_t* db) {
    int32_t i = 0;
    for (i = 0; i < db->capacity_bam_rec; ++i) {
        bam_destroy1(db->bam_rec[i]);
    }
    free(db->bam_rec);
    free(db->fasta_cache);
    free(db->read);
    free(db->et);
    free(db->f5);
    free(db);
}

void init_opt(opt_t* opt) {
    memset(opt, 0, sizeof(opt_t));
    opt->min_mapq = 30;
    opt->batch_size = 500;
}
