#include <assert.h>
#include <math.h>
#include <pthread.h>
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
    MALLOC_CHK(db->read);
    db->read_len = (int32_t*)(malloc(sizeof(int32_t) * db->capacity_bam_rec));
    MALLOC_CHK(db->read_len);

    db->f5 = (fast5_t**)malloc(sizeof(fast5_t*) * db->capacity_bam_rec);
    MALLOC_CHK(db->f5);

    db->et = (event_table*)malloc(sizeof(event_table) * db->capacity_bam_rec);
    MALLOC_CHK(db->et);

    db->scalings =
        (scalings_t*)malloc(sizeof(scalings_t) * db->capacity_bam_rec);
    MALLOC_CHK(db->scalings);

    db->event_align_pairs =
        (AlignedPair**)malloc(sizeof(AlignedPair*) * db->capacity_bam_rec);
    MALLOC_CHK(db->event_align_pairs);
    db->n_event_align_pairs =
        (int32_t*)malloc(sizeof(int32_t) * db->capacity_bam_rec);
    MALLOC_CHK(db->n_event_align_pairs);

    db->event_alignment = (event_alignment_t**)malloc(
        sizeof(event_alignment_t*) * db->capacity_bam_rec);
    MALLOC_CHK(db->event_alignment);
    db->n_event_alignment =
        (int32_t*)malloc(sizeof(int32_t*) * db->capacity_bam_rec);
    MALLOC_CHK(db->n_event_alignment);

    db->events_per_base =
        (double*)malloc(sizeof(double*) * db->capacity_bam_rec);
    MALLOC_CHK(db->events_per_base);

    db->base_to_event_map =
        (index_pair_t**)malloc(sizeof(index_pair_t*) * db->capacity_bam_rec);
    MALLOC_CHK(db->base_to_event_map);

    db->read_stat_flag = (int32_t *)malloc(sizeof(int32_t) * db->capacity_bam_rec);
    MALLOC_CHK(db->read_stat_flag);

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
        db->read_len[i] = strlen(db->read[i]);

        db->read_stat_flag[i] = 0; //reset the flag
    }
    // fprintf(stderr,"%s:: %d fast5 read\n",__func__,db->n_bam_rec);

    return db->n_bam_rec;
}

void* align_pthread(void* voidargs) {
    int i;
    pthread_arg_t* args = (pthread_arg_t*)voidargs;
    db_t* db = args->db;
    core_t* core = args->core;

    for (i = args->starti; i < args->endi; i++) {
        db->n_event_align_pairs[i] = align(
            db->event_align_pairs[i], db->read[i], db->read_len[i], db->et[i],
            core->model, db->scalings[i], db->f5[i]->sample_rate);
        //fprintf(stderr,"readlen %d,n_events %d\n",db->read_len[i],n_event_align_pairs);
    }

    //fprintf(stderr,"Thread %d done\n",(myargs->position)/THREADS);
    pthread_exit(0);
}

void align_db(core_t* core, db_t* db) {
#ifdef HAVE_CUDA
    if (!(core->opt.flag & F5C_DISABLE_CUDA)) {
        fprintf(stderr, "cuda\n");
        align_cuda(core, db);
    }
#endif

    if (core->opt.flag & F5C_DISABLE_CUDA) {
        fprintf(stderr, "cpu\n");
        if (core->opt.num_thread == 1) {
            int i;
            for (i = 0; i < db->n_bam_rec; i++) {
                db->n_event_align_pairs[i] =
                    align(db->event_align_pairs[i], db->read[i],
                          db->read_len[i], db->et[i], core->model,
                          db->scalings[i], db->f5[i]->sample_rate);
                //fprintf(stderr,"readlen %d,n_events %d\n",db->read_len[i],n_event_align_pairs);
            }
        } else {
            //create threads
            pthread_t tids[core->opt.num_thread];
            pthread_arg_t pt_args[core->opt.num_thread];
            int32_t t, ret;
            int32_t i = 0;
            int32_t num_thread = core->opt.num_thread;
            int32_t step = (db->n_bam_rec + num_thread - 1) / num_thread;
            //todo : check for higher num of threads than the data
            for (t = 0; t < num_thread; t++) {
                pt_args[t].core = core;
                pt_args[t].db = db;
                pt_args[t].starti = i;
                i += step;
                if (i > db->n_bam_rec) {
                    pt_args[t].endi = db->n_bam_rec;
                } else {
                    pt_args[t].endi = i;
                }
                //fprintf(stderr,"t%d : %d-%d\n",t,pt_args[t].starti,pt_args[t].endi);
                ret = pthread_create(&tids[t], NULL, align_pthread,
                                     (void*)(&pt_args[t]));
                NEG_CHK(ret);
            }

            //pthread joining
            for (t = 0; t < core->opt.num_thread; t++) {
                int ret = pthread_join(tids[t], NULL);
                NEG_CHK(ret);
            }
        }
    }
}

void process_db(core_t* core, db_t* db, double realtime0) {
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

        //get the scalings
        db->scalings[i] = estimate_scalings_using_mom(
            db->read[i], db->read_len[i], core->model, db->et[i]);
    }
    fprintf(stderr, "[%s::%.3f*%.2f] Events computed\n", __func__,
            realtime() - realtime0, cputime() / (realtime() - realtime0));

    for (i = 0; i < db->n_bam_rec; i++) {
        db->event_align_pairs[i] = (AlignedPair*)malloc(
            sizeof(AlignedPair) * db->et[i].n *
            2); //todo : find a good heuristic to save memory
        MALLOC_CHK(db->event_align_pairs[i]);
    }

    align_db(core, db);

    fprintf(stderr, "[%s::%.3f*%.2f] Banded alignment done\n", __func__,
            realtime() - realtime0, cputime() / (realtime() - realtime0));

    for (i = 0; i < db->n_bam_rec; i++) {
        db->event_alignment[i] = NULL;
        db->n_event_alignment[i] = 0;
        db->events_per_base[i] = 0; //todo : is double needed? not just int8?

        int32_t n_kmers = db->read_len[i] - KMER_SIZE + 1;
        db->base_to_event_map[i]=(index_pair_t*)(malloc(sizeof(index_pair_t) * n_kmers));
        MALLOC_CHK(db->base_to_event_map[i]);

        if (db->n_event_align_pairs[i] > 0) {
            // prepare data structures for the final calibration
            
            db->event_alignment[i] = (event_alignment_t*)malloc(
                sizeof(event_alignment_t) * db->n_event_align_pairs[i]);
            MALLOC_CHK(db->event_alignment[i]);

            // for (int j = 0; j < n_event_align_pairs; ++j) {
            //     fprintf(stderr, "%d-%d\n",event_align_pairs[j].ref_pos,event_align_pairs[j].read_pos);
            // }


            //todo : verify if this n is needed is needed
            db->n_event_alignment[i] = postalign(
                db->event_alignment[i],db->base_to_event_map[i], &db->events_per_base[i], db->read[i],
                n_kmers, db->event_align_pairs[i], db->n_event_align_pairs[i]);

            //fprintf(stderr,"n_event_alignment %d\n",n_events);

            // run recalibration to get the best set of scaling parameters and the residual
            // between the (scaled) event levels and the model.

            // internally this function will set shift/scale/etc of the pore model
            bool calibrated = recalibrate_model(
                core->model, db->et[i], &db->scalings[i],
                db->event_alignment[i], db->n_event_alignment[i], 1);

            // QC calibration
            if (!calibrated || db->scalings[i].var > MIN_CALIBRATION_VAR) {
                //     events[strand_idx].clear();
                free(db->event_alignment[i]);
                //free(db->event_align_pairs[i]);
                //     g_failed_calibration_reads += 1; //todo : add these stats
                db->read_stat_flag[i] |= FAILED_CALIBRATION;
                continue;
            }

            free(db->event_alignment[i]);

        } else {
            // Could not align, fail this read
            // this->events[strand_idx].clear();
            // this->events_per_base[strand_idx] = 0.0f;
            //free(db->event_align_pairs[i]);
            db->read_stat_flag[i] |= FAILED_ALIGNMENT;
            // g_failed_alignment_reads += 1; //todo : add these stats
            continue;
        }

        // Filter poor quality reads that have too many "stays"

        if (db->events_per_base[i] > 5.0) {
            //     g_qc_fail_reads += 1; //todo : add these stats
            //     events[0].clear();
            //     events[1].clear();
            //free(db->event_align_pairs[i]);
            db->read_stat_flag[i] |= FAILED_QUALITY_CHK;
            continue;
        }

        calculate_methylation_for_read(db->fasta_cache[i], db->bam_rec[i], db->read_len[i], db->et[i].event, db->base_to_event_map[i],
db->scalings[i], core->cpgmodel);
    }

    return;
}

void output_db(core_t* core, db_t* db) {
    if (core->opt.flag & F5C_PRINT_EVENTS) {
        int32_t i = 0;
        for (i = 0; i < db->n_bam_rec; i++) {
            printf(">%s\tLN:%d\tEVENTSTART:%d\tEVENTEND:%d\n",
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
    if (core->opt.flag & F5C_PRINT_BANDED_ALN) {
        int32_t i = 0;
        for (i = 0; i < db->n_bam_rec; i++) {
            if((db->read_stat_flag[i]) & FAILED_ALIGNMENT){
                continue;
            }
            printf(">%s\tN_ALGN_PAIR:%d\t{ref_os,read_pos}\n",
                   bam_get_qname(db->bam_rec[i]),
                   (int)db->n_event_align_pairs[i]);
            AlignedPair* event_align_pairs = db->event_align_pairs[i];
            int32_t j = 0;
            for (j = 0; j < db->n_event_align_pairs[i]; j++) {
                printf("{%d,%d}\t", event_align_pairs[j].ref_pos,
                       event_align_pairs[j].read_pos);
            }
            printf("\n");
        }
    }

    if (core->opt.flag & F5C_PRINT_SCALING) {
        int32_t i = 0;
        printf("read\tshift\tscale\tvar\n");

        for (i = 0; i < db->n_bam_rec; i++) {
            if((db->read_stat_flag[i])&(FAILED_ALIGNMENT|FAILED_CALIBRATION)){
                continue;
            }
            printf("%s\t%.2lf\t%.2lf\t%.2lf\n", bam_get_qname(db->bam_rec[i]),
                   db->scalings[i].shift, db->scalings[i].scale,
                   db->scalings[i].var);
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
        free(db->event_align_pairs[i]);
        free(db->base_to_event_map[i]);
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
    free(db->read_len);
    free(db->et);
    free(db->f5);
    free(db->scalings);
    free(db->event_align_pairs);
    free(db->n_event_align_pairs);
    free(db->event_alignment);
    free(db->n_event_alignment);
    free(db->events_per_base);
    free(db->base_to_event_map);
    free(db->read_stat_flag);

    free(db);
}

void init_opt(opt_t* opt) {
    memset(opt, 0, sizeof(opt_t));
    opt->min_mapq = 30;
    opt->batch_size = 512;
    opt->num_thread = 8;
#ifndef HAVE_CUDA
    opt->flag |= F5C_DISABLE_CUDA;
#endif
    opt->cuda_block_size=8;    
}
