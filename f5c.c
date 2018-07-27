#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "f5c.h"
#include "f5cmisc.h"

#define m_min_mapping_quality 30

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

    //reference file
    core->fai = fai_load(fastafile);
    NULL_CHK(core->fai);

    //readbb
    core->readbb = new ReadDB;
    core->readbb->load(fastqfile);

    core->opt = opt;
    return core;
}

void free_core(core_t* core) {
    delete core->readbb;
    fai_destroy(core->fai);
    sam_close(core->m_bam_fh);
    free(core);
}

db_t* init_db() {
    db_t* db = (db_t*)(malloc(sizeof(db_t)));
    MALLOC_CHK(db);

    db->capacity_bam_rec = 512;
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
    db->f5 = (fast5_t**)malloc(sizeof(fast5_t*) * db->capacity_bam_rec);
    MALLOC_CHK(db->f5);

    return db;
}

int32_t load_db(core_t* core, db_t* db) {
    //get bams
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
                    m_min_mapping_quality) { //remove secondraies?
                //printf("%s\t%d\n",bam_get_qname(db->bam_rec[db->n_bam_rec]),result);
                db->n_bam_rec++;
            }
        }
    }
    //fprintf(stderr,"%s:: %d queries read\n",__func__,db->n_bam_rec);

    //get ref sequences (can make efficient by taking the the start and end of the sorted bam)
    for (i = 0; i < db->n_bam_rec; i++) {
        bam1_t* record = db->bam_rec[i];
        char* ref_name = core->m_hdr->target_name[record->core.tid];
        //printf("refname : %s\n",ref_name);
        int32_t ref_start_pos = record->core.pos;
        int32_t ref_end_pos = bam_endpos(record);
        assert(ref_end_pos >= ref_start_pos);

        // Extract the reference sequence for this region
        int32_t fetched_len = 0;
        char* refseq =
            faidx_fetch_seq(core->fai, ref_name, ref_start_pos, ref_end_pos,
                            &fetched_len); //error handle?
        db->fasta_cache[i] = refseq;
        // printf("seq : %s\n",db->fasta_cache[i]);

        //get the fast5

        // Get the read type from the fast5 file
        std::string qname = bam_get_qname(db->bam_rec[i]);
        char* fast5_path =
            (char*)malloc(core->readbb->get_signal_path(qname).size() +
                          10); //is +10 needed? do errorcheck
        strcpy(fast5_path, core->readbb->get_signal_path(qname).c_str());

        hid_t hdf5_file = fast5_open(fast5_path);
        if (hdf5_file >= 0) {
            db->f5[i] = (fast5_t*)calloc(1, sizeof(fast5_t));
            MALLOC_CHK(db->f5[i]);
            fast5_read(hdf5_file, db->f5[i]); //todo : errorhandle
            fast5_close(hdf5_file);
        } else {
            WARNING("Fast5 file is unreadable and will be skipped: %s",
                    fast5_path);
        }

        if (core->opt.print_raw) {
            printf("@%s\t%s\t%llu\n", qname.c_str(), fast5_path,
                   db->f5[i]->nsample);
            uint32_t j = 0;
            for (j = 0; j < db->f5[i]->nsample; j++) {
                printf("%d\t", (int)db->f5[i]->rawptr[j]);
            }
            printf("\n");
        }
    }
    //fprintf(stderr,"%s:: %d fast5 read\n",__func__,db->n_bam_rec);

    return db->n_bam_rec;
}

//todo :
// destroy_databatch();

void* process_db(core_t* core, db_t* db) {
    event_table* et = (event_table*)malloc(sizeof(event_table) * db->n_bam_rec);
    MALLOC_CHK(et);

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
        et[i] = getevents(db->f5[i]->nsample, rawptr);
    }

    return (void*)et;
}

void free_db(db_t* db) {

}

void init_opt(opt_t* opt) {
    memset(opt, 0, sizeof(opt_t));

    opt->print_raw = 0;
    opt->min_mapq = 30;
    opt->con_sec = 0;
}
