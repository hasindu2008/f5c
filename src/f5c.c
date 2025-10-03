/* @file f5c.c
**
** f5c interface implementation
** @author: Hasindu Gamaarachchi (hasindu@unsw.edu.au)
** @@
******************************************************************************/

#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "f5c.h"
#include "f5cmisc.h"

#include <sys/wait.h>
#include <unistd.h>

/*
todo :
Error counter for consecutive failures in the skip unreadable mode
not all the memory allocations are needed for eventalign mode
*/

//read bed file
static inline char **read_bed_regions(char *bedfile, int64_t *count){

    FILE *bedfp = fopen(bedfile,"r");
    F_CHK(bedfp,bedfile);

    char* buffer = (char*)malloc(sizeof(char) * (100)); //READ+newline+nullcharacter
    MALLOC_CHK(buffer);

    int64_t reg_capcacity = 1024;
    int64_t reg_i = 0;
    char **reg_list = (char **)malloc(reg_capcacity * sizeof(char *));
    MALLOC_CHK(reg_list);


    size_t bufferSize = 100;
    ssize_t readlinebytes = 0;
    int64_t line_no = 0;


    while ((readlinebytes = getline(&buffer, &bufferSize, bedfp)) != -1) {

        char *ref = (char *)malloc(sizeof(char)*readlinebytes);
        MALLOC_CHK(ref);
        int64_t beg=-1;
        int64_t end=-1;

        //TODO can optimised though strtok etc later
        int ret=sscanf(buffer,"%s\t%ld\t%ld",ref,&beg, &end);
        if(ret!=3 || end<beg){
            ERROR("Malformed bed entry at line %ld",line_no);
            exit(EXIT_FAILURE);
        }

        if(reg_i>=reg_capcacity){
            if(reg_capcacity>1000000){
                WARNING("The region bed file has over %ld regions. To reduce memory usage, you may consider merging bed regions.",reg_i);
            }
            reg_capcacity=reg_capcacity*2;
            reg_list = (char **)realloc((void *)reg_list,reg_capcacity * sizeof(char *));
            MALLOC_CHK(reg_list);

        }

        reg_list[reg_i] = (char *)malloc(sizeof(char)*readlinebytes);
        sprintf(reg_list[reg_i],"%s:%ld-%ld",ref, beg, end);
        reg_i++;


        free(ref);



        line_no++;
    }

    fclose(bedfp);
    free(buffer);
    *count = reg_i;

    return reg_list;
}


int8_t drna_detect(slow5_file_t *sp){

    const slow5_hdr_t* hdr = sp->header;
    int8_t rna = 0;
    char *exp =slow5_hdr_get("experiment_type", 0, hdr);
    if(exp==NULL){
        WARNING("%s","experiment_type not found in SLOW5 header. Assuming genomic_dna");
        return 0;
    }
    if (strcmp(exp,"genomic_dna")==0){
        rna = 0;
    }else if (strcmp(exp,"rna")==0){
        rna = 1;
    } else {
        WARNING("Unknown experiment type: %s. Assuming genomic_dna", exp);
    }

    for(uint32_t  i=1; i < hdr->num_read_groups; i++){
        char *curr =slow5_hdr_get("experiment_type", i, hdr);
        if (strcmp(curr, exp)){
            WARNING("Experiment type mismatch: %s != %s in read group %d. Defaulted to %s", curr, exp, i, exp);
        }
    }
    return rna;
}


int8_t pore_detect(slow5_file_t *sp){

    const slow5_hdr_t* hdr = sp->header;
    int8_t pore = OPT_PORE_R9;
    char *kit =slow5_hdr_get("sequencing_kit", 0, hdr);
    if(kit==NULL){
        WARNING("%s","sequencing_kit not found in SLOW5 header. Assuming R9.4.1");
        return 0;
    }
    if (strstr(kit,"114")!=NULL){
        pore = OPT_PORE_R10;
    } else if (strstr(kit,"rna004")!=NULL){
        pore = OPT_PORE_RNA004;
    } else {
        pore = OPT_PORE_R9;
    }

    for(uint32_t  i=1; i < hdr->num_read_groups; i++){
        char *curr =slow5_hdr_get("sequencing_kit", i, hdr);
        if (strcmp(curr, kit)){
            WARNING("sequencing_kit type mismatch: %s != %s in read group %d. Defaulted to %s", curr, kit, i, kit);
        }
    }
    return pore;
}

/* initialise the core data structure */
core_t* init_core(const char* bamfilename, const char* fastafile,
                  const char* fastqfile, const char* tmpfile, opt_t opt,double realtime0, int8_t mode, char *eventalignsummary, char *slow5file) {

    core_t* core = (core_t*)malloc(sizeof(core_t));
    MALLOC_CHK(core);

    if(opt.num_iop > 1){
        init_iop(core,opt);
    }

    // load bam file
    core->m_bam_fh = sam_open(bamfilename, "r");
    NULL_CHK(core->m_bam_fh);

    // load bam index file
    core->m_bam_idx = sam_index_load(core->m_bam_fh, bamfilename);
    if(core->m_bam_idx==NULL){
        ERROR("could not load the .bai index file for %s", bamfilename);
        fprintf(stderr, "Please run 'samtools index %s'\n", bamfilename);
        exit(EXIT_FAILURE);
    }

    // read the bam header
    core->m_hdr = sam_hdr_read(core->m_bam_fh);
    NULL_CHK(core->m_hdr);

    // If processing a region of the genome, get clipping coordinates
    core->clip_start = -1;
    core->clip_end = -1;

    core->reg_list=NULL; //region list is NULL by default
    core->reg_i=0;
    core->reg_n=0;
    core->itr = NULL;

    if(opt.region_str == NULL){
        core->itr = sam_itr_queryi(core->m_bam_idx, HTS_IDX_START, 0, 0);
        if(core->itr==NULL){
            ERROR("%s","sam_itr_queryi failed. A problem with the BAM index?");
            exit(EXIT_FAILURE);
        }
    }
    else{
        //determine if .bed
        int region_str_len = strlen(opt.region_str);
        if(region_str_len>=4 && strcmp(&(opt.region_str[region_str_len-4]),".bed")==0 ){
            STDERR("Fetching the list of regions from file: %s", opt.region_str);
            WARNING("%s", "When loading windows from a bed file, output is based on reads that are unclipped. Also, there may be repeated entries when regions overlap.");
            int64_t count=0;
            char **reg_l = read_bed_regions(opt.region_str, &count);
            core->reg_list = reg_l;
            core->reg_i = 0;
            core->reg_n = count;
        }
        else{
            STDERR("Iterating over region: %s\n", opt.region_str);
            core->itr = sam_itr_querys(core->m_bam_idx, core->m_hdr, opt.region_str);
            if(core->itr==NULL){
                ERROR("sam_itr_querys failed. Please check if the region string you entered [%s] is valid",opt.region_str);
                exit(EXIT_FAILURE);
            }
            hts_parse_reg(opt.region_str, &(core->clip_start) , &(core->clip_end));
        }
    }


    //open the bam file for writing skipped ultra long reads
    core->ultra_long_tmp=NULL; //todo :  at the moment this is used to detect if the load balance mode is enabled. A better method in the opt flags.
    if(tmpfile!=NULL){
        core->ultra_long_tmp = sam_open(tmpfile, "wb");
        NULL_CHK(core->ultra_long_tmp);

        //write the header to the temporary file
        int ret_sw=sam_hdr_write(core->ultra_long_tmp,core->m_hdr);
        NEG_CHK(ret_sw);
    }

    if(opt.flag & F5C_WR_RAW_DUMP){
        core->raw_dump = fopen("f5c.tmp.bin","wb");
        F_CHK(core->raw_dump,"f5c.tmp.bin");
    }
    if(opt.flag & F5C_RD_RAW_DUMP){
        core->raw_dump = fopen("f5c.tmp.bin","rb");
        F_CHK(core->raw_dump,"f5c.tmp.bin");
    }
    if(opt.flag & F5C_RD_SLOW5){
        if(opt.flag & F5C_RD_RAW_DUMP){
            STDERR("%s","Option --slow5 is incompatible with --read-dump");
            exit(EXIT_FAILURE);
        }
       if(opt.flag & F5C_WR_RAW_DUMP){
            STDERR("%s","Option --slow5 is incompatible with --write-dump");
            exit(EXIT_FAILURE);
       }
        core->sf = slow5_open(slow5file,"r");
        if (core->sf == NULL) {
            STDERR("Error opening SLOW5 file %s\n",slow5file);
            exit(EXIT_FAILURE);
        }
        int ret=slow5_idx_load(core->sf);
        if(ret<0){
            STDERR("Error in loading SLOW5 index for %s\n",slow5file);
            exit(EXIT_FAILURE);
        }
        slow5_set_skip_rid();

        if(drna_detect(core->sf)) {
            opt.flag |= F5C_RNA;
            if(opt.verbosity > 1) fprintf(stderr,"Detected RNA data. --rna was set automatically.\n");
            if(mode == 0){
                ERROR("%s","Methylation calling is not supported for RNA data.");
                exit(EXIT_FAILURE);
            }
        }


        if(opt.pore==NULL){
            int8_t pore = pore_detect(core->sf);
            if(pore){
                opt.flag |= F5C_R10;
                if(opt.verbosity > 1) {
                    if (pore == OPT_PORE_R10) fprintf(stderr,"Detected R10 data. --pore r10 was set automatically.\n");
                    if (pore == OPT_PORE_RNA004) fprintf(stderr,"Detected RNA004 data. --pore rna004 was set automatically.\n");
                    if (pore == OPT_PORE_R10 && (opt.flag & F5C_RNA)){
                        ERROR("%s","R10 RNA data does not exist! But header header indicates that the data is R10 RNA.");
                        exit(EXIT_FAILURE);
                    }
                }
            }
        }
    }

    // reference file
    core->fai = fai_load(fastafile);
    NULL_CHK(core->fai);

    // readbb
    core->readbb = new ReadDB;
    core->readbb->load(fastqfile, (opt.flag & F5C_RD_SLOW5)? 1 : 0);

    //model
    core->model = (model_t*)malloc(sizeof(model_t) * MAX_NUM_KMER);
    MALLOC_CHK(core->model);
    if(mode == 0){
        core->cpgmodel = (model_t*)malloc(sizeof(model_t) * MAX_NUM_KMER_METH);
        MALLOC_CHK(core->cpgmodel);
    } else {
        core->cpgmodel = NULL;
    }

    //load the model from files
    uint32_t kmer_size=0;
    uint32_t kmer_size_meth=0;
    if (opt.model_file) {
        kmer_size=read_model(core->model, opt.model_file, MODEL_TYPE_NUCLEOTIDE);
    } else {
        if(opt.flag & F5C_RNA){
            if(opt.flag & F5C_R10){
                INFO("%s","builtin RNA004 nucleotide model loaded");
                kmer_size=set_model(core->model, MODEL_ID_RNA_RNA004_NUCLEOTIDE);
            } else {
                INFO("%s","builtin RNA R9 nucleotide model loaded");
                kmer_size=set_model(core->model, MODEL_ID_RNA_R9_NUCLEOTIDE);
            }
        }
        else{
            if(opt.flag & F5C_R10){
                INFO("%s","builtin DNA R10 nucleotide model loaded");
                kmer_size=set_model(core->model, MODEL_ID_DNA_R10_NUCLEOTIDE);
            } else{
                INFO("%s","builtin DNA R9 nucleotide model loaded");
                kmer_size=set_model(core->model, MODEL_ID_DNA_R9_NUCLEOTIDE);
            }
        }
    }
    if(mode == 0){
        if (opt.meth_model_file) {
            kmer_size_meth=read_model(core->cpgmodel, opt.meth_model_file, MODEL_TYPE_METH);
        } else {
            if(opt.flag & F5C_R10){
                INFO("%s","builtin DNA R10 cpg model loaded");
                kmer_size_meth=set_model(core->cpgmodel, MODEL_ID_DNA_R10_CPG);
            } else {
                INFO("%s","builtin DNA R9 cpg model loaded");
                kmer_size_meth=set_model(core->cpgmodel, MODEL_ID_DNA_R9_CPG);
            }
        }
        if(kmer_size != kmer_size_meth){
            ERROR("The k-mer size of the nucleotide model (%d) and the methylation model (%d) should be the same.",kmer_size,kmer_size_meth);
            exit(EXIT_FAILURE);
        }
    }
    core->kmer_size = kmer_size;

    core->opt = opt;

    //realtime0
    core->realtime0=realtime0;

    core->load_db_time=0;
    core->process_db_time=0;

    core->db_bam_time=0;
    core->db_fasta_time=0;
    core->db_fast5_time=0;
    core->db_fast5_open_time=0;
    core->db_fast5_read_time=0;

    core->event_time=0;
    core->align_time=0;
    core->est_scale_time=0;
    core->meth_time=0;

    core->output_time=0;

    //cuda stuff
#ifdef HAVE_CUDA
    if (!(core->opt.flag & F5C_DISABLE_CUDA)) {
        init_cuda(core);
    }
#endif

    core->sum_bases=0;
    core->total_reads=0; //total number mapped entries in the bam file (after filtering based on flags, mapq etc)
    core->bad_fast5_file=0; //empty fast5 path returned by readdb, could not open fast5
    core->ultra_long_skipped=0;
    core->qc_fail_reads=0;
    core->failed_calibration_reads=0;
    core->failed_alignment_reads=0;
    core->skip_mapq_reads=0;
    core->unmapped_reads=0;
    core->skip_sec_reads=0;

    //eventalign related
    core->mode = mode;
    core->read_index=0;
    if(mode==1){
        if(eventalignsummary!=NULL){
            core->event_summary_fp = fopen(eventalignsummary,"w");
            F_CHK(core->event_summary_fp,eventalignsummary);
        }
        else{
            core->event_summary_fp =NULL;
        }

        if(core->opt.flag & F5C_SAM){
            core->sam_output = hts_open("-", "w");
        }
    }

    return core;
}

/* free the core data structure */
void free_core(core_t* core,opt_t opt) {
    free(core->model);
    free(core->cpgmodel);
    delete core->readbb;
    fai_destroy(core->fai);

    if(core->itr){
        sam_itr_destroy(core->itr);
    }
    if(core->reg_list){
        for(int64_t i=0;i<core->reg_n;i++){
            free(core->reg_list[i]);
        }
        free(core->reg_list);
    }

    bam_hdr_destroy(core->m_hdr);
    hts_idx_destroy(core->m_bam_idx);
    sam_close(core->m_bam_fh);

    if(core->ultra_long_tmp!=NULL){
        sam_close(core->ultra_long_tmp);
    }
    if(core->opt.flag&F5C_WR_RAW_DUMP || core->opt.flag&F5C_RD_RAW_DUMP){
        fclose(core->raw_dump);
    }
    if(core->opt.flag & F5C_RD_SLOW5){
        slow5_idx_unload(core->sf);
        slow5_close(core->sf);
    }
#ifdef HAVE_CUDA
    if (!(core->opt.flag & F5C_DISABLE_CUDA)) {
        free_cuda(core);
    }
#endif
    //eventalign related
    if(core->mode==1 && core->event_summary_fp!=NULL){
        fclose(core->event_summary_fp);
    }
    if(core->mode==1 && core->opt.flag & F5C_SAM){
        hts_close(core->sam_output);
    }
    if(opt.num_iop > 1){
        free_iop(core,opt);
    }
    free(core);
}

/* initialise a data batch */
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
    db->read_idx = (int64_t*)(malloc(sizeof(int64_t) * db->capacity_bam_rec));
    MALLOC_CHK(db->read_idx);

    db->sig = (signal_t**)malloc(sizeof(signal_t*) * db->capacity_bam_rec);
    MALLOC_CHK(db->sig);

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
        (int32_t*)malloc(sizeof(int32_t) * db->capacity_bam_rec);
    MALLOC_CHK(db->n_event_alignment);

    db->events_per_base =
        (double*)malloc(sizeof(double) * db->capacity_bam_rec);
    MALLOC_CHK(db->events_per_base);

    db->base_to_event_map =
        (index_pair_t**)malloc(sizeof(index_pair_t*) * db->capacity_bam_rec);
    MALLOC_CHK(db->base_to_event_map);

    db->read_stat_flag = (int32_t *)malloc(sizeof(int32_t) * db->capacity_bam_rec);
    MALLOC_CHK(db->read_stat_flag);

    db->site_score_map = (std::map<int, ScoredSite> **)malloc(sizeof(std::map<int, ScoredSite> *) * db->capacity_bam_rec);
    MALLOC_CHK(db->site_score_map);

    for (i = 0; i < db->capacity_bam_rec; ++i) {
        db->site_score_map[i] = new std::map<int, ScoredSite>;
        NULL_CHK(db->site_score_map[i]);
    }

    db->total_reads=0;
    db->bad_fast5_file=0;
    db->ultra_long_skipped=0;
    db->skip_mapq_reads = 0;
    db->unmapped_reads = 0;
    db->skip_sec_reads = 0;

    //eventalign related
    if(core->mode==1){
        db->eventalign_summary = (EventalignSummary *)malloc(sizeof(EventalignSummary) * db->capacity_bam_rec);
        MALLOC_CHK(db->eventalign_summary);

        db->event_alignment_result = (std::vector<event_alignment_t> **)malloc(sizeof(std::vector<event_alignment_t> *) * db->capacity_bam_rec);
        MALLOC_CHK(db->event_alignment_result);

        db->event_alignment_result_str = (char **)malloc(sizeof(char *) * db->capacity_bam_rec);
        MALLOC_CHK(db->event_alignment_result_str);

        for (i = 0; i < db->capacity_bam_rec; ++i) {
            db->event_alignment_result[i] = new std::vector<event_alignment_t> ;
            NULL_CHK(db->event_alignment_result[i]);
            (db->eventalign_summary[i]).num_events=0; //done here in the same loop for efficiency
            db->event_alignment_result_str[i] = NULL;
        }

    }
    else{
        db->eventalign_summary = NULL;
        db->event_alignment_result = NULL;
        db->event_alignment_result_str = NULL;
    }

    return db;
}

/* load a data batch from disk */
ret_status_t load_db(core_t* core, db_t* db) {
    if(core->opt.num_iop == 1 || (core->opt.flag & F5C_RD_SLOW5) ){
        return load_db1(core,db);
    }
    else{
        if (core->opt.flag & F5C_PRINT_RAW) {
            ERROR("%s","Printing data unsupported with --iop");
            exit(EXIT_FAILURE);
        }
        if (core->opt.flag & F5C_RD_RAW_DUMP){
            ERROR("%s","Reading from raw dump is unsupported with --iop");
            assert(0);
        }
        if(core->opt.flag & F5C_WR_RAW_DUMP){
            ERROR("%s","Writing to raw dump is unsupported with --iop");
            exit(EXIT_FAILURE);
        }
        return load_db2(core,db);
    }
}


#ifdef WORK_STEAL
static inline int32_t steal_work(pthread_arg_t* all_args, int32_t n_threads)
{

	int32_t i, c_i = -1;
	int32_t k;
	for (i = 0; i < n_threads; ++i){
        pthread_arg_t args = all_args[i];
        //fprintf(stderr,"endi : %d, starti : %d\n",args.endi,args.starti);
		if (args.endi-args.starti > STEAL_THRESH) {
            //fprintf(stderr,"gap : %d\n",args.endi-args.starti);
            c_i = i;
            break;
        }
    }
    if(c_i<0){
        return -1;
    }
	k = __sync_fetch_and_add(&(all_args[c_i].starti), 1);
    //fprintf(stderr,"k : %d, end %d, start %d\n",k,all_args[c_i].endi,all_args[c_i].starti);
	return k >= all_args[c_i].endi ? -1 : k;
}
#endif

void* pthread_single(void* voidargs) {
    int32_t i;
    pthread_arg_t* args = (pthread_arg_t*)voidargs;
    db_t* db = args->db;
    core_t* core = args->core;

#ifndef WORK_STEAL
    for (i = args->starti; i < args->endi; i++) {
        args->func(core,db,i);
    }
#else
    pthread_arg_t* all_args = (pthread_arg_t*)(args->all_pthread_args);
    //adapted from kthread.c in minimap2
	for (;;) {
		i = __sync_fetch_and_add(&args->starti, 1);
		if (i >= args->endi) {
            break;
        }
		args->func(core,db,i);
	}
	while ((i = steal_work(all_args,core->opt.num_thread)) >= 0){
		args->func(core,db,i);
    }
#endif

    //fprintf(stderr,"Thread %d done\n",(myargs->position)/THREADS);
    pthread_exit(0);
}


void pthread_db(core_t* core, db_t* db, void (*func)(core_t*,db_t*,int)){

    if (core->opt.num_thread == 1) {
        int i;
        for (i = 0; i < db->n_bam_rec; i++) {
            func(core,db,i);
        }
    }
    else{
        //create threads
        pthread_t tids[core->opt.num_thread];
        pthread_arg_t pt_args[core->opt.num_thread];
        int32_t t, ret;
        int32_t i = 0;
        int32_t num_thread = core->opt.num_thread;
        int32_t step = (db->n_bam_rec + num_thread - 1) / num_thread;
        //todo : check for higher num of threads than the data
        //current works but many threads are created despite

        //set the data structures
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
            pt_args[t].func=func;
        #ifdef WORK_STEAL
            pt_args[t].all_pthread_args =  (void *)pt_args;
        #endif
            //fprintf(stderr,"t%d : %d-%d\n",t,pt_args[t].starti,pt_args[t].endi);

        }

        //create threads
        for(t = 0; t < core->opt.num_thread; t++){
            ret = pthread_create(&tids[t], NULL, pthread_single,
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


void event_single(core_t* core,db_t* db, int32_t i) {

    if(db->sig[i]->nsample>0){

        float* rawptr = db->sig[i]->rawptr;
        float range = db->sig[i]->range;
        float digitisation = db->sig[i]->digitisation;
        float offset = db->sig[i]->offset;
        int32_t nsample = db->sig[i]->nsample;

        // convert to pA
        float raw_unit = range / digitisation;
        for (int32_t j = 0; j < nsample; j++) {
            rawptr[j] = (rawptr[j] + offset) * raw_unit;
        }

        int8_t rna=0;
        if (core->opt.flag & F5C_RNA){
            rna=1;
        }
        db->et[i] = getevents(db->sig[i]->nsample, rawptr, rna);

        // if(db->et[i].n/(float)db->read_len[i] > 20){
        //     fprintf(stderr,"%s\tevents_per_base\t%f\tread_len\t%d\n",bam_get_qname(db->bam_rec[i]), db->et[i].n/(float)db->read_len[i],db->read_len[i]);
        // }

        //get the scalings
        db->scalings[i] = estimate_scalings_using_mom(
            db->read[i], db->read_len[i], core->model, core->kmer_size, db->et[i]);

        //If sequencing RNA, reverse the events to be 3'->5'
        if (rna){
            event_t *events = db->et[i].event;
            size_t n_events = db->et[i].n;
            for (size_t i = 0; i < n_events/2; ++i) {
                event_t tmp_event = events[i];
                events[i]=events[n_events-1-i];
                events[n_events-1-i]=tmp_event;
            }
        }

        //allocate memory for the next alignment step
        db->event_align_pairs[i] = (AlignedPair*)malloc(
                    sizeof(AlignedPair) * (db->et[i].n + db->read_len[i]));
        MALLOC_CHK(db->event_align_pairs[i]);
    }
    else{
        db->et[i].n = 0;
        db->et[i].event = NULL;
        db->event_align_pairs[i] = NULL;
    }

}

void scaling_single(core_t* core, db_t* db, int32_t i){

    db->event_alignment[i] = NULL;
    db->n_event_alignment[i] = 0;
    db->events_per_base[i] = 0; //todo : is double needed? not just int8?

    int32_t n_kmers = db->read_len[i] - core->kmer_size + 1;

    if (db->n_event_align_pairs[i] > 0) {

        db->base_to_event_map[i]=(index_pair_t*)(malloc(sizeof(index_pair_t) * n_kmers));
        MALLOC_CHK(db->base_to_event_map[i]);

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
            n_kmers, db->event_align_pairs[i], db->n_event_align_pairs[i], core->kmer_size);

        //fprintf(stderr,"n_event_alignment %d\n",n_events);

        // run recalibration to get the best set of scaling parameters and the residual
        // between the (scaled) event levels and the model.

        // internally this function will set shift/scale/etc of the pore model
        bool calibrated = recalibrate_model(
            core->model, core->kmer_size, db->et[i], &db->scalings[i],
            db->event_alignment[i], db->n_event_alignment[i], 1, core->opt.min_num_events_to_rescale);

        // QC calibration
        if (!calibrated || db->scalings[i].var > MIN_CALIBRATION_VAR) {
            //     events[strand_idx].clear();
            free(db->event_alignment[i]);
            //free(db->event_align_pairs[i]);
            db->read_stat_flag[i] |= FAILED_CALIBRATION;
            return;
        }

        free(db->event_alignment[i]);

    } else {
        db->base_to_event_map[i] = NULL;
        // Could not align, fail this read
        // this->events[strand_idx].clear();
        // this->events_per_base[strand_idx] = 0.0f;
        //free(db->event_align_pairs[i]);
        db->read_stat_flag[i] |= FAILED_ALIGNMENT;
        return;
    }

    // Filter poor quality reads that have too many "stays"

    if (db->events_per_base[i] > 5.0) {
        //     events[0].clear();
        //     events[1].clear();
        //free(db->event_align_pairs[i]);
        db->read_stat_flag[i] |= FAILED_QUALITY_CHK;
        return;
    }


}

/* align a single read specified by index i (perform ABEA for a single read) */
//note that this is used in f5c.cu and thus modifications must be done with care
void align_single(core_t* core, db_t* db, int32_t i) {

    if(db->sig[i]->nsample>0){ //if a good read
        if (db->sig[i]->nsample && (db->et[i].n)/(float)(db->read_len[i]) < AVG_EVENTS_PER_KMER_MAX){
            db->n_event_align_pairs[i] = align(
                    db->event_align_pairs[i], db->read[i], db->read_len[i], db->et[i],
                    core->model, core->kmer_size, db->scalings[i], db->sig[i]->sample_rate);
                //fprintf(stderr,"readlen %d,n_events %d\n",db->read_len[i],n_event_align_pairs);
        }
        else{//todo : too many avg events per base - oversegmented
            db->n_event_align_pairs[i]=0;
            if(core->opt.verbosity > 0){
                STDERR("Skipping over-segmented read %s with %f events per base",bam_get_qname(db->bam_rec[i]), (db->et[i].n)/(float)(db->read_len[i]));
            }
        }
    }
    else{ //if a bad read (corrupted/missing slow5 record)
       db->n_event_align_pairs[i]=0;
    }
}

/* align a data batch (perform ABEA for a data batch) */
void align_db(core_t* core, db_t* db) {
#ifdef HAVE_CUDA
    if (!(core->opt.flag & F5C_DISABLE_CUDA)) {
        //STDERR("%s","Performing on cuda");
        align_cuda(core, db);
    }
#endif

    if (core->opt.flag & F5C_DISABLE_CUDA) {
        //fprintf(stderr, "cpu\n");
        pthread_db(core, db, align_single);
    }
}


void eventalign_single(core_t* core, db_t* db, int32_t i){
    realign_read(db->event_alignment_result[i], &(db->eventalign_summary[i]),core->event_summary_fp, db->fasta_cache[i],core->m_hdr,
                  db->bam_rec[i],db->read_len[i],
                  i,
                  core->clip_start,
                  core->clip_end,
                  &(db->et[i]), core->model,core->kmer_size, db->base_to_event_map[i],db->scalings[i],db->events_per_base[i],db->sig[i]->sample_rate);

    char* qname = bam_get_qname(db->bam_rec[i]);
    char* contig = core->m_hdr->target_name[db->bam_rec[i]->core.tid];
    std::vector<event_alignment_t> *event_alignment_result = db->event_alignment_result[i];
    int8_t print_read_names = (core->opt.flag & F5C_PRINT_RNAME) ? 1 : 0;
    int8_t scale_events = (core->opt.flag & F5C_SCALE_EVENTS) ? 1 : 0;
    int8_t collapse_events = (core->opt.flag & F5C_COLLAPSE_EVENTS) ? 1 : 0;
    int8_t write_samples = (core->opt.flag & F5C_PRINT_SAMPLES) ? 1 : 0;
    int8_t write_signal_index = (core->opt.flag & F5C_PRINT_SIGNAL_INDEX) ? 1 : 0;
    int8_t sam_output = (core->opt.flag & F5C_SAM) ? 1 : 0;
    int8_t paf_output = (core->opt.flag & F5C_PAF) ? 1 : 0;
    int8_t m6anet_output = (core->opt.flag & F5C_M6ANET) ? 1 : 0;
    int8_t rna = (core->opt.flag & F5C_RNA) ? 1 : 0;

    if(paf_output){
        int64_t ref_len = core->m_hdr->target_len[db->bam_rec[i]->core.tid];
        db->event_alignment_result_str[i] = emit_event_alignment_paf(&(db->et[i]),db->sig[i]->nsample, ref_len,core->kmer_size, db->scalings[i],*event_alignment_result, db->bam_rec[i], qname, contig, rna);
    } else if(sam_output){
        int8_t sam_out_version = core->opt.sam_out_version;
        int64_t ref_len = core->m_hdr->target_len[db->bam_rec[i]->core.tid];
        db->event_alignment_result_str[i] = emit_event_alignment_sam(qname, core->m_hdr, db->bam_rec[i], *event_alignment_result, sam_out_version, &(db->et[i]), db->sig[i]->nsample, ref_len, rna, db->scalings[i]);
    } else if (m6anet_output){
        db->event_alignment_result_str[i] = emit_event_alignment_tsv_m6anet(0,&(db->et[i]),core->model,core->kmer_size, db->scalings[i],*event_alignment_result, print_read_names, scale_events, write_samples, write_signal_index, collapse_events,
                   db->read_idx[i], qname, contig, db->sig[i]->sample_rate, db->sig[i]->rawptr);
    } else {
        db->event_alignment_result_str[i] = emit_event_alignment_tsv(0,&(db->et[i]),core->model,core->kmer_size, db->scalings[i],*event_alignment_result, print_read_names, scale_events, write_samples, write_signal_index, collapse_events,
                   db->read_idx[i], qname, contig, db->sig[i]->sample_rate, db->sig[i]->rawptr);
    }
}

void meth_single(core_t* core, db_t* db, int32_t i){
    if(!db->read_stat_flag[i]){
        if(core->mode==0){
            calculate_methylation_for_read(db->site_score_map[i], db->fasta_cache[i], db->bam_rec[i], db->read_len[i], db->et[i].event, db->base_to_event_map[i],
            db->scalings[i], core->cpgmodel, core->kmer_size, db->events_per_base[i]);
        }
        else if (core->mode==1){
            eventalign_single(core,db,i);
        }
    }
}


void process_single(core_t* core, db_t* db,int32_t i) {
    event_single(core,db,i);
    align_single(core, db,i);
    scaling_single(core,db,i);
    meth_single(core,db,i);
}

/* completely process a data batch
   (all steps: event detection, adaptive banded event alignment, ...., HMM) */
void process_db(core_t* core, db_t* db) {

    double process_start = realtime();

    if((core->opt.flag&F5C_SEC_PROF) || (!(core->opt.flag & F5C_DISABLE_CUDA))){

        double realtime0=core->realtime0;

        double event_start = realtime();
        pthread_db(core,db,event_single);
        double event_end = realtime();
        core->event_time += (event_end-event_start);

        fprintf(stderr, "[%s::%.3f*%.2f] Events computed\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0));

        double align_start = realtime();
        align_db(core, db);
        double align_end = realtime();
        core->align_time += (align_end-align_start);

        fprintf(stderr, "[%s::%.3f*%.2f] Banded alignment done\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0));

        double est_scale_start = realtime();
        pthread_db(core,db,scaling_single);
        double est_scale_end = realtime();
        core->est_scale_time += (est_scale_end-est_scale_start);

        fprintf(stderr, "[%s::%.3f*%.2f] Scaling calibration done\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0));

        double meth_start = realtime();
        pthread_db(core, db, meth_single);
        double meth_end = realtime();
        core->meth_time += (meth_end-meth_start);

        fprintf(stderr, "[%s::%.3f*%.2f] HMM done\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0));


    }
    else{
        if (core->opt.num_thread == 1) {
            int32_t i=0;
            for (i = 0; i < db->n_bam_rec; i++) {
                process_single(core,db,i);
            }

        }
        else {
            pthread_db(core,db,process_single);
        }

    }

    double process_end= realtime();
    core->process_db_time += (process_end-process_start);

    return;
}

/* write the output for a processed data batch */
void output_db(core_t* core, db_t* db) {

    double output_start = realtime();

    if (core->opt.flag & F5C_PRINT_EVENTS) {
        int32_t i = 0;
        for (i = 0; i < db->n_bam_rec; i++) {
            printf(">%s\tLN:%d\tEVENTSTART:%d\tEVENTEND:%d\n",
                   bam_get_qname(db->bam_rec[i]), (int)db->et[i].n,
                   (int)db->et[i].start, (int)db->et[i].end);
            uint32_t j = 0;
            for (j = 0; j < db->et[i].n; j++) {
                printf("{%d,%f,%f,%f}\t", (int)db->et[i].event[j].start,
                       db->et[i].event[j].length, db->et[i].event[j].mean,
                       db->et[i].event[j].stdv);
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
            printf(">%s\tN_ALGN_PAIR:%d\t{ref_pos,read_pos}\n",
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

    core->sum_bases += db->sum_bases;
    core->total_reads += db->total_reads;
    core->bad_fast5_file += db->bad_fast5_file;
    core->ultra_long_skipped += db->ultra_long_skipped;
    core->skip_mapq_reads += db->skip_mapq_reads;
    core->skip_sec_reads += db->skip_sec_reads;
    core->unmapped_reads += db->unmapped_reads;

    int32_t i = 0;
    for (i = 0; i < db->n_bam_rec; i++){
        if(!db->read_stat_flag[i]){
            char* qname = bam_get_qname(db->bam_rec[i]);
            char* contig = core->m_hdr->target_name[db->bam_rec[i]->core.tid];

            if(core->mode==0) {
                std::map<int, ScoredSite> *site_score_map = db->site_score_map[i];
                // write all sites for this read
                for(auto iter = site_score_map->begin(); iter != site_score_map->end(); ++iter) {

                    const ScoredSite& ss = iter->second;
                    double sum_ll_m = ss.ll_methylated[0]; //+ ss.ll_methylated[1];
                    double sum_ll_u = ss.ll_unmethylated[0]; //+ ss.ll_unmethylated[1];
                    double diff = sum_ll_m - sum_ll_u;

                    // output only if inside the window boundaries
                    if( !( (core->clip_start != -1 && ss.start_position < core->clip_start) ||
                        (core->clip_end != -1 && ss.end_position >= core->clip_end) ) ) {
                        if(core->opt.meth_out_version==1){
                            printf("%s\t%d\t%d\t", contig, ss.start_position, ss.end_position);
                        }
                        else if(core->opt.meth_out_version==2){
                            printf("%s\t%c\t%d\t%d\t", contig, bam_is_rev(db->bam_rec[i]) ? '-' : '+', ss.start_position, ss.end_position);
                        }
                        printf("%s\t%.2lf\t", qname, diff);
                        printf("%.2lf\t%.2lf\t", sum_ll_m, sum_ll_u);
                        printf("%d\t%d\t%s\n", ss.strands_scored, ss.n_cpg, ss.sequence.c_str());
                    }

                }
            }

            else if(core->mode==1){
                FILE* summary_fp = core->event_summary_fp;
                EventalignSummary summary = db->eventalign_summary[i];
                scalings_t scalings = db->scalings[i];
                if(summary_fp != NULL && summary.num_events > 0) {
                    size_t strand_idx = 0;
                    std::string fast5_path_str = core->readbb->get_signal_path(qname);
                    const char *path = (core->opt.flag & F5C_RD_SLOW5)  ? "slow5" : fast5_path_str.c_str();
                    fprintf(summary_fp, "%ld\t%s\t", (long)(db->read_idx[i]), qname);
                    fprintf(summary_fp, "%s\t%s\t%s\t",path, (core->opt.flag & F5C_RNA) ? "rna" : "dna", strand_idx == 0 ? "template" : "complement" );
                    fprintf(summary_fp, "%d\t%d\t%d\t%d\t", summary.num_events, summary.num_steps, summary.num_skips, summary.num_stays);
                    fprintf(summary_fp, "%.2lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\n", summary.sum_duration/(db->sig[i]->sample_rate), scalings.shift, scalings.scale, 0.0, scalings.var);
                }

                char *event_alignment_result_str = db->event_alignment_result_str[i];
                fputs(event_alignment_result_str,stdout);

            }
        }
        else{
            if((db->read_stat_flag[i])&FAILED_CALIBRATION){
                core->failed_calibration_reads++;
            }
            else if ((db->read_stat_flag[i])&FAILED_ALIGNMENT){
                core->failed_alignment_reads++;
            }
            else if ((db->read_stat_flag[i])&FAILED_QUALITY_CHK){
                core->qc_fail_reads++;
            }
            else{
                assert(0);
            }
        }
    }
    fflush(stdout);

    //core->read_index = core->read_index + db->n_bam_rec;
    double output_end = realtime();
    core->output_time += (output_end-output_start);

}

/* partially free a data batch - only the read dependent allocations are freed */
void free_db_tmp(db_t* db) {
    int32_t i = 0;
    for (i = 0; i < db->n_bam_rec; ++i) {
        bam_destroy1(db->bam_rec[i]);
        db->bam_rec[i] = bam_init1();
        free(db->fasta_cache[i]);
        free(db->read[i]);
        free(db->sig[i]->rawptr);
        free(db->sig[i]);
        free(db->et[i].event);
        free(db->event_align_pairs[i]);
        free(db->base_to_event_map[i]);
        delete db->site_score_map[i];
        db->site_score_map[i] = new std::map<int, ScoredSite>;

        if(db->event_alignment_result){ //eventalign related
            delete db->event_alignment_result[i];
            db->event_alignment_result[i] = new std::vector<event_alignment_t>;
        }
        if(db->event_alignment_result_str){ //eventalign related
            free(db->event_alignment_result_str[i]);
            db->event_alignment_result_str[i]=NULL;
        }
    }
}

/* completely free a data batch */
void free_db(db_t* db) {
    int32_t i = 0;
    for (i = 0; i < db->capacity_bam_rec; ++i) {
        bam_destroy1(db->bam_rec[i]);
    }
    free(db->bam_rec);
    free(db->fasta_cache);
    free(db->read);
    free(db->read_len);
    free(db->read_idx);
    free(db->et);
    free(db->sig);
    free(db->scalings);
    free(db->event_align_pairs);
    free(db->n_event_align_pairs);
    free(db->event_alignment);
    free(db->n_event_alignment);
    free(db->events_per_base);
    free(db->base_to_event_map);
    free(db->read_stat_flag);
    for (i = 0; i < db->capacity_bam_rec; ++i) {
        delete db->site_score_map[i];
    }
    free(db->site_score_map);
    //eventalign related
    if(db->eventalign_summary){
        free(db->eventalign_summary);
    }
    if(db->event_alignment_result){
        for (i = 0; i < db->capacity_bam_rec; ++i) {
            delete db->event_alignment_result[i];
        }
        free(db->event_alignment_result);
    }
    if(db->event_alignment_result_str){
        free(db->event_alignment_result_str);
    }

    free(db);
}

/* initialise user specified options */
void init_opt(opt_t* opt) {
    memset(opt, 0, sizeof(opt_t));
    opt->min_mapq = 20;
    opt->batch_size = 512;
    opt->batch_size_bases = 2*1000*1000;
    opt->pore = NULL;
    opt->num_thread = 8;
    opt->num_iop = 1;       //if changed, the SLOW5 mode must be handled by default in the arg parsing
    opt->region_str = NULL; //whole genome processing if null

    opt->min_num_events_to_rescale = 200;

#ifndef HAVE_CUDA
    opt->flag |= F5C_DISABLE_CUDA;
    opt->batch_size_bases = 5*1000*1000;
#endif

    opt->flag |= F5C_SKIP_UNREADABLE;
    opt->debug_break=-1;
    opt->ultra_thresh=100000;

    opt->meth_out_version=2;
    opt->sam_out_version=2;

    opt->cuda_block_size=64;
    opt->cuda_dev_id=0;
    opt->cuda_mem_frac=1.0f; //later set by cuda_init()

    //effective only if  CPU_GPU_PROC  is set
    opt->cuda_max_readlen=3.0f;
    opt->cuda_avg_events_per_kmer=2.0f; //only if CUDA_DYNAMIC_MALLOC is unset
    opt->cuda_max_avg_events_per_kmer=5.0f;
}
