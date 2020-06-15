/* @f5c
**
** f5c interface
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
#include "profiles.h"

#include <sys/wait.h>
#include <unistd.h>

/*
todo :
Error counter for consecutive failures in the skip unreadable mode
not all the memory allocations are needed for eventalign mode
*/


//make this inline for performance reasons
void f5write(FILE* fp, void *buf, size_t element_size, size_t num_elements){
	size_t ret=fwrite(buf,element_size,num_elements,fp);
	if(ret!=num_elements){
		fprintf(stderr,"Writing error has occurred :%s\n",strerror(errno));
		exit(EXIT_FAILURE);
	}
}

static inline void f5read(FILE* fp, void *buf, size_t element_size, size_t num_elements){
	size_t ret=fread(buf,element_size,num_elements,fp);
	if(ret!=num_elements){
		fprintf(stderr,"Reading error has occurred :%s\n",strerror(errno));
		exit(EXIT_FAILURE);
	}
}

static inline int read_from_fast5_files2(char *qname_str, char *fast5_path, FILE *pipefp){

    //fprintf(stderr,"readname : %s\n",qname.c_str());,
    std::string qname = qname_str;
    int8_t success=0;


    //INFO("Opening %s",fast5_path);
    fast5_file_t fast5_file = fast5_open(fast5_path);

    if (fast5_file.hdf5_file >= 0) {
        fast5_t *f5 = (fast5_t*)calloc(1, sizeof(fast5_t));
        MALLOC_CHK(f5);

        //INFO("Reading %s",qname.c_str());
        int32_t ret=fast5_read(fast5_file, f5,qname);

        if(ret<0){
            ERROR("%s","Fast5 causes crashes with --iop");
            exit(EXIT_FAILURE);
        }

        fast5_close(fast5_file);
        success=1;

        //write to the pipe
        //STDERR("writing to pipe %s",qname_str);
        f5write(pipefp,&(f5->nsample), sizeof(hsize_t), 1);
        f5write(pipefp,f5->rawptr, sizeof(float), f5->nsample);
        f5write(pipefp,&(f5->digitisation), sizeof(float), 1);
        f5write(pipefp,&(f5->offset), sizeof(float), 1);
        f5write(pipefp,&(f5->range), sizeof(float), 1);
        f5write(pipefp,&(f5->sample_rate), sizeof(float), 1);

        ret=fflush(pipefp);
        if(ret!=0){
            ERROR("%s","Flushing the pipe failed");
            exit(EXIT_FAILURE);
        }

        free(f5->rawptr);
        free(f5);
        //STDERR("wrote to pipe %s : %lld samples",qname_str,f5->nsample);
    } else {
        ERROR("%s","Fast5 causes crashes with --iop");
        exit(EXIT_FAILURE);
    }
    assert(success==1);


    return 1;
}


//this is the child handler
void iop_handler(FILE *pipefp_c2p, FILE *pipefp_p2c){

    //buffers for geline
    char* buffer =
        (char*)malloc(sizeof(char) * (100)); //READ+newline+nullcharacter
    MALLOC_CHK(buffer);
    size_t bufferSize = 100;
    ssize_t readlinebytes = 0;

    //STDERR("IOP handler %d",0);

    int i=0;
    while ((readlinebytes = getline(&buffer, &bufferSize, pipefp_p2c)) != -1) {

        //STDERR("%d th read reading",i);
        //fprintf(stderr,"buffer: %s\n",buffer);
        char *qname = strtok(buffer,"\t");
        char *fast5_path_str = strtok(NULL,"\n");
        //fprintf(stderr,"%s\t%s\n",qname,fast5_path_str);
        if(qname==NULL || fast5_path_str==NULL){
            ERROR("%s","Malformed pipe format");
            exit(EXIT_FAILURE);
        }

        int8_t read_status = 0;
        read_status=read_from_fast5_files2(qname,fast5_path_str,pipefp_c2p);
        if(read_status!=1){
            assert(0);
        }
        //STDERR("%d th read %s read",i,qname);

        i++;
    }

    //STDERR("IOP handler exiting%d",0);

    free(buffer);
}

//raise the children
void init_iop(core_t* core,opt_t opt){

    STDERR("Spawning %d I/O processes to circumvent HDF hell",opt.num_iop);

    core->pids = (pid_t *) malloc(sizeof(pid_t)*opt.num_iop);
    MALLOC_CHK(core->pids);

    core->pipefd_p2c = (int *) malloc(sizeof(int)*opt.num_iop);
    MALLOC_CHK(core->pipefd_p2c);
    core->pipefd_c2p = (int *) malloc(sizeof(int)*opt.num_iop);
    MALLOC_CHK(core->pipefd_c2p);

    core->pipefp_p2c = (FILE **) malloc(sizeof(FILE *)*opt.num_iop);
    MALLOC_CHK(core->pipefp_p2c);
    core->pipefp_c2p = (FILE **) malloc(sizeof(FILE *)*opt.num_iop);
    MALLOC_CHK(core->pipefp_c2p);

    int pipefd_p2c[opt.num_iop][2]; //parent to child
    int pipefd_c2p[opt.num_iop][2]; //child to parent

    //create processes
    int i;
    for(i=0;i<opt.num_iop;i++){
        if (pipe(pipefd_p2c[i]) == -1) { //parent to child pipe
            perror("pipe");
            exit(EXIT_FAILURE);
        }
        if (pipe(pipefd_c2p[i]) == -1) { //child to parent pipe
            perror("pipe");
            exit(EXIT_FAILURE);
        }

        core->pids[i] = fork();
        if(core->pids[i]==-1){
            ERROR("%s","Fork failed");
            perror("");
            exit(EXIT_FAILURE);
        }
        if(core->pids[i]==0){ //child
            if(opt.verbosity>1){
                STDERR("%s","child : child spawned");
            }
            close(pipefd_c2p[i][0]); //close read end of child to parent
            close(pipefd_p2c[i][1]); //close write end of parent to child

            int j = 0;
            for(j=0;j<i;j++){
                fclose(core->pipefp_p2c[j]);
                fclose(core->pipefp_c2p[j]);
                close(core->pipefd_p2c[j]);
                close(core->pipefd_c2p[j]);
            }

            FILE *pipefp_c2p = fdopen(pipefd_c2p[i][1],"w");
            F_CHK(pipefp_c2p,"child : pipefp child to parent");
            FILE *pipefp_p2c = fdopen(pipefd_p2c[i][0],"r");
            F_CHK(pipefp_p2c,"child : pipefp parent to child");

            //This does the actual work of accepting requests to fast5 and serving them
            iop_handler(pipefp_c2p, pipefp_p2c);

            fclose(pipefp_c2p);
            fclose(pipefp_p2c);

            close(pipefd_c2p[i][1]);
            close(pipefd_p2c[i][0]);

            free(core->pids);
            free(core->pipefd_p2c);
            free(core->pipefd_c2p);
            free(core->pipefp_p2c);
            free(core->pipefp_c2p);
            free(core);

            if(opt.verbosity>1){
                INFO("%s","child : child is exiting");
            }

            //TODO : free the datastructures allocated in the children such as core
            //TODO : error handling of parent in case a child crashes and vice versa

            exit(EXIT_SUCCESS);
        }
        if(core->pids[i]>0){ //parent
            close(pipefd_c2p[i][1]); //close write end of child to parent
            close(pipefd_p2c[i][0]); //close read end of parent to child
            if(opt.verbosity>1){
                STDERR("parent : child process with pid %d created", core->pids[i]);
            }

            core->pipefd_p2c[i] = pipefd_p2c[i][1];
            core->pipefd_c2p[i] = pipefd_c2p[i][0];

            core->pipefp_p2c[i] = fdopen(pipefd_p2c[i][1],"w");;
            F_CHK(core->pipefp_p2c[i],"parent : pipefp parent to child");

            core->pipefp_c2p[i] = fdopen(pipefd_c2p[i][0],"r");
            F_CHK(core->pipefp_c2p[i],"parent : pipefp child to parent");
        }

    }
}

//get rid of the children
void free_iop(core_t* core,opt_t opt){

    int i;
    for(i=0;i<opt.num_iop;i++){
        int ret = fclose(core->pipefp_p2c[i]);
        if(ret!=0){
            WARNING("%s","Closeing the pipe p2c failed");
            perror("");
        }
        ret = fclose(core->pipefp_c2p[i]);
        if(ret!=0){
            WARNING("%s","Closeing the pipe c2p failed");
            perror("");
        }
        close(core->pipefd_p2c[i]);
        close(core->pipefd_c2p[i]);
    }

    int status,w;
    for(i=0;i<opt.num_iop;i++){

        if(core->opt.verbosity>1){
            STDERR("parent : Waiting for child with pid %d",core->pids[i]);
        }

        do {
            w = waitpid(core->pids[i], &status, WUNTRACED | WCONTINUED);
            if (w == -1) {
                ERROR("%s","waitpid failed");
                perror("");
                exit(EXIT_FAILURE);
            }

            if(core->opt.verbosity>1){
                if (WIFEXITED(status)) {
                    STDERR("child process %d exited, status=%d", core->pids[i], WEXITSTATUS(status));
                } else if (WIFSIGNALED(status)) {
                    STDERR("child process %d killed by signal %d", core->pids[i], WTERMSIG(status));
                } else if (WIFSTOPPED(status)) {
                    STDERR("child process %d stopped by signal %d", core->pids[i], WSTOPSIG(status));
                } else if (WIFCONTINUED(status)) {
                    STDERR("child process %d continued",core->pids[i]);
                }
            }
        } while (!WIFEXITED(status) && !WIFSIGNALED(status));
    }

    //TODO incase graceful exit fails kill the children

    // //int i;
    // for(i=0;i<opt.num_iop;i++){
    //     kill(core->pids[i], SIGTERM);
    // }

    //TODO : handle the kill if termination fails


    free(core->pipefd_c2p);
    free(core->pipefd_p2c);
    free(core->pipefp_c2p);
    free(core->pipefp_p2c);
    free(core->pids);


}

core_t* init_core(const char* bamfilename, const char* fastafile,
                  const char* fastqfile, const char* tmpfile, opt_t opt,double realtime0, int8_t mode, char *eventalignsummary) {
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
    if(opt.region_str == NULL){
        core->itr = sam_itr_queryi(core->m_bam_idx, HTS_IDX_START, 0, 0);
        if(core->itr==NULL){
            ERROR("%s","sam_itr_queryi failed. A problem with the BAM index?");
            exit(EXIT_FAILURE);
        }
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

    // reference file
    core->fai = fai_load(fastafile);
    NULL_CHK(core->fai);

    // readbb
    core->readbb = new ReadDB;
    core->readbb->load(fastqfile);

    //model
    core->model = (model_t*)malloc(sizeof(model_t) * NUM_KMER); //4096 is 4^6 which is hardcoded now
    MALLOC_CHK(core->model);
    core->cpgmodel = (model_t*)malloc(sizeof(model_t) * NUM_KMER_METH); //15625 is 4^6 which os hardcoded now
    MALLOC_CHK(core->cpgmodel);

    //load the model from files
    if (opt.model_file) {
        read_model(core->model, opt.model_file);
    } else {
        set_model(core->model);
    }

    //todo (low priority) : load the cpg model from file
    set_cpgmodel(core->cpgmodel);

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

void free_core(core_t* core,opt_t opt) {
    free(core->model);
    free(core->cpgmodel);
    delete core->readbb;
    fai_destroy(core->fai);
    sam_itr_destroy(core->itr);
    bam_hdr_destroy(core->m_hdr);
    hts_idx_destroy(core->m_bam_idx);
    sam_close(core->m_bam_fh);
    if(core->ultra_long_tmp!=NULL){
        sam_close(core->ultra_long_tmp);
    }
    if(core->opt.flag&F5C_WR_RAW_DUMP || core->opt.flag&F5C_RD_RAW_DUMP){
        fclose(core->raw_dump);
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

    //eventalign related
    if(core->mode==1){
        db->eventalign_summary = (EventalignSummary *)malloc(sizeof(EventalignSummary) * db->capacity_bam_rec);
        MALLOC_CHK(db->eventalign_summary);

        db->event_alignment_result = (std::vector<event_alignment_t> **)malloc(sizeof(std::vector<event_alignment_t> *) * db->capacity_bam_rec);
        MALLOC_CHK(db->event_alignment_result);
        for (i = 0; i < db->capacity_bam_rec; ++i) {
            db->event_alignment_result[i] = new std::vector<event_alignment_t> ;
            NULL_CHK(db->event_alignment_result[i]);
            (db->eventalign_summary[i]).num_events=0; //done here in the same loop for efficiency
        }
    }
    else{
        db->eventalign_summary = NULL;
        db->event_alignment_result = NULL;
    }

    return db;
}

static inline void handle_bad_fast5(core_t* core, db_t* db,std::string fast5_path_str, std::string qname){
    db->bad_fast5_file++;
    if (core->opt.flag & F5C_SKIP_UNREADABLE) {
        WARNING("Fast5 file [%s] for read [%s] is unreadable and will be skipped",
                fast5_path_str.c_str(),qname.c_str());
    } else {
        ERROR("Fast5 file [%s] could not be opened for read [%s]", fast5_path_str.c_str(), qname.c_str());
        exit(EXIT_FAILURE);
    }
    return;
}


static inline int read_from_fast5_dump(core_t *core, db_t *db , int32_t i){

    //return 1 if success, 0 if failed
    db->f5[i] = (fast5_t*)calloc(1, sizeof(fast5_t));
    MALLOC_CHK(db->f5[i]);

    f5read(core->raw_dump,&(db->f5[i]->nsample), sizeof(hsize_t), 1);

    if(db->f5[i]->nsample>0){
        db->f5[i]->rawptr = (float*)calloc(db->f5[i]->nsample, sizeof(float));
        MALLOC_CHK( db->f5[i]->rawptr);
        f5read(core->raw_dump,db->f5[i]->rawptr, sizeof(float), db->f5[i]->nsample);
        f5read(core->raw_dump,&(db->f5[i]->digitisation), sizeof(float), 1);
        f5read(core->raw_dump,&(db->f5[i]->offset), sizeof(float), 1);
        f5read(core->raw_dump,&(db->f5[i]->range), sizeof(float), 1);
        f5read(core->raw_dump,&(db->f5[i]->sample_rate), sizeof(float), 1);
        return 1;
    }
    else{
        return 0;
    }


}

static inline int read_from_fast5_files(core_t *core, db_t *db, std::string qname, std::string fast5_path_str, int32_t i){
    char* fast5_path =
        (char*)malloc(fast5_path_str.size() + 10); // is +10 needed? do errorcheck
    strcpy(fast5_path, fast5_path_str.c_str());

    //fprintf(stderr,"readname : %s\n",qname.c_str());
    int8_t success=0;

    double t = realtime();
    fast5_file_t fast5_file = fast5_open(fast5_path);
    double ot = realtime() - t;
    core->db_fast5_open_time += ot;
    core->db_fast5_time += ot;
    if (fast5_file.hdf5_file >= 0) {
        db->f5[i] = (fast5_t*)calloc(1, sizeof(fast5_t));
        MALLOC_CHK(db->f5[i]);
        t = realtime();
        int32_t ret=fast5_read(fast5_file, db->f5[i],qname);
        double rt = realtime() - t;
        core->db_fast5_read_time += rt;
        core->db_fast5_time += rt;
        if(ret<0){
            handle_bad_fast5(core, db,fast5_path,qname);
            if(core->opt.flag & F5C_WR_RAW_DUMP){
                hsize_t tmp_nsample = 0;
                f5write(core->raw_dump,&tmp_nsample, sizeof(hsize_t), 1);
            }
            free(fast5_path);
            return 0;
        }
        t = realtime();
        fast5_close(fast5_file);
        core->db_fast5_time += realtime() - t;

        if (core->opt.flag & F5C_PRINT_RAW) {
            printf(">%s\tPATH:%s\tLN:%llu\n", qname.c_str(), fast5_path,
                db->f5[i]->nsample);
            uint32_t j = 0;
            for (j = 0; j < db->f5[i]->nsample; j++) {
                printf("%d\t", (int)db->f5[i]->rawptr[j]);
            }
            printf("\n");
        }
        if(core->opt.flag & F5C_WR_RAW_DUMP){
            //write the fast5 dump to the binary file pointer core->raw_dump
            f5write(core->raw_dump,&(db->f5[i]->nsample), sizeof(hsize_t), 1);
            f5write(core->raw_dump,db->f5[i]->rawptr, sizeof(float), db->f5[i]->nsample);
            f5write(core->raw_dump,&(db->f5[i]->digitisation), sizeof(float), 1);
            f5write(core->raw_dump,&(db->f5[i]->offset), sizeof(float), 1);
            f5write(core->raw_dump,&(db->f5[i]->range), sizeof(float), 1);
            f5write(core->raw_dump,&(db->f5[i]->sample_rate), sizeof(float), 1);
        }

        //db->n_bam_rec++;
        //t = realtime();
        //status.num_bases += read_length;
        //core->db_fasta_time += realtime() - t;
        success=1;
    } else {
        handle_bad_fast5(core, db,fast5_path,qname);
        if(core->opt.flag & F5C_WR_RAW_DUMP){
            hsize_t tmp_nsample = 0;
            f5write(core->raw_dump,&tmp_nsample, sizeof(hsize_t), 1);
        }
        return 0;
    }
    free(fast5_path);
    assert(success==1);
    return 1;
}

ret_status_t load_db1(core_t* core, db_t* db) { //old method

    double load_start = realtime();

    // get bams
    bam1_t* record;
    int32_t result = 0;
    db->n_bam_rec = 0;
    db->sum_bases = 0;
    db->total_reads = 0;
    db->bad_fast5_file = 0;
    db->ultra_long_skipped =0;

    ret_status_t status={0,0};
    int32_t i = 0;
    double t = 0;
    while (db->n_bam_rec < db->capacity_bam_rec && status.num_bases<core->opt.batch_size_bases) {
        i=db->n_bam_rec;
        record = db->bam_rec[i];
        t = realtime();
        result = sam_itr_next(core->m_bam_fh, core->itr, record);
        core->db_bam_time += realtime() - t;

        //set read index
        db->read_idx[i]= core->read_index;
        core->read_index +=1;

        if (result < 0) {
            break;
        } else {
            if ((record->core.flag & BAM_FUNMAP) == 0 &&
                record->core.qual >= core->opt.min_mapq) {
                // printf("%s\t%d\n",bam_get_qname(db->bam_rec[db->n_bam_rec]),result);

                if(!(core->opt.flag & F5C_SECONDARY_YES)){
                    if((record->core.flag & BAM_FSECONDARY)){
                        continue;
                    }
                }

                db->total_reads++; // candidate read

                std::string qname = bam_get_qname(record);
                t = realtime();
                //todo : make efficient (redudantly accessed below, can be combined with it?)
                int64_t read_length=core->readbb->get_read_sequence(qname).size();
                std::string fast5_path_str = core->readbb->get_signal_path(qname);
                core->db_fasta_time += realtime() - t;

                //skipping ultra-long-reads
                if(core->ultra_long_tmp!=NULL && read_length > core->opt.ultra_thresh){
                    db->ultra_long_skipped++;
                    int ret_wr=sam_write1(core->ultra_long_tmp,core->m_hdr,record);
                    NEG_CHK(ret_wr);
                    continue;
                }

                if(fast5_path_str==""){
                    handle_bad_fast5(core, db,fast5_path_str,qname);
                    continue;
                }

                int8_t read_status = 0;
                if (core->opt.flag & F5C_RD_RAW_DUMP){
                    t = realtime();
                    read_status=read_from_fast5_dump(core, db,i);
                    double rt = realtime() - t;
                    core->db_fast5_read_time += rt;
                    core->db_fast5_time += rt;
                }
                else{
                   read_status=read_from_fast5_files(core, db, qname,fast5_path_str,i);
                }
                if(read_status==1){
                    db->n_bam_rec++;
                    status.num_bases += read_length;
                }

            }
        }
    }
    // fprintf(stderr,"%s:: %d queries read\n",__func__,db->n_bam_rec);

    // get ref sequences (todo can make efficient by taking the the start and end of the sorted bam)
    for (i = 0; i < db->n_bam_rec; i++) {
        bam1_t* record = db->bam_rec[i];
        char* ref_name = core->m_hdr->target_name[record->core.tid];
        // printf("refname : %s\n",ref_name);
        int32_t ref_start_pos = record->core.pos;
        int32_t ref_end_pos = bam_endpos(record);
        assert(ref_end_pos >= ref_start_pos);

        // Extract the reference sequence for this region
        int32_t fetched_len = 0;
        t = realtime();
        char* refseq = faidx_fetch_seq(core->fai, ref_name, ref_start_pos, ref_end_pos, &fetched_len); // todo : error handle?
        core->db_fasta_time += realtime() - t;
        db->fasta_cache[i] = refseq;
        // printf("seq : %s\n",db->fasta_cache[i]);

        // get the fast5

        std::string qname = bam_get_qname(db->bam_rec[i]);
        t = realtime();
        std::string read_seq = core->readbb->get_read_sequence(qname);
        core->db_fasta_time += realtime() - t;

        //get the read in ascci
        db->read[i] =
            (char*)malloc(read_seq.size() + 1); // todo : is +1 needed? do errorcheck
        strcpy(db->read[i], read_seq.c_str());
        db->read_len[i] = strlen(db->read[i]);
        db->sum_bases += db->read_len[i];

        db->read_stat_flag[i] = 0; //reset the flag
    }
    // fprintf(stderr,"%s:: %d fast5 read\n",__func__,db->n_bam_rec);
    if(core->opt.verbosity>1){
        STDERR("Average read len %.0f",db->sum_bases/(float)db->n_bam_rec);
    }
    status.num_reads=db->n_bam_rec;
    assert(status.num_bases==db->sum_bases);

    double load_end = realtime();
    core->load_db_time += (load_end-load_start);

    return status;
}




void *get_fast5_from_pipe(void *voidargs){

    pthread_arg_t *args = (pthread_arg_t *)voidargs;
    core_t* core = args->core;
    db_t* db = args->db;
    int starti = args->starti;
    int endi = args->endi;
    int32_t t = args->thread_index;


    FILE *pipefp = core->pipefp_c2p[t];

    //read from the pipe
    int i;
    for(i = starti; i < endi; i++){

        if(core->opt.verbosity>2){
            STDERR("tid %d : %d reading from pipe",t,i);
        }
        db->f5[i] = (fast5_t*)calloc(1, sizeof(fast5_t));
        MALLOC_CHK(db->f5[i]);
        f5read(pipefp,&(db->f5[i]->nsample), sizeof(hsize_t), 1);
        //STDERR("tid %d : %d th read samples %lld",t,i,db->f5[i]->nsample);
        if(db->f5[i]->nsample>0){
            db->f5[i]->rawptr = (float*)calloc(db->f5[i]->nsample, sizeof(float));
            MALLOC_CHK( db->f5[i]->rawptr);
            f5read(pipefp,db->f5[i]->rawptr, sizeof(float), db->f5[i]->nsample);
            //STDERR("tid %d : %d th read samples read",t,i);
            f5read(pipefp,&(db->f5[i]->digitisation), sizeof(float), 1);
            f5read(pipefp,&(db->f5[i]->offset), sizeof(float), 1);
            f5read(pipefp,&(db->f5[i]->range), sizeof(float), 1);
            f5read(pipefp,&(db->f5[i]->sample_rate), sizeof(float), 1);
        }
        if(core->opt.verbosity>2){
            STDERR("tid %d : %d read from pipe",t,i);
        }

    }
    if(core->opt.verbosity>1){
        STDERR("%d-%d read from pipe",starti,endi);
    }

    pthread_exit(0);

}



void *request_fast5(void *voidargs){

    pthread_arg_t *args = (pthread_arg_t *)voidargs;
    core_t* core = args->core;
    db_t* db = args->db;
    int starti = args->starti;
    int endi = args->endi;
    int32_t t = args->thread_index;


    FILE *pipefp = core->pipefp_p2c[t];

    //write to pipe
    int i;
    for(i = starti; i < endi; i++){

        std::string qname = bam_get_qname(db->bam_rec[i]);
        std::string fast5_path_str = core->readbb->get_signal_path(qname);

        int ret = fprintf(pipefp,"%s\t%s\n",qname.c_str(),fast5_path_str.c_str());
        if(ret <=0){
            ERROR("%s","Could not write to pipe.");
            exit(EXIT_FAILURE);
        }
        //fprintf(stderr,"%d::%s\t%s\n",i,qname.c_str(),fast5_path_str.c_str());


    }
    if(core->opt.verbosity>1){
        STDERR("%d-%d written to pipe",starti,endi);
    }
    int ret=fflush(pipefp);
    if(ret !=0){
        ERROR("%s", "Flushing the pipe failed");
        exit(EXIT_FAILURE);
    }

    pthread_exit(0);

}

void fetch_fast5_multi_iop(core_t* core, db_t* db){


    int32_t num_io_proc = core->opt.num_iop;
    //INFO("Running with %d IO procs",num_io_proc);

    //create threads
    pthread_t tids_p2c[num_io_proc];
    pthread_t tids_c2p[num_io_proc];
    pthread_arg_t pt_args[num_io_proc];
    int32_t t, ret;
    int32_t i = 0;

    int32_t step = (db->n_bam_rec + num_io_proc - 1) / num_io_proc;
    //todo : check for higher num of threads than the data
    //current works but many threads are created despite

    //set the data structures
    for (t = 0; t < num_io_proc; t++) {
        pt_args[t].core = core;
        pt_args[t].db = db;
        pt_args[t].thread_index = t;
        pt_args[t].starti = i;
        i += step;
        if (i > db->n_bam_rec) {
            pt_args[t].endi = db->n_bam_rec;
        } else {
            pt_args[t].endi = i;
        }
        //pt_args[t].func=func;
        //fprintf(stderr,"t%d : %d-%d\n",t,pt_args[t].starti,pt_args[t].endi);

    }

    double rt;
    rt = realtime();
    //create threads
    for(t = 0; t < num_io_proc; t++){
        ret = pthread_create(&tids_p2c[t], NULL, request_fast5,
                                (void*)(&pt_args[t]));
        NEG_CHK(ret);
        ret = pthread_create(&tids_c2p[t], NULL, get_fast5_from_pipe,
                                (void*)(&pt_args[t]));
        NEG_CHK(ret);
    }
    if(core->opt.verbosity>1){
        STDERR("%s","created pipeline handling threads");
    }

    core->db_fast5_open_time += realtime() - rt;

    rt = realtime();
    //pthread joining
    for (t = 0; t < num_io_proc; t++) {
        int ret = pthread_join(tids_p2c[t], NULL);
        NEG_CHK(ret);
        ret = pthread_join(tids_c2p[t], NULL);
        NEG_CHK(ret);
    }
    core->db_fast5_read_time += realtime() - rt;
    if(core->opt.verbosity>1){
        STDERR("%s","joined threads");
    }
}


ret_status_t load_db2(core_t* core, db_t* db) { //separately load fast5 for multiple I/O procs

    double load_start = realtime();

    // get bams
    bam1_t* record;
    int32_t result = 0;
    db->n_bam_rec = 0;
    db->sum_bases = 0;
    db->total_reads = 0;
    db->bad_fast5_file = 0;
    db->ultra_long_skipped =0;

    ret_status_t status={0,0};
    int32_t i = 0;
    double t = 0;
    while (db->n_bam_rec < db->capacity_bam_rec && status.num_bases<core->opt.batch_size_bases) {
        i=db->n_bam_rec;
        record = db->bam_rec[i];
        t = realtime();
        result = sam_itr_next(core->m_bam_fh, core->itr, record);
        core->db_bam_time += realtime() - t;

        //set read index
        db->read_idx[i]= core->read_index;
        core->read_index +=1;

        if (result < 0) {
            break;
        } else {
            if ((record->core.flag & BAM_FUNMAP) == 0 &&
                record->core.qual >= core->opt.min_mapq) {
                // printf("%s\t%d\n",bam_get_qname(db->bam_rec[db->n_bam_rec]),result);

                if(!(core->opt.flag & F5C_SECONDARY_YES)){
                    if((record->core.flag & BAM_FSECONDARY)){
                        continue;
                    }
                }

                db->total_reads++; // candidate read

                std::string qname = bam_get_qname(record);
                t = realtime();
                //todo : make efficient (redudantly accessed below, can be combined with it?)
                int64_t read_length=core->readbb->get_read_sequence(qname).size();
                std::string fast5_path_str = core->readbb->get_signal_path(qname);
                core->db_fasta_time += realtime() - t;

                //skipping ultra-long-reads
                if(core->ultra_long_tmp!=NULL && read_length > core->opt.ultra_thresh){
                    db->ultra_long_skipped++;
                    int ret_wr=sam_write1(core->ultra_long_tmp,core->m_hdr,record);
                    NEG_CHK(ret_wr);
                    continue;
                }

                if(fast5_path_str==""){
                    handle_bad_fast5(core, db,fast5_path_str,qname);
                    continue;
                }

                db->n_bam_rec++;
                status.num_bases += read_length;
            }

        }
    }
    // fprintf(stderr,"%s:: %d queries read\n",__func__,db->n_bam_rec);

    // get ref sequences (todo can make efficient by taking the the start and end of the sorted bam)
    for (i = 0; i < db->n_bam_rec; i++) {
        bam1_t* record = db->bam_rec[i];
        char* ref_name = core->m_hdr->target_name[record->core.tid];
        // printf("refname : %s\n",ref_name);
        int32_t ref_start_pos = record->core.pos;
        int32_t ref_end_pos = bam_endpos(record);
        assert(ref_end_pos >= ref_start_pos);

        // Extract the reference sequence for this region
        int32_t fetched_len = 0;
        t = realtime();
        char* refseq = faidx_fetch_seq(core->fai, ref_name, ref_start_pos, ref_end_pos, &fetched_len); // todo : error handle?
        core->db_fasta_time += realtime() - t;
        db->fasta_cache[i] = refseq;
        // printf("seq : %s\n",db->fasta_cache[i]);

        // get the fast5

        std::string qname = bam_get_qname(db->bam_rec[i]);
        t = realtime();
        std::string read_seq = core->readbb->get_read_sequence(qname);
        core->db_fasta_time += realtime() - t;

        //get the read in ascci
        db->read[i] =
            (char*)malloc(read_seq.size() + 1); // todo : is +1 needed? do errorcheck
        strcpy(db->read[i], read_seq.c_str());
        db->read_len[i] = strlen(db->read[i]);
        db->sum_bases += db->read_len[i];

        db->read_stat_flag[i] = 0; //reset the flag
    }
    // fprintf(stderr,"%s:: %d fast5 read\n",__func__,db->n_bam_rec);
    if(core->opt.verbosity>1){
        STDERR("Average read len %.0f",db->sum_bases/(float)db->n_bam_rec);
    }
    status.num_reads=db->n_bam_rec;
    assert(status.num_bases==db->sum_bases);

    //read the fast5 batch
    t=realtime();
    fetch_fast5_multi_iop(core,db);
    core->db_fast5_time += realtime() - t;

    double load_end = realtime();
    core->load_db_time += (load_end-load_start);

    return status;
}

ret_status_t load_db(core_t* core, db_t* db) {
    if(core->opt.num_iop == 1){
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


void event_single(core_t* core,db_t* db, int32_t i) {

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

    // if(db->et[i].n/(float)db->read_len[i] > 20){
    //     fprintf(stderr,"%s\tevents_per_base\t%f\tread_len\t%d\n",bam_get_qname(db->bam_rec[i]), db->et[i].n/(float)db->read_len[i],db->read_len[i]);
    // }

    //get the scalings
    db->scalings[i] = estimate_scalings_using_mom(
        db->read[i], db->read_len[i], core->model, db->et[i]);

}

void event_db(core_t* core, db_t* db){

    if (core->opt.num_thread == 1) {
        int32_t i=0;
        for (i = 0; i < db->n_bam_rec; i++) {
            event_single(core,db,i);
        }

    }

    else {
        pthread_db(core,db,event_single);
    }

}



void scaling_single(core_t* core, db_t* db, int32_t i){

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
            db->read_stat_flag[i] |= FAILED_CALIBRATION;
            return;
        }

        free(db->event_alignment[i]);

    } else {
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

void scaling_db(core_t* core, db_t* db){
    if (core->opt.num_thread == 1) {
        int32_t i=0;
        for (i = 0; i < db->n_bam_rec; i++) {
            scaling_single(core,db,i);
        }

    }
    else {
        pthread_db(core,db,scaling_single);
    }
}

void align_single(core_t* core, db_t* db, int32_t i) {
    db->n_event_align_pairs[i] = align(
            db->event_align_pairs[i], db->read[i], db->read_len[i], db->et[i],
            core->model, db->scalings[i], db->f5[i]->sample_rate);
        //fprintf(stderr,"readlen %d,n_events %d\n",db->read_len[i],n_event_align_pairs);
}


void align_db(core_t* core, db_t* db) {
#ifdef HAVE_CUDA
    if (!(core->opt.flag & F5C_DISABLE_CUDA)) {
        //STDERR("%s","Performing on cuda");
        align_cuda(core, db);
    }
#endif

    if (core->opt.flag & F5C_DISABLE_CUDA) {
        //fprintf(stderr, "cpu\n");
        if (core->opt.num_thread == 1) {
            int i;
            for (i = 0; i < db->n_bam_rec; i++) {
                align_single(core, db, i);
            }
        } else {
            pthread_db(core, db, align_single);
        }
    }
}


void meth_single(core_t* core, db_t* db, int32_t i){
    if(!db->read_stat_flag[i]){
        if(core->mode==0){
            calculate_methylation_for_read(db->site_score_map[i], db->fasta_cache[i], db->bam_rec[i], db->read_len[i], db->et[i].event, db->base_to_event_map[i],
            db->scalings[i], core->cpgmodel,db->events_per_base[i]);
        }
        else if (core->mode==1){
            realign_read(db->event_alignment_result[i], &(db->eventalign_summary[i]),core->event_summary_fp, db->fasta_cache[i],core->m_hdr,
                  db->bam_rec[i],db->read_len[i],
                  i,
                  core->clip_start,
                  core->clip_end,
                  &(db->et[i]), core->model,db->base_to_event_map[i],db->scalings[i],db->events_per_base[i], db->f5[i]->sample_rate);
        }
    }
}

void meth_db(core_t* core, db_t* db) {
    if (core->opt.num_thread == 1) {
        int i;
        for (i = 0; i < db->n_bam_rec; i++) {
            meth_single(core, db, i);
        }
    }
    else {
        pthread_db(core, db, meth_single);
    }
}



void process_single(core_t* core, db_t* db,int32_t i) {

    event_single(core,db,i);

    db->event_align_pairs[i] = (AlignedPair*)malloc(
        sizeof(AlignedPair) * db->et[i].n * 2); //todo : find a good heuristic to save memory //todo : save memory by freeing here itself
    MALLOC_CHK(db->event_align_pairs[i]);

    align_single(core, db,i);

    db->event_alignment[i] = NULL;
    db->n_event_alignment[i] = 0;
    db->events_per_base[i] = 0; //todo : is double needed? not just float?

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
            db->read_stat_flag[i] |= FAILED_CALIBRATION;
            return;
        }

        free(db->event_alignment[i]);

    } else {
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

    if(core->mode==0){
        calculate_methylation_for_read(db->site_score_map[i], db->fasta_cache[i], db->bam_rec[i], db->read_len[i], db->et[i].event, db->base_to_event_map[i],
            db->scalings[i], core->cpgmodel,db->events_per_base[i]);
    }

    else if(core->mode==1){
        //hack
        realign_read(db->event_alignment_result[i], &(db->eventalign_summary[i]),core->event_summary_fp,db->fasta_cache[i],core->m_hdr,
                  db->bam_rec[i],db->read_len[i],
                  i,
                  core->clip_start,
                  core->clip_end,
                  &(db->et[i]), core->model,db->base_to_event_map[i],db->scalings[i],db->events_per_base[i],db->f5[i]->sample_rate);
    }
}

void process_db(core_t* core, db_t* db) {

    double process_start = realtime();

    if((core->opt.flag&F5C_SEC_PROF) || (!(core->opt.flag & F5C_DISABLE_CUDA))){

        double realtime0=core->realtime0;
        int32_t i;

        double event_start = realtime();
        event_db(core,db);
        double event_end = realtime();
        core->event_time += (event_end-event_start);

        fprintf(stderr, "[%s::%.3f*%.2f] Events computed\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0));

        for (i = 0; i < db->n_bam_rec; i++) {
            db->event_align_pairs[i] = (AlignedPair*)malloc(
                sizeof(AlignedPair) * db->et[i].n * 2); //todo : find a good heuristic to save memory
            MALLOC_CHK(db->event_align_pairs[i]);
        }

        double align_start = realtime();
        align_db(core, db);
        double align_end = realtime();
        core->align_time += (align_end-align_start);

        fprintf(stderr, "[%s::%.3f*%.2f] Banded alignment done\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0));

        double est_scale_start = realtime();
        scaling_db(core,db);
        double est_scale_end = realtime();
        core->est_scale_time += (est_scale_end-est_scale_start);

        fprintf(stderr, "[%s::%.3f*%.2f] Scaling calibration done\n", __func__,
                realtime() - realtime0, cputime() / (realtime() - realtime0));

        double meth_start = realtime();
        meth_db(core,db);
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

void output_db(core_t* core, db_t* db) {
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

    core->sum_bases += db->sum_bases;
    core->total_reads += db->total_reads;
    core->bad_fast5_file += db->bad_fast5_file;
    core->ultra_long_skipped += db->ultra_long_skipped;

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

                    // fprintf(stderr, "%s\t%d\t%d\t", ss.chromosome.c_str(), ss.start_position, ss.end_position);
                    // fprintf(stderr, "%s\t%.2lf\t", qname, diff);
                    // fprintf(stderr, "%.2lf\t%.2lf\t", sum_ll_m, sum_ll_u);
                    // fprintf(stderr, "%d\t%d\t%s\n", ss.strands_scored, ss.n_cpg, ss.sequence.c_str());

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
                    fprintf(summary_fp, "%ld\t%s\t", (long)(db->read_idx[i]), qname);
                    fprintf(summary_fp, "%s\t%s\t%s\t",fast5_path_str.c_str(), "dna", strand_idx == 0 ? "template" : "complement");
                    fprintf(summary_fp, "%d\t%d\t%d\t%d\t", summary.num_events, summary.num_steps, summary.num_skips, summary.num_stays);
                    fprintf(summary_fp, "%.2lf\t%.3lf\t%.3lf\t%.3lf\t%.3lf\n", summary.sum_duration/(db->f5[i]->sample_rate), scalings.shift, scalings.scale, 0.0, scalings.var);
                }
                std::vector<event_alignment_t> *event_alignment_result = db->event_alignment_result[i];
                int8_t print_read_names = (core->opt.flag & F5C_PRINT_RNAME) ? 1 : 0;
                int8_t scale_events = (core->opt.flag & F5C_SCALE_EVENTS) ? 1 : 0;
                int8_t write_samples = (core->opt.flag & F5C_PRINT_SAMPLES) ? 1 : 0;
                int8_t sam_output = (core->opt.flag & F5C_SAM) ? 1 : 0;

                if(sam_output==0){
                    emit_event_alignment_tsv(stdout,0,&(db->et[i]),core->model,db->scalings[i],*event_alignment_result, print_read_names, scale_events, write_samples,
                              db->read_idx[i], qname, contig, db->f5[i]->sample_rate);
                }
                else{
                    emit_event_alignment_sam(core->sam_output , qname, core->m_hdr, db->bam_rec[i], *event_alignment_result);
                }
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
    //core->read_index = core->read_index + db->n_bam_rec;

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
        delete db->site_score_map[i];
        db->site_score_map[i] = new std::map<int, ScoredSite>;

        if(db->event_alignment_result){ //eventalign related
            delete db->event_alignment_result[i];
            db->event_alignment_result[i] = new std::vector<event_alignment_t>;
        }
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
    free(db->read_idx);
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

    free(db);
}

void init_opt(opt_t* opt) {
    memset(opt, 0, sizeof(opt_t));
    opt->min_mapq = 20;
    opt->batch_size = 512;
    opt->batch_size_bases = 2*1000*1000;
    opt->num_thread = 8;
    opt->num_iop = 1;
    opt->region_str = NULL; //whole genome processing if null
#ifndef HAVE_CUDA
    opt->flag |= F5C_DISABLE_CUDA;
    opt->batch_size_bases = 5*1000*1000;
#endif

    opt->flag |= F5C_SKIP_UNREADABLE;
    opt->debug_break=-1;
    opt->ultra_thresh=100000;

    opt->meth_out_version=1;

    opt->cuda_block_size=64;
    opt->cuda_dev_id=0;
    opt->cuda_mem_frac=1.0f; //later set by cuda_init()

    //effective only if  CPU_GPU_PROC  is set
    opt->cuda_max_readlen=3.0f;
    opt->cuda_avg_events_per_kmer=2.0f; //only if CUDA_DYNAMIC_MALLOC is unset
    opt->cuda_max_avg_events_per_kmer=5.0f;
}


//CHANGE: function that sets the parameter values based on specified profile name or config file.
int set_profile(char *profile, opt_t *opt){
    //Load preset value if argument passed is the name of a machine for which there is a default profile
    if(strcmp(profile,"jetson-nano") == 0){
        set_opt_profile(opt,Nanojet);
    }else if(strcmp(profile,"jetson-tx2") == 0){
        set_opt_profile(opt,JetsonTx2);
    }else if(strcmp(profile,"jetson-xavier") == 0){
        set_opt_profile(opt,Xavier);
    }else{
        //Try to read from .profile file
        //profile name specifies a file from which to read values from.
        FILE *fptr = fopen(profile, "r");
        F_CHK(fptr,profile);
        int32_t batch_size, num_thread;
        int64_t batch_size_bases, ultra_thresh;
        float cuda_max_readlen,cuda_avg_events_per_kmer,cuda_max_avg_events_per_kmer;

        //read file and set parameter values
        int result = fscanf(fptr, "%f %f %f %d %" PRId64 " %d %" PRId64,
        &cuda_max_readlen,&cuda_avg_events_per_kmer,&cuda_max_avg_events_per_kmer,
        &batch_size,&batch_size_bases,&num_thread,&ultra_thresh);

        fprintf(stderr,"[%s] Profile loaded\n",__func__);
        fprintf(stderr,"[%s] batch_size: %d\nbatch_size_bases: %ld\nnum_thread: %d\nultra_thresh: %ld\ncuda_max_readlen: %f\ncuda_avg_events_per_kmer: %.2f\ncuda_max_avg_events_per_kmer: %.2f\n",
        __func__,batch_size,(long)batch_size_bases,num_thread,(long)ultra_thresh,cuda_max_readlen,cuda_avg_events_per_kmer,cuda_max_avg_events_per_kmer);

        if(result < 7){
            ERROR("%s","Malformed profile config file.");
            exit(EXIT_FAILURE);
        }

        set_opts(opt,batch_size,batch_size_bases,num_thread,ultra_thresh,cuda_max_readlen,
        cuda_avg_events_per_kmer,cuda_max_avg_events_per_kmer);
    }
    return 0;
}

//CHANGE: helper function to set profile to a preset value.
void set_opt_profile(opt_t *opt, parameters machine){
    opt->cuda_max_readlen = machine.cuda_max_readlen;
    opt->cuda_avg_events_per_kmer = machine.cuda_avg_events_per_kmer;
    opt->cuda_max_avg_events_per_kmer = machine.cuda_max_events_per_kmer;
    opt->batch_size = machine.batch_size;
    opt->batch_size_bases = machine.batch_size_bases;
    opt->num_thread = machine.num_thread;
    opt->ultra_thresh = machine.ultra_thresh;
}

//CHANGE: helper function to set user specified options. Pass -1 to corresponding arg if using default parameter value.
void set_opts(opt_t *opt, int32_t batch_size, int64_t batch_size_bases, int32_t num_thread, int64_t ultra_thresh, float cuda_max_readlen, float cuda_avg_events_per_kmer, float cuda_max_avg_events_per_kmer){
    if( batch_size != -1) opt->batch_size = batch_size;
    if( batch_size_bases != -1) opt->batch_size_bases = batch_size_bases;
    if( num_thread != -1) opt->num_thread = num_thread;
    if( ultra_thresh != -1) opt->ultra_thresh = ultra_thresh;
    if( cuda_max_readlen != -1) opt->cuda_max_readlen = cuda_max_readlen;
    if( cuda_avg_events_per_kmer != -1) opt->cuda_avg_events_per_kmer = cuda_avg_events_per_kmer;
    if( cuda_max_avg_events_per_kmer != -1) opt->cuda_max_avg_events_per_kmer = cuda_max_avg_events_per_kmer;
}
