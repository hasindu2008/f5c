/* @file f5cio.c
**
** f5c interface implementation - input related functions
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
#include "fast5lite.h"
#include "f5cmisc.h"

#include <sys/wait.h>
#include <unistd.h>

//in f5c.c
void pthread_db(core_t* core, db_t* db, void (*func)(core_t*,db_t*,int));

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
        signal_t *sig = (signal_t*)calloc(1, sizeof(signal_t));
        MALLOC_CHK(sig);

        //INFO("Reading %s",qname.c_str());
        int32_t ret=fast5_read(fast5_file, sig,qname);

        if(ret<0){
            ERROR("%s","Fast5 causes crashes with --iop");
            exit(EXIT_FAILURE);
        }

        fast5_close(fast5_file);
        success=1;

        //write to the pipe
        //STDERR("writing to pipe %s",qname_str);
        f5write(pipefp,&(sig->nsample), sizeof(uint64_t), 1);
        f5write(pipefp,sig->rawptr, sizeof(float), sig->nsample);
        f5write(pipefp,&(sig->digitisation), sizeof(float), 1);
        f5write(pipefp,&(sig->offset), sizeof(float), 1);
        f5write(pipefp,&(sig->range), sizeof(float), 1);
        f5write(pipefp,&(sig->sample_rate), sizeof(float), 1);

        ret=fflush(pipefp);
        if(ret!=0){
            ERROR("%s","Flushing the pipe failed");
            exit(EXIT_FAILURE);
        }

        free(sig->rawptr);
        free(sig);
        //STDERR("wrote to pipe %s : %lld samples",qname_str,sig->nsample);
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
    db->sig[i] = (signal_t*)calloc(1, sizeof(signal_t));
    MALLOC_CHK(db->sig[i]);

    f5read(core->raw_dump,&(db->sig[i]->nsample), sizeof(uint64_t), 1);

    if(db->sig[i]->nsample>0){
        db->sig[i]->rawptr = (float*)calloc(db->sig[i]->nsample, sizeof(float));
        MALLOC_CHK( db->sig[i]->rawptr);
        f5read(core->raw_dump,db->sig[i]->rawptr, sizeof(float), db->sig[i]->nsample);
        f5read(core->raw_dump,&(db->sig[i]->digitisation), sizeof(float), 1);
        f5read(core->raw_dump,&(db->sig[i]->offset), sizeof(float), 1);
        f5read(core->raw_dump,&(db->sig[i]->range), sizeof(float), 1);
        f5read(core->raw_dump,&(db->sig[i]->sample_rate), sizeof(float), 1);
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
        db->sig[i] = (signal_t*)calloc(1, sizeof(signal_t));
        MALLOC_CHK(db->sig[i]);
        t = realtime();
        int32_t ret=fast5_read(fast5_file, db->sig[i],qname);
        double rt = realtime() - t;
        core->db_fast5_read_time += rt;
        core->db_fast5_time += rt;
        if(ret<0){
            handle_bad_fast5(core, db,fast5_path,qname);
            if(core->opt.flag & F5C_WR_RAW_DUMP){
                uint64_t tmp_nsample = 0;
                f5write(core->raw_dump,&tmp_nsample, sizeof(uint64_t), 1);
            }
            free(fast5_path);
            return 0;
        }
        t = realtime();
        fast5_close(fast5_file);
        core->db_fast5_time += realtime() - t;

        if (core->opt.flag & F5C_PRINT_RAW) {
            printf(">%s\tPATH:%s\tLN:%lu\n", qname.c_str(), fast5_path,
                db->sig[i]->nsample);
            uint32_t j = 0;
            for (j = 0; j < db->sig[i]->nsample; j++) {
                printf("%d\t", (int)db->sig[i]->rawptr[j]);
            }
            printf("\n");
        }
        if(core->opt.flag & F5C_WR_RAW_DUMP){
            //write the fast5 dump to the binary file pointer core->raw_dump
            f5write(core->raw_dump,&(db->sig[i]->nsample), sizeof(uint64_t), 1);
            f5write(core->raw_dump,db->sig[i]->rawptr, sizeof(float), db->sig[i]->nsample);
            f5write(core->raw_dump,&(db->sig[i]->digitisation), sizeof(float), 1);
            f5write(core->raw_dump,&(db->sig[i]->offset), sizeof(float), 1);
            f5write(core->raw_dump,&(db->sig[i]->range), sizeof(float), 1);
            f5write(core->raw_dump,&(db->sig[i]->sample_rate), sizeof(float), 1);
        }

        //db->n_bam_rec++;
        //t = realtime();
        //status.num_bases += read_length;
        //core->db_fasta_time += realtime() - t;
        success=1;
    } else {
        handle_bad_fast5(core, db,fast5_path,qname);
        if(core->opt.flag & F5C_WR_RAW_DUMP){
            uint64_t tmp_nsample = 0;
            f5write(core->raw_dump,&tmp_nsample, sizeof(uint64_t), 1);
        }
        return 0;
    }
    free(fast5_path);
    assert(success==1);
    return 1;
}


/*************** Start of multiple I/O thread based SLOW5 reading *****************************/

//read a single read from a SLOW5 file
void read_slow5_single(core_t* core, db_t* db, int i){

    // int ret = 0;
    std::string qname = bam_get_qname(db->bam_rec[i]);

    db->sig[i] = (signal_t*)calloc(1, sizeof(signal_t));
    MALLOC_CHK(db->sig[i]);

    int len=0;

    slow5_rec_t *record=NULL;
    len = slow5_get(qname.c_str(), &record, core->sf);


    if(record==NULL || len <0){ //todo : should we free if len<0
        db->bad_fast5_file++;
        if (core->opt.flag & F5C_SKIP_UNREADABLE) {
            WARNING("Slow5 record for read [%s] is unavailable/unreadable and will be skipped", qname.c_str());
            db->sig[i]->nsample = 0;
            db->sig[i]->rawptr = NULL;
        } else {
            ERROR("Slow5 record for read [%s] is unavailable/unreadable", qname.c_str());
            exit(EXIT_FAILURE);
        }
        // ret=0;
    }
    else{
        assert(strcmp(qname.c_str(),record->read_id)==0);
        db->sig[i]->nsample = record->len_raw_signal;  //n_samples
        assert(db->sig[i]->nsample>0);
        db->sig[i]->rawptr = (float*)calloc(db->sig[i]->nsample, sizeof(float));
        MALLOC_CHK( db->sig[i]->rawptr);

        db->sig[i]->digitisation = (float)record->digitisation;
        db->sig[i]->offset = (float)record->offset;
        db->sig[i]->range = (float)record->range;
        db->sig[i]->sample_rate = (float)record->sampling_rate;

        for (int j = 0; j < (int)db->sig[i]->nsample; j++) { //check for int overflow
            db->sig[i]->rawptr[j] = (float)record->raw_signal[j];
        }

        // ret=1;
        slow5_rec_free(record);
    }

    // if(ret!=1){
    //     assert(0);
    // }
}


/******************************* SLOW5 end ************************************/

int f5c_sam_itr_next(core_t* core, bam1_t* record){

    if(core->reg_list==NULL){
        return(sam_itr_next(core->m_bam_fh, core->itr, record));
    }
    else{
        while(1){
            if(core->itr==NULL){ //if the iterator for the current region is NULL, go to the next region
                if(core->reg_i<core->reg_n){
                    if(core->opt.verbosity > 0){
                        STDERR("Iterating over region: %s", core->reg_list[core->reg_i]);
                    }
                    core->itr = sam_itr_querys(core->m_bam_idx, core->m_hdr, core->reg_list[core->reg_i]);
                    if(core->itr==NULL){
                        ERROR("sam_itr_querys failed. Please check if the region string in bed file [%s] is valid",core->reg_list[core->reg_i]);
                        exit(EXIT_FAILURE);
                    }
                    core->reg_i++;
                }
                else{
                    return -1; //no more regions left
                }
            }
            else{ //the current region still left
                int ret = sam_itr_next(core->m_bam_fh, core->itr, record);
                if(ret>=0){
                    return ret;
                }
                else{ //the region is done
                    sam_itr_destroy(core->itr);
                    core->itr=NULL;
                }
            }
        }
    }

    assert(0);

}


ret_status_t load_db1(core_t* core, db_t* db) { //no iop - used for slow5

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
        result = f5c_sam_itr_next(core, record);
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
                std::string fast5_path_str;
                t = realtime();
                //todo : make efficient (redudantly accessed below, can be combined with it?)
                int64_t read_length=core->readbb->get_read_sequence(qname).size();
                if(!(core->opt.flag & F5C_RD_SLOW5)){ //if not slow5
                    fast5_path_str = core->readbb->get_signal_path(qname);
                }
                core->db_fasta_time += realtime() - t;

                //skipping ultra-long-reads
                if(core->ultra_long_tmp!=NULL && read_length > core->opt.ultra_thresh){
                    db->ultra_long_skipped++;
                    int ret_wr=sam_write1(core->ultra_long_tmp,core->m_hdr,record);
                    NEG_CHK(ret_wr);
                    continue;
                }

                if(!(core->opt.flag & F5C_RD_SLOW5) && fast5_path_str==""){
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
                else if (core->opt.flag & F5C_RD_SLOW5){
                    read_status = 1; //we do this later with multiple threaded
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

        // get the read in ASCII
        std::string qname = bam_get_qname(db->bam_rec[i]);

        t = realtime();
        std::string read_seq = core->readbb->get_read_sequence(qname);
        core->db_fasta_time += realtime() - t;

        db->read[i] = (char*)malloc(read_seq.size() + 1); // todo : is +1 needed? do errorcheck
        strcpy(db->read[i], read_seq.c_str());
        db->read_len[i] = strlen(db->read[i]);
        if(core->opt.flag & F5C_RNA){
            replace_char(db->read[i], 'U', 'T');
        }
        db->sum_bases += db->read_len[i];

        db->read_stat_flag[i] = 0; //reset the flag
    }
    // fprintf(stderr,"%s:: %d fast5 read\n",__func__,db->n_bam_rec);
    if(core->opt.verbosity>1){
        STDERR("Average read len %.0f",db->sum_bases/(float)db->n_bam_rec);
    }
    status.num_reads=db->n_bam_rec;
    assert(status.num_bases==db->sum_bases);

    //read the slow5 batch
    if(core->opt.flag & F5C_RD_SLOW5){
        t=realtime();
        pthread_db(core, db, read_slow5_single);
        core->db_fast5_time += realtime() - t;
    }

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
        db->sig[i] = (signal_t*)calloc(1, sizeof(signal_t));
        MALLOC_CHK(db->sig[i]);
        f5read(pipefp,&(db->sig[i]->nsample), sizeof(uint64_t), 1);
        //STDERR("tid %d : %d th read samples %lld",t,i,db->sig[i]->nsample);
        if(db->sig[i]->nsample>0){
            db->sig[i]->rawptr = (float*)calloc(db->sig[i]->nsample, sizeof(float));
            MALLOC_CHK( db->sig[i]->rawptr);
            f5read(pipefp,db->sig[i]->rawptr, sizeof(float), db->sig[i]->nsample);
            //STDERR("tid %d : %d th read samples read",t,i);
            f5read(pipefp,&(db->sig[i]->digitisation), sizeof(float), 1);
            f5read(pipefp,&(db->sig[i]->offset), sizeof(float), 1);
            f5read(pipefp,&(db->sig[i]->range), sizeof(float), 1);
            f5read(pipefp,&(db->sig[i]->sample_rate), sizeof(float), 1);
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
        result = f5c_sam_itr_next(core, record);
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

        // get the read in ASCII
        std::string qname = bam_get_qname(db->bam_rec[i]);
        t = realtime();
        std::string read_seq = core->readbb->get_read_sequence(qname);
        core->db_fasta_time += realtime() - t;

        db->read[i] = (char*)malloc(read_seq.size() + 1); // todo : is +1 needed? do errorcheck
        strcpy(db->read[i], read_seq.c_str());
        db->read_len[i] = strlen(db->read[i]);
        if(core->opt.flag & F5C_RNA){
            replace_char(db->read[i], 'U', 'T');
        }
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
