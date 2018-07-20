
/*********************** Some error checks *********************/
/*Die on error. Print the error and exit if the return value of the previous function NULL*/
#define errorCheckNULL(ret) ({\
    if (ret==NULL){ \
        fprintf(stderr,"Error at File %s line number %d : %s\n",__FILE__, __LINE__,strerror(errno));\
        exit(EXIT_FAILURE);\
    }\
    })

/*Die on error. Print the error and exit if the return value of the previous function is -1*/
#define errorCheck(ret) ({\
    if (ret<0){ \
        fprintf(stderr,"Error at File %s line number %d : %s\n",__FILE__, __LINE__,strerror(errno));\
        exit(EXIT_FAILURE);\
    }\
    })


struct data_batch_t{


	//region string
	char *region;

	//bam records
	bam1_t** bam_rec;
	int capacity_bam_rec;
	int n_bam_rec;

	//fasta cache
	char** fasta_cache;


	//fast5 cache
	fast5** f5;
};

typedef struct data_batch_t data_batch;

struct data_global_t{
	//bam file related
	htsFile* m_bam_fh;
    hts_idx_t* m_bam_idx;
    bam_hdr_t* m_hdr;
    hts_itr_t* itr;

    //fa related
    faidx_t *fai;

    //readbb
    ReadDB *readbb;
};


typedef struct data_global_t data_global;

