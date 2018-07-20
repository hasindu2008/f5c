/*********************** Some value definitions *********************/

//We do some static memory allocations. Some definitions for maximum sizes. Increase them if they are not adequate
#define MAX_READNAME_LEN 100  //The maximun size of a read name (qname size in bytes +1 )
#define MAX_READ_LEN 150      //maximum size of a read (number of bases +1)   
#define MAX_N_CIGAR 16        //no idea what this number of CIGAR ops mean at the moment  


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
    

/**************************** The data structure that stores a read ******************/
    
struct alignedRead {
    
  char qname[MAX_READNAME_LEN];         //Query template NAME (the name of the read)
  uint32_t flag;                        //bitwise FLAG
  int32_t chromID;                      //References sequence ID (Note that this is not the chromosome name. The chromosome name that maps to this ID must be found through the BAM header, See printRead function in common.c to see how to do this)
  uint32_t pos;                         //0-based leftmost mapping POSition (note that this is not 1-based coordinates)
  uint8_t mapq;                         //MAPping Quality 
  uint32_t cigarOps[2*MAX_N_CIGAR];     //CIGAR ops            
  uint32_t mateChromID;                 //Reference ID of the mate/next read (Note that this is not the chromosome name. The chromosome name that maps to this ID must be found through the BAM header, See printRead function in common.c to see how to do this)
  uint32_t matePos;                     //Position of the mate/next read (0-based)
  uint32_t tlen;                        //observed Template LENgth (I just later set this to 0)
  char seq[MAX_READ_LEN];               //segment SEQuence (the actual read sequence)
  uint8_t qual[MAX_READ_LEN];           //quality string   

  //some other dtuff
  uint32_t cigarLen;
  uint32_t rlen;                        //Length of SEQuence
  uint32_t end;
  uint32_t insertSize;
 
  //need to get some other fields in the BAM

};    


/* The function that gets a read to the alignedRead structure 
    Arguments    : The destination structure for the read to be stored, and the bam1_t pointer from htslib (See sequentialaccess.c for example usage)
    Return value : The same pointer of the input argument theRead */
struct alignedRead* getRead(struct alignedRead*  theRead, bam1_t *b/*, int storeRgID, char** rgID*/);


/*A function that prints a read to the stdout. First call getRead function and then this*/
void printRead(struct alignedRead* theRead, bam_hdr_t *header);
   

//an internally used function   
char _getBase(uint8_t *s, int i);   