#include "f5c.h"
#include "f5cmisc.h"

//#define INPUT_DEBUG 1

float profile_hmm_score(const char *m_seq,const char *m_rc_seq, event_t* event, scalings_t scaling,  model_t* cpgmodel, uint32_t event_start_idx,
    uint32_t event_stop_idx,
    uint8_t strand,
    int8_t event_stride,
    uint8_t rc, double events_per_base,uint32_t hmm_flags 
){

#ifdef INPUT_DEBUG

    fprintf(stderr,"m_seq : %s\n",m_seq);
    fprintf(stderr,"m_rc_seq : %s\n",m_rc_seq);
    fprintf(stderr,"event_start_idx %d, event_stop_idx %d, event_stride %d, rc %d\n",event_start_idx,event_stop_idx,event_stride,rc);
    fprintf(stderr,"events_per_base %f\n",events_per_base);
#endif

    return 0;


}