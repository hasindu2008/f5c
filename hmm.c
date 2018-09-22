float profile_hmm_score(const char *m_seq,const char *m_rc_seq, event_t* event, scalings_t scaling,  model_t* cpgmodel, uint32_t event_start_idx,
    uint32_t event_stop_idx,
    uint8_t strand,
    int8_t event_stride,
    uint8_t rc
);