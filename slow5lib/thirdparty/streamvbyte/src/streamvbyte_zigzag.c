#include "streamvbyte_zigzag.h"

static inline
uint32_t _zigzag_encode_32 (int32_t val) {
	return (val + val) ^ (val >> 31);
}

// NOTE: appending __slow5_
void __slow5_zigzag_encode(const int32_t * in, uint32_t * out, size_t N) {
    for(size_t i = 0; i < N; i++)
      out[i] = _zigzag_encode_32(in[i]);
}

// NOTE: appending __slow5_
void __slow5_zigzag_delta_encode(const int32_t * in, uint32_t * out, size_t N, int32_t prev) {
    for (size_t i = 0; i < N; i++) {
      out[i] = _zigzag_encode_32(in[i] - prev);
      prev = in[i];
    }
}

static inline
int32_t _zigzag_decode_32 (uint32_t val) {
	return (val >> 1) ^ -(val & 1);
}

// NOTE: appending __slow5_
void __slow5_zigzag_decode(const uint32_t * in, int32_t * out, size_t N) {
    for(size_t i = 0; i < N; i++)
      out[i] = _zigzag_decode_32(in[i]);
}

// NOTE: appending __slow5_
void __slow5_zigzag_delta_decode(const uint32_t * in, int16_t * out, size_t N, int32_t prev) {
    for(size_t i = 0; i < N; i++) {
      int32_t val =_zigzag_decode_32(in[i]);
      out[i] = val + prev;
      prev += val;
    }
}
