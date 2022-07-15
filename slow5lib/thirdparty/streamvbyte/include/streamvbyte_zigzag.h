// NOTE: This is a modified version of the original Makefile.
// appended with slow5 to prevent nane clashes

#ifndef SLOW5_INCLUDE_STREAMVBYTE_ZIGZAG_H_
#define SLOW5_INCLUDE_STREAMVBYTE_ZIGZAG_H_
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <stdint.h>// please use a C99-compatible compiler
#include <stddef.h>

#if defined(__cplusplus)
extern "C" {
#endif

// NOTE: __slow5_ is appended to each function to prevent name clashes with the original library

/**
 * Convert N signed integers to N unsigned integers, using zigzag
 * encoding.
 */
void __slow5_zigzag_encode(const int32_t * in, uint32_t * out, size_t N);

/**
 * Convert N signed integers to N unsigned integers, using zigzag
 * delta encoding.
 */
void __slow5_zigzag_delta_encode(const int32_t * in, uint32_t * out, size_t N, int32_t prev);

/**
 * Convert N unsigned integers to N signed integers, using zigzag
 * encoding.
 */
void __slow5_zigzag_decode(const uint32_t * in, int32_t * out, size_t N);

/**
 * Convert N unsigned integers to N signed integers, using zigzag
 * delta encoding.
 */
void __slow5_zigzag_delta_decode(const uint32_t * in, int16_t * out, size_t N, int32_t prev);


#if defined(__cplusplus)
};
#endif

#endif /* INCLUDE_STREAMVBYTE_ZIGZAG_H_ */
