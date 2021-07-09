/**
 * @file slow5_press.h
 * @brief SLOW5 compression and decompression functions
 * @author Sasha Jenner (jenner.sasha@gmail.com)
 * @date 27/02/2021
 */

/**************************************************************************************************
 ***  Low-level API *******************************************************************************
 **************************************************************************************************/

/*
IMPORTANT: The low-level API is not yet stable. Subject to changes in the future.
Function proptotypes can be changed without notice or completely removed
So do NOT use these functions in your code
these functions are used by slow5tools and pyslow5 - so any change to a function here means slow5tools and pyslow5 must be fixed
*/

/*
MIT License

Copyright (c) 2020 Sasha Jenner

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef SLOW5_PRESS_H
#define SLOW5_PRESS_H

#include <zlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
//#include "slow5_misc.h"

#ifdef __cplusplus
extern "C" {
#endif

// Compression methods
enum slow5_press_method {
    SLOW5_COMPRESS_NONE,
    SLOW5_COMPRESS_GZIP
};
typedef uint8_t slow5_press_method_t;

#define SLOW5_Z_MEM_DEFAULT (8)
#define SLOW5_Z_OUT_CHUNK (16384) // 2^14

// Gzip stream
struct slow5_gzip_stream {
    z_stream strm_inflate;
    z_stream strm_deflate;
    int flush;
};


// Compression streams
union slow5_press_stream {
    struct slow5_gzip_stream *gzip;
};

// Compression object
struct slow5_press {
    slow5_press_method_t method;
    union slow5_press_stream *stream;
};
typedef struct slow5_press slow5_press_t;

/* --- Init / free slow5_press structure --- */
struct slow5_press *slow5_press_init(slow5_press_method_t method);
void slow5_press_free(struct slow5_press *comp);

/* --- Compress / decompress a ptr to some memory --- */
void *slow5_ptr_compress(struct slow5_press *comp, const void *ptr, size_t count, size_t *n);
static inline void *slow5_str_compress(struct slow5_press *comp, const char *str, size_t *n) {
    return slow5_ptr_compress(comp, str, strlen(str) + 1, n); // Include '\0'
}
void *slow5_ptr_depress(struct slow5_press *comp, const void *ptr, size_t count, size_t *n);
void *slow5_ptr_depress_multi(slow5_press_method_t method, const void *ptr, size_t count, size_t *n);

/* --- Compress / decompress a ptr to some file --- */
size_t slow5_fwrite_compress(struct slow5_press *comp, const void *ptr, size_t size, size_t nmemb, FILE *fp);
size_t slow5_fwrite_depress(struct slow5_press *comp, const void *ptr, size_t size, size_t nmemb, FILE *fp); // TODO
static inline size_t slow5_print_compress(struct slow5_press *comp, const void *ptr, size_t size, size_t nmemb) {
    return slow5_fwrite_compress(comp, ptr, size, nmemb, stdout);
}
static inline size_t slow5_print_depress(struct slow5_press *comp, const void *ptr, size_t size, size_t nmemb) {
    return slow5_fwrite_depress(comp, ptr, size, nmemb, stdout);
}
static inline size_t slow5_fwrite_str_compress(struct slow5_press *comp, const char *str, FILE *fp) {
    return slow5_fwrite_compress(comp, str, sizeof *str, strlen(str), fp); // Don't include '\0'
}
static inline size_t slow5_print_str_compress(struct slow5_press *comp, const char *str) {
    return slow5_fwrite_str_compress(comp, str, stdout);
}

/* --- Decompress to a ptr from some file --- */
void *slow5_fread_depress(struct slow5_press *comp, size_t count, FILE *fp, size_t *n);
void *slow5_pread_depress(struct slow5_press *comp, int fd, size_t count, off_t offset, size_t *n);
void *slow5_pread_depress_multi(slow5_press_method_t method, int fd, size_t count, off_t offset, size_t *n);

/* --- Compress with format string to some file --- */
int slow5_fprintf_compress(struct slow5_press *comp, FILE *fp, const char *format, ...);
int slow5_printf_compress(struct slow5_press *comp, const char *format, ...);

/* --- Write compression footer on immediate next compression call --- */
void slow5_compress_footer_next(struct slow5_press *comp);

#ifdef __cplusplus
}
#endif

#endif
