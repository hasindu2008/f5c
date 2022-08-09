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
IMPORTANT: The following low-level APIfunctions  are not yet finalised or documented.
If anyone is interested, please open a GitHub issue, rather than trying to figure out from the code.
Function proptotypes can be changed without notice or completely removed. So do NOT use these functions in your code.
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

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* zlib macros */
#define SLOW5_ZLIB_MEM_DEFAULT (8)
#define SLOW5_ZLIB_COMPRESS_CHUNK (131072) /* compression buffer size -- 128KB (2^17) */
#define SLOW5_ZLIB_DEPRESS_CHUNK (262144) /* decompression buffer size -- 256KB (2^18) */

#define SLOW5_ZSTD_COMPRESS_LEVEL (1)

/* (de)compression methods */
enum slow5_press_method {
    SLOW5_COMPRESS_NONE,
    SLOW5_COMPRESS_ZLIB,
    SLOW5_COMPRESS_SVB_ZD, /* streamvbyte zigzag delta */
    SLOW5_COMPRESS_ZSTD,
};
typedef struct{
    enum slow5_press_method record_method;
    enum slow5_press_method signal_method;
} slow5_press_method_t;

/* zlib stream */
struct slow5_zlib_stream {
    z_stream strm_inflate;
    z_stream strm_deflate;
    int flush;
};

/* (de)compression streams */
union slow5_press_stream {
    struct slow5_zlib_stream *zlib;
};

/* (de)compression object */
struct __slow5_press {
    enum slow5_press_method method;
    union slow5_press_stream *stream;
};

typedef struct slow5_press {
    struct __slow5_press *record_press;
    struct __slow5_press *signal_press;
} slow5_press_t;

/* init or free for multiple (de)compress calls */
struct slow5_press *slow5_press_init(slow5_press_method_t method);
struct __slow5_press *__slow5_press_init(enum slow5_press_method method);
void slow5_press_free(struct slow5_press *comp);
void __slow5_press_free(struct __slow5_press *comp);

/* (de)compress ptr */
void *slow5_ptr_compress(struct __slow5_press *comp, const void *ptr, size_t count, size_t *n);
void *slow5_ptr_compress_solo(enum slow5_press_method method, const void *ptr, size_t count, size_t *n);
void *slow5_ptr_depress(struct __slow5_press *comp, const void *ptr, size_t count, size_t *n);
void *slow5_ptr_depress_solo(enum slow5_press_method method, const void *ptr, size_t count, size_t *n);
static inline void *slow5_str_compress(struct __slow5_press *comp, const char *str, size_t *n);

/* (de)compress ptr and write */
ssize_t slow5_fwrite_compress(struct __slow5_press *comp, const void *ptr, size_t size, size_t nmemb, FILE *fp);
ssize_t slow5_fwrite_depress(struct __slow5_press *comp, const void *ptr, size_t size, size_t nmemb, FILE *fp); /* TODO */
static inline ssize_t slow5_fwrite_str_compress(struct __slow5_press *comp, const char *str, FILE *fp);
static inline ssize_t slow5_print_compress(struct __slow5_press *comp, const void *ptr, size_t size, size_t nmemb);
static inline ssize_t slow5_print_depress(struct __slow5_press *comp, const void *ptr, size_t size, size_t nmemb);
static inline ssize_t slow5_print_str_compress(struct __slow5_press *comp, const char *str);
int slow5_fprintf_compress(struct __slow5_press *comp, FILE *fp, const char *format, ...);
int slow5_printf_compress(struct __slow5_press *comp, const char *format, ...);

/* read and (de)compress */
void *slow5_fread_depress(struct __slow5_press *comp, size_t count, FILE *fp, size_t *n);
void *slow5_pread_depress(struct __slow5_press *comp, int fd, size_t count, off_t offset, size_t *n);
void *slow5_pread_depress_solo(enum slow5_press_method method, int fd, size_t count, off_t offset, size_t *n);

/* write compression footer on next compress call */
void slow5_compress_footer_next(struct __slow5_press *comp);

/* static inline function definitions */
static inline void *slow5_str_compress(struct __slow5_press *comp, const char *str, size_t *n) {
    return slow5_ptr_compress(comp, str, strlen(str) + 1, n); /* include '\0' */
}
static inline ssize_t slow5_fwrite_str_compress(struct __slow5_press *comp, const char *str, FILE *fp) {
    return slow5_fwrite_compress(comp, str, sizeof *str, strlen(str), fp); /* don't include '\0' */
}
static inline ssize_t slow5_print_compress(struct __slow5_press *comp, const void *ptr, size_t size, size_t nmemb) {
    return slow5_fwrite_compress(comp, ptr, size, nmemb, stdout);
}
static inline ssize_t slow5_print_depress(struct __slow5_press *comp, const void *ptr, size_t size, size_t nmemb) {
    return slow5_fwrite_depress(comp, ptr, size, nmemb, stdout);
}
static inline ssize_t slow5_print_str_compress(struct __slow5_press *comp, const char *str) {
    return slow5_fwrite_str_compress(comp, str, stdout);
}

uint8_t slow5_encode_record_press(enum slow5_press_method method);
enum slow5_press_method slow5_decode_record_press(uint8_t method);
uint8_t slow5_encode_signal_press(enum slow5_press_method method);
enum slow5_press_method slow5_decode_signal_press(uint8_t method);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* slow5_press.h */
