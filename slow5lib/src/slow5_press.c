#define _XOPEN_SOURCE 700
#include <zlib.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <slow5/slow5.h>
#include <slow5/slow5_error.h>
#include <slow5/slow5_press.h>
#include <slow5/slow5_defs.h>
#include <streamvbyte.h>
#include <streamvbyte_zigzag.h>
#ifdef SLOW5_USE_ZSTD
#include <zstd.h>
#endif /* SLOW5_USE_ZSTD */
#include "slow5_misc.h"

extern enum slow5_log_level_opt  slow5_log_level;
extern enum slow5_exit_condition_opt  slow5_exit_condition;

/* zlib */
static int zlib_init_deflate(z_stream *strm);
static int zlib_init_inflate(z_stream *strm);
static void *ptr_compress_zlib(struct slow5_zlib_stream *zlib, const void *ptr, size_t count, size_t *n);
static void *ptr_compress_zlib_solo(const void *ptr, size_t count, size_t *n);
static void *ptr_depress_zlib(struct slow5_zlib_stream *zlib, const void *ptr, size_t count, size_t *n);
static void *ptr_depress_zlib_solo(const void *ptr, size_t count, size_t *n);
static ssize_t fwrite_compress_zlib(struct slow5_zlib_stream *zlib, const void *ptr, size_t size, size_t nmemb, FILE *fp);

/* streamvbyte */
static uint8_t *ptr_compress_svb(const uint32_t *ptr, size_t count, size_t *n);
static uint8_t *ptr_compress_svb_zd(const int16_t *ptr, size_t count, size_t *n);
static uint32_t *ptr_depress_svb(const uint8_t *ptr, size_t count, size_t *n);
static int16_t *ptr_depress_svb_zd(const uint8_t *ptr, size_t count, size_t *n);

#ifdef SLOW5_USE_ZSTD
/* zstd */
static void *ptr_compress_zstd(const void *ptr, size_t count, size_t *n);
static void *ptr_depress_zstd(const void *ptr, size_t count, size_t *n);
#endif /* SLOW5_USE_ZSTD */

/* other */
static int vfprintf_compress(struct __slow5_press *comp, FILE *fp, const char *format, va_list ap);


/* Convert the press method back an forth from the library format to spec format
   Implementation is horrible, but for now.
 */

// convert the record compression method from library format to the spec format
uint8_t slow5_encode_record_press(enum slow5_press_method method){
    uint8_t ret = 0;
    switch(method){
        case SLOW5_COMPRESS_NONE:
            ret = 0;
            break;
        case SLOW5_COMPRESS_ZLIB:
            ret = 1;
            break;
        case SLOW5_COMPRESS_ZSTD:
            ret = 2;
            break;
        case SLOW5_COMPRESS_SVB_ZD:  //hidden feature hack for devs
            SLOW5_WARNING("You are using a hidden dev features (record compression in %s). Output files may be useless.","svb-zd");
            ret = 250;
            break;
        default:  //todo: Proper error?
            ret = 255;
            SLOW5_WARNING("Unknown record compression method %d",method);
            break;
    }
    return ret;
}

// convert the record compression method from spec format to the library format
enum slow5_press_method slow5_decode_record_press(uint8_t method){
    enum slow5_press_method ret = SLOW5_COMPRESS_NONE;
    switch(method){
        case 0:
            ret = SLOW5_COMPRESS_NONE;
            break;
        case 1:
            ret = SLOW5_COMPRESS_ZLIB;
            break;
        case 2:
            ret = SLOW5_COMPRESS_ZSTD;
            break;
        case 250:  //hidden feature hack for devs
            ret = SLOW5_COMPRESS_SVB_ZD;
            break;
        default:    //todo Proper error
            ret = 255;
            SLOW5_WARNING("Unknown record compression method %d",method);
            break;
    }
    return ret;
}

// convert the signal compression from library format to the spec format
uint8_t slow5_encode_signal_press(enum slow5_press_method method){
    uint8_t ret = 0;
    switch(method){
        case SLOW5_COMPRESS_NONE:
            ret = 0;
            break;
        case SLOW5_COMPRESS_SVB_ZD:
            ret = 1;
            break;
        case SLOW5_COMPRESS_ZLIB: //hidden feature hack for devs
            SLOW5_WARNING("You are using a hidden dev features (signal compression in %s). Output files may be useless.", "zlib");
            ret = 250;
            break;
        case SLOW5_COMPRESS_ZSTD: //hidden feature hack for devs
            SLOW5_WARNING("You are using a hidden dev features (signal compression in %s). Output files may be useless.","zstd");
            ret = 251;
            break;
        default:    //todo: Proper error?
            ret = 255;
            SLOW5_WARNING("Unknown signal compression method %d",method);
            break;
    }
    return ret;
}

// convert the signal compression from spec format to library format
enum slow5_press_method slow5_decode_signal_press(uint8_t method){
    enum slow5_press_method ret = 0;
    switch(method){
        case 0:
            ret = SLOW5_COMPRESS_NONE;
            break;
        case 1:
            ret = SLOW5_COMPRESS_SVB_ZD;
            break;
        case 250: //hidden feature hack for devs
            ret = SLOW5_COMPRESS_ZLIB;
            break;
        case 251: //hidden feature hack for devs
            ret = SLOW5_COMPRESS_ZSTD;
            break;
        default: //todo: Proper error?
            ret = 255;
            SLOW5_WARNING("Unknown signal compression method %d",method);
            break;
    }
    return ret;
}


/* --- Init / free (__)slow5_press structures --- */

/*
 * init slow5_press struct
 * returns NULL on error and sets errno
 * errors
 * SLOW5_ERR_MEM
 * SLOW5_ERR_ARG    method is bad
 * SLOW5_ERR_PRESS  (de)compression failure
 */
struct slow5_press *slow5_press_init(slow5_press_method_t method) {

    struct __slow5_press *record_comp = __slow5_press_init(method.record_method);
    if (!record_comp) {
        return NULL;
    }

    struct __slow5_press *signal_comp = __slow5_press_init(method.signal_method);
    if (!signal_comp) {
        __slow5_press_free(record_comp);
        return NULL;
    }

    struct slow5_press *comp = (struct slow5_press *) calloc(1, sizeof *comp);
    if (!comp) {
        SLOW5_MALLOC_ERROR();
        __slow5_press_free(record_comp);
        __slow5_press_free(signal_comp);
        slow5_errno = SLOW5_ERR_MEM;
        return NULL;
    }
    comp->record_press = record_comp;
    comp->signal_press = signal_comp;

    return comp;
}

void slow5_press_free(struct slow5_press *comp) {

    if (comp) {
        __slow5_press_free(comp->record_press);
        __slow5_press_free(comp->signal_press);
        free(comp);
    }
}

/*
 * init __slow5_press struct
 * returns NULL on error and sets errno
 * errors
 * SLOW5_ERR_MEM
 * SLOW5_ERR_ARG    method is bad
 * SLOW5_ERR_PRESS  (de)compression failure
 */
struct __slow5_press *__slow5_press_init(enum slow5_press_method method) {

    struct __slow5_press *comp = NULL;

    comp = (struct __slow5_press *) calloc(1, sizeof *comp);
    if (!comp) {
        SLOW5_MALLOC_ERROR();
        slow5_errno = SLOW5_ERR_MEM;
        return NULL;
    }
    comp->method = method;

    switch (method) {
        case SLOW5_COMPRESS_NONE: break;
        case SLOW5_COMPRESS_ZLIB:
            {
            struct slow5_zlib_stream *zlib;

            zlib = (struct slow5_zlib_stream *) malloc(sizeof *zlib);
            if (!zlib) {
                SLOW5_MALLOC_ERROR();
                free(comp);
                slow5_errno = SLOW5_ERR_MEM;
                return NULL;
            }

            if (zlib_init_deflate(&(zlib->strm_deflate)) != Z_OK) {
                SLOW5_ERROR("zlib deflate init failed: %s.", zlib->strm_deflate.msg);
                free(zlib);
                free(comp);
                slow5_errno = SLOW5_ERR_PRESS;
                return NULL;
            }
            if (zlib_init_inflate(&(zlib->strm_inflate)) != Z_OK) {
                SLOW5_ERROR("zlib inflate init failed: %s.", zlib->strm_inflate.msg);
                if (deflateEnd(&(zlib->strm_deflate)) != Z_OK) {
                    SLOW5_ERROR("zlib deflate end failed: %s.", zlib->strm_deflate.msg);
                }
                free(zlib);
                free(comp);
                slow5_errno = SLOW5_ERR_PRESS;
                return NULL;
            }

            zlib->flush = Z_NO_FLUSH;
            comp->stream = (union slow5_press_stream *) malloc(sizeof *comp->stream);
            if (!comp->stream) {
                SLOW5_MALLOC_ERROR();
                if (deflateEnd(&(zlib->strm_deflate)) != Z_OK) {
                    SLOW5_ERROR("zlib deflate end failed: %s.", zlib->strm_deflate.msg);
                }
                if (inflateEnd(&(zlib->strm_inflate)) != Z_OK) {
                    SLOW5_ERROR("zlib inflate end failed: %s.", zlib->strm_inflate.msg);
                }
                free(zlib);
                free(comp);
                slow5_errno = SLOW5_ERR_PRESS;
                return NULL;
            }

            comp->stream->zlib = zlib;
            } break;

        case SLOW5_COMPRESS_SVB_ZD: break;
        case SLOW5_COMPRESS_ZSTD:
#ifdef SLOW5_USE_ZSTD
            break;
#else
            SLOW5_ERROR("%s","slow5lib has not been compiled with zstd support to read/write zstd compressed BLOW5 files.");
            free(comp);
            slow5_errno = SLOW5_ERR_ARG;
            return NULL;
#endif /* SLOW5_USE_ZSTD */

        default:
            SLOW5_ERROR("Invalid or unsupported (de)compression method '%d'.", method);
            free(comp);
            slow5_errno = SLOW5_ERR_ARG;
            return NULL;
    }

    return comp;
}

void __slow5_press_free(struct __slow5_press *comp) {

    if (comp) {

        switch (comp->method) {

            case SLOW5_COMPRESS_NONE: break;
            case SLOW5_COMPRESS_ZLIB:
                (void) deflateEnd(&(comp->stream->zlib->strm_deflate));
                (void) inflateEnd(&(comp->stream->zlib->strm_inflate));
                free(comp->stream->zlib);
                free(comp->stream);
                break;
            case SLOW5_COMPRESS_SVB_ZD: break;
#ifdef SLOW5_USE_ZSTD
            case SLOW5_COMPRESS_ZSTD: break;
#endif /* SLOW5_USE_ZSTD */

            default:
                SLOW5_ERROR("Invalid or unsupported (de)compression method '%d'.", comp->method);
                slow5_errno = SLOW5_ERR_ARG;
                break;
        }

        free(comp);
    }
}

void *slow5_ptr_compress_solo(enum slow5_press_method method, const void *ptr, size_t count, size_t *n) {
    void *out = NULL;
    size_t n_tmp = 0;

    if (!ptr) {
        SLOW5_ERROR("Argument '%s' cannot be NULL.", SLOW5_TO_STR(ptr))
        slow5_errno = SLOW5_ERR_ARG;
    } else {
        switch (method) {

            case SLOW5_COMPRESS_NONE:
                out = (void *) malloc(count);
                SLOW5_MALLOC_CHK(out);
                if (!out) {
                    // Malloc failed
                    return out;
                }
                memcpy(out, ptr, count);
                n_tmp = count;
                break;

            case SLOW5_COMPRESS_ZLIB:
                out = ptr_compress_zlib_solo(ptr, count, &n_tmp);
                break;

            case SLOW5_COMPRESS_SVB_ZD:
                out = ptr_compress_svb_zd(ptr, count, &n_tmp);
                break;

#ifdef SLOW5_USE_ZSTD
            case SLOW5_COMPRESS_ZSTD:
                out = ptr_compress_zstd(ptr, count, &n_tmp);
                break;
#endif /* SLOW5_USE_ZSTD */

            default:
                SLOW5_ERROR("Invalid or unsupported (de)compression method '%d'.", method);
                slow5_errno = SLOW5_ERR_ARG;
                break;
        }
    }

    if (n) {
        *n = n_tmp;
    }

    return out;
}

/* --- Compress / decompress to some memory --- */

void *slow5_ptr_compress(struct __slow5_press *comp, const void *ptr, size_t count, size_t *n) {
    void *out = NULL;
    size_t n_tmp = 0;

    if (comp && ptr) {

        switch (comp->method) {

            case SLOW5_COMPRESS_NONE:
                out = (void *) malloc(count);
                SLOW5_MALLOC_CHK(out);
                if (!out) {
                    // Malloc failed
                    return out;
                }
                memcpy(out, ptr, count);
                n_tmp = count;
                break;

            case SLOW5_COMPRESS_ZLIB:
                if (comp->stream && comp->stream->zlib) {
                    out = ptr_compress_zlib(comp->stream->zlib, ptr, count, &n_tmp);
                }
                break;

            case SLOW5_COMPRESS_SVB_ZD:
                out = ptr_compress_svb_zd(ptr, count, &n_tmp);
                break;

#ifdef SLOW5_USE_ZSTD
            case SLOW5_COMPRESS_ZSTD:
                out = ptr_compress_zstd(ptr, count, &n_tmp);
                break;
#endif /* SLOW5_USE_ZSTD */

            default:
                SLOW5_ERROR("Invalid or unsupported (de)compression method '%d'.", comp->method);
                slow5_errno = SLOW5_ERR_ARG;
                break;
        }
    }

    if (n) {
        *n = n_tmp;
    }

    return out;
}

/*
 * ptr cannot be NULL
 */
void *slow5_ptr_depress_solo(enum slow5_press_method method, const void *ptr, size_t count, size_t *n) {
    void *out = NULL;
    size_t n_tmp = 0;

    if (!ptr) {
        SLOW5_ERROR("Argument '%s' cannot be NULL.", SLOW5_TO_STR(ptr))
        slow5_errno = SLOW5_ERR_ARG;
    } else {

        switch (method) {

            case SLOW5_COMPRESS_NONE:
                out = (void *) malloc(count);
                SLOW5_MALLOC_CHK(out);
                if (!out) {
                    /* Malloc failed */
                    return out;
                }
                memcpy(out, ptr, count);
                n_tmp = count;
                break;

            case SLOW5_COMPRESS_ZLIB:
                out = ptr_depress_zlib_solo(ptr, count, &n_tmp);
                break;

            case SLOW5_COMPRESS_SVB_ZD:
                out = ptr_depress_svb_zd(ptr, count, &n_tmp);
                break;

#ifdef SLOW5_USE_ZSTD
            case SLOW5_COMPRESS_ZSTD:
                out = ptr_depress_zstd(ptr, count, &n_tmp);
                break;
#endif /* SLOW5_USE_ZSTD */

            default:
                SLOW5_ERROR("Invalid or unsupported (de)compression method '%d'.", method);
                slow5_errno = SLOW5_ERR_ARG;
                break;
        }
    }

    if (n) {
        *n = n_tmp;
    }

    return out;
}

/*
 * decompress count bytes of a ptr to compressed memory
 * returns pointer to decompressed memory of size *n bytes to be later freed
 * returns NULL on error and *n set to 0
 * DONE
 */
void *slow5_ptr_depress(struct __slow5_press *comp, const void *ptr, size_t count, size_t *n) {
    void *out = NULL;

    if (!comp || !ptr) {
        if (!comp) {
            SLOW5_ERROR("Argument '%s' cannot be NULL.", SLOW5_TO_STR(comp));
        }
        if (!ptr) {
            SLOW5_ERROR("Argument '%s' cannot be NULL.", SLOW5_TO_STR(ptr));
        }
        if (n) {
            *n = 0;
        }
        slow5_errno = SLOW5_ERR_ARG;
        return NULL;
    } else {
        size_t n_tmp = 0;

        switch (comp->method) {

            case SLOW5_COMPRESS_NONE:
                out = (void *) malloc(count);
                SLOW5_MALLOC_CHK(out);
                if (!out) { /* malloc failed */
                    if (n) {
                        *n = 0;
                    }
                    return out;
                }
                memcpy(out, ptr, count);
                if (n) {
                    *n = count;
                }
                break;

            case SLOW5_COMPRESS_ZLIB:
                if (!comp->stream) {
                    SLOW5_ERROR("%s", "Decompression stream cannot be NULL.")
                } else {
                    out = ptr_depress_zlib(comp->stream->zlib, ptr, count, &n_tmp);
                    if (!out) {
                        SLOW5_ERROR("%s", "zlib decompression failed.")
                    }
                }
                break;

            case SLOW5_COMPRESS_SVB_ZD:
                out = ptr_depress_svb_zd(ptr, count, &n_tmp);
                break;

#ifdef SLOW5_USE_ZSTD
            case SLOW5_COMPRESS_ZSTD:
                out = ptr_depress_zstd(ptr, count, &n_tmp);
                break;
#endif /* SLOW5_USE_ZSTD */

            default:
                SLOW5_ERROR("Invalid or unsupported (de)compression method '%d'.", comp->method);
                slow5_errno = SLOW5_ERR_ARG;
                break;
        }

        if (n) {
            *n = n_tmp;
        }
    }

    return out;
}

/* --- Compress / decompress a ptr to some file --- */

/* return -1 on failure */
ssize_t slow5_fwrite_compress(struct __slow5_press *comp, const void *ptr, size_t size, size_t nmemb, FILE *fp) {
    ssize_t bytes = -1;
    size_t bytes_tmp = 0;
    void *out = NULL;

    if (comp) {
        switch (comp->method) {

            case SLOW5_COMPRESS_NONE:
                bytes = fwrite(ptr, size, nmemb, fp);
                if (bytes != size * nmemb || ferror(fp)) {
                    if (bytes != size * nmemb) {
                        SLOW5_ERROR("Expected to write '%zu' bytes, instead wrote '%zu' bytes.",
                                size * nmemb, bytes);
                    } else {
                        SLOW5_ERROR("%s", "File error after trying to write.");
                    }
                    slow5_errno = SLOW5_ERR_IO;
                    return -1;
                }
                break;

            case SLOW5_COMPRESS_ZLIB:
                if (comp->stream && comp->stream->zlib) {
                    bytes = fwrite_compress_zlib(comp->stream->zlib, ptr, size, nmemb, fp);
                }
                break;

            case SLOW5_COMPRESS_SVB_ZD:
                out = ptr_compress_svb_zd(ptr, size * nmemb, &bytes_tmp);
                if (!out) {
                    return -1;
                }
                break;

#ifdef SLOW5_USE_ZSTD
            case SLOW5_COMPRESS_ZSTD:
                out = ptr_compress_zstd(ptr, size * nmemb, &bytes_tmp);
                if (!out) {
                    return -1;
                }
                break;
#endif /* SLOW5_USE_ZSTD */

            default:
                SLOW5_ERROR("Invalid or unsupported (de)compression method '%d'.", comp->method);
                slow5_errno = SLOW5_ERR_ARG;
                return -1;
        }
    }

    if (out) {
        size_t bytes_written = fwrite(out, 1, bytes_tmp, fp);
        free(out);
        if (bytes_written != bytes_tmp || ferror(fp)) {
            if (bytes_written != bytes_tmp) {
                SLOW5_ERROR("Expected to write '%zu' compressed bytes, instead wrote '%zu' bytes.",
                        bytes_tmp, bytes_written);
            } else {
                SLOW5_ERROR("%s", "File error after trying to write.");
            }
            slow5_errno = SLOW5_ERR_IO;
            return -1;
        }
    }

    return bytes = (ssize_t) bytes_tmp;
}

/*
 * decompress count bytes from some file
 * returns pointer to decompressed memory of size *n bytes to be later freed
 * returns NULL on error and *n set to 0
 * DONE
 */
void *slow5_fread_depress(struct __slow5_press *comp, size_t count, FILE *fp, size_t *n) {
    void *raw = (void *) malloc(count);
    SLOW5_MALLOC_CHK(raw);
    if (!raw) {
        return NULL;
    }

    if (fread(raw, count, 1, fp) != 1) {
        SLOW5_ERROR("Failed to read '%zu' bytes from file.", count);
        free(raw);
        return NULL;
    }

    void *out = slow5_ptr_depress(comp, raw, count, n);
    if (!out) {
        SLOW5_ERROR("%s", "Decompression failed.")
    }
    free(raw);

    return out;
}
/*
static void *fread_depress_solo(slow5_press_method_t method, size_t count, FILE *fp, size_t *n) {
    void *raw = (void *) malloc(count);
    SLOW5_MALLOC_CHK(raw);

    if (fread(raw, count, 1, fp) != 1) {
        free(raw);
        return NULL;
    }

    void *out = slow5_ptr_depress_solo(method, raw, count, n);
    free(raw);

    return out;
}
*/
void *slow5_pread_depress(struct __slow5_press *comp, int fd, size_t count, off_t offset, size_t *n) {
    void *raw = (void *) malloc(count);
    SLOW5_MALLOC_CHK(raw);

    if (pread(fd, raw, count, offset) == -1) {
        free(raw);
        return NULL;
    }

    void *out = slow5_ptr_depress(comp, raw, count, n);
    free(raw);

    return out;
}

/*
 * pread size bytes from offset then decompress the data according to the method
 */
void *slow5_pread_depress_solo(enum slow5_press_method method, int fd, size_t count, off_t offset, size_t *n) {
    void *raw = (void *) malloc(count);
    SLOW5_MALLOC_CHK(raw);
    if (!raw) {
        return NULL;
    }

    ssize_t ret;
    if ((ret = pread(fd, raw, count, offset)) != count) {
        if (ret == -1) {
            SLOW5_ERROR("pread failed to read '%zu' bytes: %s", count, strerror(errno));
        } else if (ret == 0) {
            SLOW5_ERROR("End of file reached. pread failed to read '%zu' bytes.", count);
        } else {
            SLOW5_ERROR("pread read less bytes '%zd' than expected '%zu'.", ret, count);
        }
        free(raw);
        return NULL;
    }

    void *out = slow5_ptr_depress_solo(method, raw, count, n);
    free(raw);

    return out;
}


/* --- Compress with printf to some file --- */

int slow5_fprintf_compress(struct __slow5_press *comp, FILE *fp, const char *format, ...) {
    int ret = -1;

    va_list ap;
    va_start(ap, format);

    ret = vfprintf_compress(comp, fp, format, ap);

    va_end(ap);

    return ret;
}

int slow5_printf_compress(struct __slow5_press *comp, const char *format, ...) {
    int ret = -1;

    va_list ap;
    va_start(ap, format);

    ret = vfprintf_compress(comp, stdout, format, ap);

    va_end(ap);

    return ret;
}

static int vfprintf_compress(struct __slow5_press *comp, FILE *fp, const char *format, va_list ap) {
    int ret = -1;

    if (comp) {

        if (comp->method == SLOW5_COMPRESS_NONE) {
            ret = vfprintf(fp, format, ap);
        } else {
            char *buf;
            if (slow5_vasprintf(&buf, format, ap) != -1) {
                ret = slow5_fwrite_str_compress(comp, buf, fp);
                free(buf);
            }
        }

    }

    return ret;
}


/* --- Write compression footer on immediate next compression call --- */

void slow5_compress_footer_next(struct __slow5_press *comp) {

    if (comp && comp->stream) {
        switch (comp->method) {
            case SLOW5_COMPRESS_NONE: break;
            case SLOW5_COMPRESS_ZLIB: {
                struct slow5_zlib_stream *zlib = comp->stream->zlib;

                if (zlib) {
                    zlib->flush = Z_FINISH;
                }
            } break;
            case SLOW5_COMPRESS_SVB_ZD: break;
#ifdef SLOW5_USE_ZSTD
            case SLOW5_COMPRESS_ZSTD: break;
#endif /* SLOW5_USE_ZSTD */

            default:
                SLOW5_ERROR("Invalid or unsupported (de)compression method '%d'.", comp->method);
                slow5_errno = SLOW5_ERR_ARG;
                break;
        }
    }
}




/********
 * ZLIB *
 ********/

static int zlib_init_deflate(z_stream *strm) {
    strm->zalloc = Z_NULL;
    strm->zfree = Z_NULL;
    strm->opaque = Z_NULL;

    return deflateInit2(strm,
            Z_DEFAULT_COMPRESSION,
            Z_DEFLATED,
            MAX_WBITS,
            SLOW5_ZLIB_MEM_DEFAULT,
            Z_DEFAULT_STRATEGY);
}

static int zlib_init_inflate(z_stream *strm) {
    strm->zalloc = Z_NULL;
    strm->zfree = Z_NULL;
    strm->opaque = Z_NULL;

    return inflateInit2(strm, MAX_WBITS);
}

static void *ptr_compress_zlib(struct slow5_zlib_stream *zlib, const void *ptr, size_t count, size_t *n) {
    uint8_t *out = NULL;

    size_t n_cur = 0;
    z_stream *strm = &(zlib->strm_deflate);

    strm->avail_in = count;
    strm->next_in = (Bytef *) ptr;

    uLong chunk_sz = SLOW5_ZLIB_COMPRESS_CHUNK;

    do {
        out = (uint8_t *) realloc(out, n_cur + chunk_sz);
        SLOW5_MALLOC_CHK(out);

        strm->avail_out = chunk_sz;
        strm->next_out = out + n_cur;

        if (deflate(strm, zlib->flush) == Z_STREAM_ERROR) {
            free(out);
            out = NULL;
            n_cur = 0;
            break;
        }

        n_cur += chunk_sz - strm->avail_out;

    } while (strm->avail_out == 0);

    *n = n_cur;

    if (zlib->flush == Z_FINISH) {
        zlib->flush = Z_NO_FLUSH;
        deflateReset(strm);
    }

    return out;
}

static void *ptr_compress_zlib_solo(const void *ptr, size_t count, size_t *n) {
    uint8_t *out = NULL;

    size_t n_cur = 0;

    z_stream strm_local;
    zlib_init_deflate(&strm_local);
    z_stream *strm = &strm_local;

    strm->avail_in = count;
    strm->next_in = (Bytef *) ptr;

    uLong chunk_sz = SLOW5_ZLIB_COMPRESS_CHUNK;

    do {
        out = (uint8_t *) realloc(out, n_cur + chunk_sz);
        SLOW5_MALLOC_CHK(out);

        strm->avail_out = chunk_sz;
        strm->next_out = out + n_cur;

        if (deflate(strm, Z_FINISH) == Z_STREAM_ERROR) {
            free(out);
            out = NULL;
            n_cur = 0;
            break;
        }

        n_cur += chunk_sz - strm->avail_out;

    } while (strm->avail_out == 0);

    *n = n_cur;

    (void) inflateEnd(strm);

    return out;
}

/*
 * zlib decompress count bytes of a ptr to compressed memory
 * returns pointer to decompressed memory of size *n bytes to be later freed
 * returns NULL on error and *n set to 0
 * DONE
 */
static void *ptr_depress_zlib(struct slow5_zlib_stream *zlib, const void *ptr, size_t count, size_t *n) {
    uint8_t *out = NULL;

    size_t n_cur = 0;
    if (!zlib) {
        SLOW5_ERROR("%s", "zlib stream cannot be NULL.")
        return NULL;
    }
    z_stream *strm = &(zlib->strm_inflate);
    if (!strm) {
        SLOW5_ERROR("%s", "zlib inflate stream cannot be NULL.")
        return NULL;
    }

    strm->avail_in = count;
    strm->next_in = (Bytef *) ptr;

    do {
        uint8_t *out_new = (uint8_t *) realloc(out, n_cur + SLOW5_ZLIB_DEPRESS_CHUNK);
        SLOW5_MALLOC_CHK(out_new);
        if (!out_new) {
            free(out);
            out = NULL;
            n_cur = 0;
            break;
        }
        out = out_new;

        strm->avail_out = SLOW5_ZLIB_DEPRESS_CHUNK;
        strm->next_out = out + n_cur;

        int ret = inflate(strm, Z_NO_FLUSH);
        if (ret == Z_STREAM_ERROR || ret == Z_DATA_ERROR) {
            SLOW5_ERROR("zlib inflate failed with error code %d.", ret)
            free(out);
            out = NULL;
            n_cur = 0;
            break;
        }

        n_cur += SLOW5_ZLIB_DEPRESS_CHUNK- strm->avail_out;

    } while (strm->avail_out == 0);

    *n = n_cur;
    if (out && inflateReset(strm) == Z_STREAM_ERROR) {
        SLOW5_WARNING("%s", "Stream state is inconsistent.")
    };

    return out;
}

static void *ptr_depress_zlib_solo(const void *ptr, size_t count, size_t *n) {
    uint8_t *out = NULL;

    size_t n_cur = 0;

    z_stream strm_local;
    zlib_init_inflate(&strm_local);
    z_stream *strm = &strm_local;

    strm->avail_in = count;
    strm->next_in = (Bytef *) ptr;

    do {
        out = (uint8_t *) realloc(out, n_cur + SLOW5_ZLIB_DEPRESS_CHUNK);
        SLOW5_MALLOC_CHK(out);

        strm->avail_out = SLOW5_ZLIB_DEPRESS_CHUNK;
        strm->next_out = out + n_cur;

        int ret = inflate(strm, Z_NO_FLUSH);
        if (ret == Z_STREAM_ERROR || ret == Z_DATA_ERROR || ret == Z_NEED_DICT || ret == Z_MEM_ERROR) {
            SLOW5_ERROR("%s","inflate failed");
            free(out);
            out = NULL;
            n_cur = 0;
            break;
        }

        n_cur += SLOW5_ZLIB_DEPRESS_CHUNK- strm->avail_out;

    } while (strm->avail_out == 0);

    *n = n_cur;

    (void) inflateEnd(strm);

    return out;
}

static ssize_t fwrite_compress_zlib(struct slow5_zlib_stream *zlib, const void *ptr, size_t size, size_t nmemb, FILE *fp) {

    ssize_t bytes = 0;
    z_stream *strm = &(zlib->strm_deflate);

    strm->avail_in = size * nmemb;
    strm->next_in = (Bytef *) ptr;

    uLong chunk_sz = SLOW5_ZLIB_COMPRESS_CHUNK;
    uint8_t *buf = (uint8_t *) malloc(sizeof *buf * chunk_sz);
    SLOW5_MALLOC_CHK(buf);
    if (!buf) {
        return -1;
    }

    do {
        strm->avail_out = chunk_sz;
        strm->next_out = buf;

        if (deflate(strm, zlib->flush) == Z_STREAM_ERROR) {
            bytes = -1;
            break;
        }

        size_t have = (sizeof *buf * chunk_sz) - strm->avail_out;
        if (fwrite(buf, sizeof *buf, have, fp) != have || ferror(fp)) {
            bytes = -1;
            break;
        }

        bytes += have;

    } while (strm->avail_out == 0);

    free(buf);

    if (zlib->flush == Z_FINISH) {
        zlib->flush = Z_NO_FLUSH;
    }

    return bytes;
}



/***************
 * STREAMVBYTE *
 ***************/

/* return NULL on malloc error, n cannot be NULL */
static uint8_t *ptr_compress_svb(const uint32_t *ptr, size_t count, size_t *n) {
    uint32_t length = count / sizeof *ptr;

    size_t max_n = __slow5_streamvbyte_max_compressedbytes(length);
    uint8_t *out = (uint8_t *) malloc(max_n + sizeof length);
    if (!out) {
        SLOW5_MALLOC_ERROR();
        slow5_errno = SLOW5_ERR_MEM;
        return NULL;
    }

    *n = __slow5_streamvbyte_encode(ptr, length, out + sizeof length);
    memcpy(out, &length, sizeof length); /* copy original length of ptr (needed for depress) */
    *n = *n + sizeof length;
    SLOW5_LOG_DEBUG("max svb bytes=%zu\nsvb bytes=%zu\n",
            max_n, *n); /* TESTING */
    return out;
}

/* return NULL on malloc error, n cannot be NULL */
static uint8_t *ptr_compress_svb_zd(const int16_t *ptr, size_t count, size_t *n) {
    uint32_t length = count / sizeof *ptr;
    int32_t *in = (int32_t *) malloc(length * sizeof *in);
    if (!in) {
        SLOW5_MALLOC_ERROR();
        slow5_errno = SLOW5_ERR_MEM;
        return NULL;
    }
    for (int64_t i = 0; i < length; ++ i) {
        in[i] = ptr[i];
    }

    uint32_t *diff = (uint32_t *) malloc(length * sizeof *diff);
    if (!diff) {
        SLOW5_MALLOC_ERROR();
        free(in);
        slow5_errno = SLOW5_ERR_MEM;
        return NULL;
    }
    __slow5_zigzag_delta_encode(in, diff, length, 0);


    SLOW5_LOG_DEBUG("orig bytes=%zu\n", count); /* TESTING */
    uint8_t *out = ptr_compress_svb(diff, length * sizeof *diff, n);

    free(in);
    free(diff);
    return out;
}

/* return NULL on malloc error, n cannot be NULL */
static uint32_t *ptr_depress_svb(const uint8_t *ptr, size_t count, size_t *n) {
    uint32_t length;
    memcpy(&length, ptr, sizeof length); /* get original array length */

    uint32_t *out = (uint32_t *) malloc(length * sizeof *out);
    if (!out) {
        SLOW5_MALLOC_ERROR();
        slow5_errno = SLOW5_ERR_MEM;
        return NULL;
    }

    size_t bytes_read;
    if ((bytes_read = __slow5_streamvbyte_decode(ptr + sizeof length, out, length)) != count - sizeof length) {
        SLOW5_ERROR("Expected streamvbyte_decode to read '%zu' bytes, instead read '%zu' bytes.",
                count - sizeof length, bytes_read);
        slow5_errno = SLOW5_ERR_PRESS;
        free(out);
        return NULL;
    }

    *n = length * sizeof *out;
    return out;
}

/* return NULL on malloc error, n cannot be NULL */
static int16_t *ptr_depress_svb_zd(const uint8_t *ptr, size_t count, size_t *n) {

    uint32_t *diff = ptr_depress_svb(ptr, count, n);
    if (!diff) {
        return NULL;
    }
    uint32_t length = *n / sizeof *diff;

    int16_t *orig = (int16_t *) malloc(length * sizeof *orig);
    if (!orig) {
        SLOW5_MALLOC_ERROR();
        free(diff);
        slow5_errno = SLOW5_ERR_MEM;
        return NULL;
    }
    __slow5_zigzag_delta_decode(diff, orig, length, 0);

    // int16_t *orig = (int16_t *) malloc(length * sizeof *orig);
    // for (int64_t i = 0; i < length; ++ i) {
    //     orig[i] = out[i];
    // }

    *n = length * sizeof *orig;
    free(diff);
    //free(out);
    return orig;
}



#ifdef SLOW5_USE_ZSTD
/********
 * ZSTD *
 ********/

/* return NULL on error */
static void *ptr_compress_zstd(const void *ptr, size_t count, size_t *n) {
    size_t max_bytes = ZSTD_compressBound(count);

    void *out = malloc(max_bytes);
    if (!out) {
        SLOW5_MALLOC_ERROR();
        slow5_errno = SLOW5_ERR_MEM;
        return NULL;
    }

    *n = ZSTD_compress(out, max_bytes, ptr, count, SLOW5_ZSTD_COMPRESS_LEVEL);
    if (ZSTD_isError(*n)) {
        SLOW5_ERROR("zstd compress failed with error code %zu.", *n);
        free(out);
        slow5_errno = SLOW5_ERR_PRESS;
        return NULL;
    }

    return out;
}

/* return NULL on error */
static void *ptr_depress_zstd(const void *ptr, size_t count, size_t *n) {
    unsigned long long depress_bytes = ZSTD_getFrameContentSize(ptr, count);
    if (depress_bytes == ZSTD_CONTENTSIZE_UNKNOWN ||
            depress_bytes == ZSTD_CONTENTSIZE_ERROR) {
        SLOW5_ERROR("zstd get decompressed size failed with error code %llu\n", depress_bytes);
        slow5_errno = SLOW5_ERR_PRESS;
        return NULL;
    }

    void *out = malloc(depress_bytes);
    if (!out) {
        SLOW5_MALLOC_ERROR();
        slow5_errno = SLOW5_ERR_MEM;
        return NULL;
    }

    *n = ZSTD_decompress(out, depress_bytes, ptr, count);
    if (ZSTD_isError(*n)) {
        SLOW5_ERROR("zstd decompress failed with error code %zu.", *n);
        free(out);
        slow5_errno = SLOW5_ERR_PRESS;
        return NULL;
    }

    return out;
}
#endif /* SLOW5_USE_ZSTD */



/* Decompress a zlib-compressed string
 *
 * @param       compressed string
 * @param       ptr to size of compressed string, updated to size of returned malloced string
 * @return      malloced string
 */
/*
unsigned char *z_inflate_buf(const char *comp_str, size_t *n) {

    z_stream strm;
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    strm.avail_in = *n;
    strm.next_in = (Bytef *) comp_str;

    *n = 0;

    uLong prev_sz = 0;
    uLong out_sz = 16328;
    unsigned char *out = (unsigned char *) malloc(sizeof *out * out_sz);
    SLOW5_MALLOC_CHK(out);

    int ret = inflateInit2(&strm, ZLIB_WBITS);

    if (ret != Z_OK) {
        free(out);
        return NULL;
    }

    do {
        strm.avail_out = out_sz;
        strm.next_out = out + prev_sz;

        ret = inflate(&strm, Z_NO_FLUSH);
        SLOW5_ASSERT(ret != Z_STREAM_ERROR);

        switch (ret) {
            case Z_NEED_DICT:
                ret = Z_DATA_ERROR;
            case Z_DATA_ERROR:
            case Z_MEM_ERROR:
                free(out);
                (void) inflateEnd(&strm);
                return NULL;
        }


        unsigned have = out_sz - strm.avail_out;
        prev_sz += have;

        if (strm.avail_out == 0) {
            out = (unsigned char *) realloc(out, sizeof *out * (prev_sz + out_sz));
            SLOW5_MALLOC_CHK(out);
        }

    } while (strm.avail_out == 0);

    *n = prev_sz;
    (void) inflateEnd(&strm);

    return out;
}

size_t z_deflate_buf(z_streamp strmp, const void *ptr, uLong size, FILE *f_out, int flush, int *err) {
unsigned char *z_inflate_buf(const char *comp_str, size_t *n) {
    size_t written = 0;

    strmp->avail_in = size;
    strmp->next_in = (Bytef *) ptr;

    uLong out_sz = SLOW5_ZLIB_COMPRESS_CHUNK;
    unsigned char *out = (unsigned char *) malloc(sizeof *out * out_sz);
    SLOW5_MALLOC_CHK(out);

    do {
        strmp->avail_out = out_sz;
        strmp->next_out = out;

        ret = deflate(strmp, flush);
        if (ret == Z_STREAM_ERROR) {
            ERROR("deflate failed\n%s", ""); // testing
            *err = ret;
            return written;
        }

        unsigned have = out_sz - strmp->avail_out;
        size_t tmp;
        if ((tmp = fwrite(out, sizeof *out, have, f_out)) != have || ferror(f_out)) {
            ERROR("fwrite\n%s", ""); // testing
            *err = Z_ERRNO;
            written += tmp * sizeof *out;
            return written;
        }
        written += tmp * sizeof *out;

    } while (strmp->avail_out == 0);

    free(out);
    out = NULL;

    // If still input to deflate
    if (strmp->avail_in != 0) {
        ERROR("still more input to deflate\n%s", ""); // testing
        *err = Z_ERRNO;
    }

    *err = Z_OK;
    return written;
}
*/
