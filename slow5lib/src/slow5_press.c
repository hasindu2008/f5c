#define _XOPEN_SOURCE 700
#include <zlib.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <slow5/slow5_error.h>
#include <slow5/slow5_press.h>
#include "slow5_misc.h"

extern enum slow5_log_level_opt  slow5_log_level;
extern enum slow5_exit_condition_opt  slow5_exit_condition;

static int gzip_init_deflate(z_stream *strm);
static int gzip_init_inflate(z_stream *strm);

static void *ptr_compress_gzip(struct slow5_gzip_stream *gzip, const void *ptr, size_t count, size_t *n);
static void *ptr_depress_gzip(struct slow5_gzip_stream *gzip, const void *ptr, size_t count, size_t *n);
static void *ptr_depress_gzip_multi(const void *ptr, size_t count, size_t *n);
static size_t fwrite_compress_gzip(struct slow5_gzip_stream *gzip, const void *ptr, size_t size, size_t nmemb, FILE *fp);
static int vfprintf_compress(struct slow5_press *comp, FILE *fp, const char *format, va_list ap);

/* --- Init / free slow5_press structure --- */

struct slow5_press *slow5_press_init(slow5_press_method_t method) {

    struct slow5_press *comp = NULL;

    comp = (struct slow5_press *) malloc(sizeof *comp);
    comp->method = method;

    switch (method) {

        case SLOW5_COMPRESS_NONE:
            comp->stream = NULL;
            break;

        case SLOW5_COMPRESS_GZIP: {
            struct slow5_gzip_stream *gzip;

            gzip = (struct slow5_gzip_stream *) malloc(sizeof *gzip);
            if (gzip_init_deflate(&(gzip->strm_deflate)) != Z_OK) {
                free(gzip);
                comp->stream = NULL;
                SLOW5_WARNING("%s","zlib deflate initialisation failed.");
            } else if (gzip_init_inflate(&(gzip->strm_inflate)) != Z_OK) {
                (void) deflateEnd(&(gzip->strm_deflate));
                free(gzip);
                comp->stream = NULL;
                SLOW5_WARNING("%s","zlib inflate initialisation failed.");
            } else {
                gzip->flush = Z_NO_FLUSH;
                comp->stream = (union slow5_press_stream *) malloc(sizeof *comp->stream);
                comp->stream->gzip = gzip;
            }

        } break;
    }

    return comp;
}

static int gzip_init_deflate(z_stream *strm) {
    strm->zalloc = Z_NULL;
    strm->zfree = Z_NULL;
    strm->opaque = Z_NULL;

    return deflateInit2(strm,
            Z_DEFAULT_COMPRESSION,
            Z_DEFLATED,
            MAX_WBITS,
            SLOW5_Z_MEM_DEFAULT,
            Z_DEFAULT_STRATEGY);
}

static int gzip_init_inflate(z_stream *strm) {
    strm->zalloc = Z_NULL;
    strm->zfree = Z_NULL;
    strm->opaque = Z_NULL;

    return inflateInit2(strm, MAX_WBITS);
}

void slow5_press_free(struct slow5_press *comp) {

    if (comp != NULL) {

        switch (comp->method) {

            case SLOW5_COMPRESS_NONE:
                break;

            case SLOW5_COMPRESS_GZIP: {
                (void) deflateEnd(&(comp->stream->gzip->strm_deflate));
                (void) inflateEnd(&(comp->stream->gzip->strm_inflate));
                free(comp->stream->gzip);
                free(comp->stream);
            } break;
        }

        free(comp);
    }
}


/* --- Compress / decompress to some memory --- */

void *slow5_ptr_compress(struct slow5_press *comp, const void *ptr, size_t count, size_t *n) {
    void *out = NULL;
    size_t n_tmp = 0;

    if (comp != NULL && ptr != NULL) {

        switch (comp->method) {

            case SLOW5_COMPRESS_NONE:
                out = (void *) malloc(count);
                if (out == NULL) {
                    // Malloc failed
                    return out;
                }
                memcpy(out, ptr, count);
                n_tmp = count;
                break;

            case SLOW5_COMPRESS_GZIP:
                if (comp->stream != NULL && comp->stream->gzip != NULL) {
                    out = ptr_compress_gzip(comp->stream->gzip, ptr, count, &n_tmp);
                }
                break;
        }
    }

    if (n != NULL) {
        *n = n_tmp;
    }

    return out;
}

static void *ptr_compress_gzip(struct slow5_gzip_stream *gzip, const void *ptr, size_t count, size_t *n) {
    uint8_t *out = NULL;

    size_t n_cur = 0;
    z_stream *strm = &(gzip->strm_deflate);

    strm->avail_in = count;
    strm->next_in = (Bytef *) ptr;

    uLong chunk_sz = SLOW5_Z_OUT_CHUNK;

    do {
        out = (uint8_t *) realloc(out, n_cur + chunk_sz);

        strm->avail_out = chunk_sz;
        strm->next_out = out + n_cur;

        if (deflate(strm, gzip->flush) == Z_STREAM_ERROR) {
            free(out);
            out = NULL;
            n_cur = 0;
            break;
        }

        n_cur += chunk_sz - strm->avail_out;

    } while (strm->avail_out == 0);

    *n = n_cur;

    if (gzip->flush == Z_FINISH) {
        gzip->flush = Z_NO_FLUSH;
        deflateReset(strm);
    }

    return out;
}

void *slow5_ptr_depress_multi(slow5_press_method_t method, const void *ptr, size_t count, size_t *n) {
    void *out = NULL;
    size_t n_tmp = 0;

    if (ptr != NULL) {

        switch (method) {

            case SLOW5_COMPRESS_NONE:
                out = (void *) malloc(count);
                if (out == NULL) {
                    // Malloc failed
                    return out;
                }
                memcpy(out, ptr, count);
                n_tmp = count;
                break;

            case SLOW5_COMPRESS_GZIP:
                out = ptr_depress_gzip_multi(ptr, count, &n_tmp);
                break;
        }
    }
    else{
        SLOW5_WARNING("%s","slow5_ptr_depress_multi got a NULL ptr");
    }

    if (n != NULL) {
        *n = n_tmp;
    }

    return out;
}

void *slow5_ptr_depress(struct slow5_press *comp, const void *ptr, size_t count, size_t *n) {
    void *out = NULL;
    size_t n_tmp = 0;

    if (comp != NULL && ptr != NULL) {

        switch (comp->method) {

            case SLOW5_COMPRESS_NONE:
                out = (void *) malloc(count);
                if (out == NULL) {
                    // Malloc failed
                    return out;
                }
                memcpy(out, ptr, count);
                n_tmp = count;
                break;

            case SLOW5_COMPRESS_GZIP:
                if (comp->stream != NULL && comp->stream->gzip != NULL) {
                    out = ptr_depress_gzip(comp->stream->gzip, ptr, count, &n_tmp);
                }
                break;
        }
    }

    if (n != NULL) {
        *n = n_tmp;
    }

    return out;
}

static void *ptr_depress_gzip(struct slow5_gzip_stream *gzip, const void *ptr, size_t count, size_t *n) {
    uint8_t *out = NULL;

    size_t n_cur = 0;
    z_stream *strm = &(gzip->strm_inflate);

    strm->avail_in = count;
    strm->next_in = (Bytef *) ptr;

    do {
        out = (uint8_t *) realloc(out, n_cur + SLOW5_Z_OUT_CHUNK);

        strm->avail_out = SLOW5_Z_OUT_CHUNK;
        strm->next_out = out + n_cur;

        int ret = inflate(strm, Z_NO_FLUSH);
        if (ret == Z_STREAM_ERROR || ret == Z_DATA_ERROR) {
            free(out);
            out = NULL;
            n_cur = 0;
            break;
        }

        n_cur += SLOW5_Z_OUT_CHUNK - strm->avail_out;

    } while (strm->avail_out == 0);

    *n = n_cur;
    inflateReset(strm);

    return out;
}

static void *ptr_depress_gzip_multi(const void *ptr, size_t count, size_t *n) {
    uint8_t *out = NULL;

    size_t n_cur = 0;

    z_stream strm_local;
    gzip_init_inflate(&strm_local);
    z_stream *strm = &strm_local;

    strm->avail_in = count;
    strm->next_in = (Bytef *) ptr;

    do {
        out = (uint8_t *) realloc(out, n_cur + SLOW5_Z_OUT_CHUNK);
        SLOW5_ASSERT(out);

        strm->avail_out = SLOW5_Z_OUT_CHUNK;
        strm->next_out = out + n_cur;

        int ret = inflate(strm, Z_NO_FLUSH);
        if (ret == Z_STREAM_ERROR || ret == Z_DATA_ERROR) {
            SLOW5_WARNING("%s","inflate failed");
            free(out);
            out = NULL;
            n_cur = 0;
            break;
        }

        n_cur += SLOW5_Z_OUT_CHUNK - strm->avail_out;

    } while (strm->avail_out == 0);

    *n = n_cur;

    (void) inflateEnd(strm);

    return out;
}


/* --- Compress / decompress a ptr to some file --- */

size_t slow5_fwrite_compress(struct slow5_press *comp, const void *ptr, size_t size, size_t nmemb, FILE *fp) {
    size_t bytes = -1;

    if (comp != NULL) {
        switch (comp->method) {

            case SLOW5_COMPRESS_NONE:
                bytes = fwrite(ptr, size, nmemb, fp);
                break;

            case SLOW5_COMPRESS_GZIP:
                if (comp->stream != NULL && comp->stream->gzip != NULL) {
                    bytes = fwrite_compress_gzip(comp->stream->gzip, ptr, size, nmemb, fp);
                }
                break;
        }
    }

    return bytes;
}

static size_t fwrite_compress_gzip(struct slow5_gzip_stream *gzip, const void *ptr, size_t size, size_t nmemb, FILE *fp) {

    size_t bytes = 0;
    z_stream *strm = &(gzip->strm_deflate);

    strm->avail_in = size * nmemb;
    strm->next_in = (Bytef *) ptr;

    uLong chunk_sz = SLOW5_Z_OUT_CHUNK;
    uint8_t *buf = (uint8_t *) malloc(sizeof *buf * chunk_sz);
    if (buf == NULL) {
        return -1;
    }

    do {
        strm->avail_out = chunk_sz;
        strm->next_out = buf;

        if (deflate(strm, gzip->flush) == Z_STREAM_ERROR) {
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

    if (gzip->flush == Z_FINISH) {
        gzip->flush = Z_NO_FLUSH;
    }

    return bytes;
}


/* --- Decompress to a ptr from some file --- */

void *slow5_fread_depress(struct slow5_press *comp, size_t count, FILE *fp, size_t *n) {
    void *raw = (void *) malloc(count);
    SLOW5_MALLOC_CHK(raw);

    if (fread(raw, count, 1, fp) != 1) {
        free(raw);
        return NULL;
    }

    void *out = slow5_ptr_depress(comp, raw, count, n);
    free(raw);

    return out;
}
/*
static void *fread_depress_multi(slow5_press_method_t method, size_t count, FILE *fp, size_t *n) {
    void *raw = (void *) malloc(count);
    SLOW5_MALLOC_CHK(raw);

    if (fread(raw, count, 1, fp) != 1) {
        free(raw);
        return NULL;
    }

    void *out = slow5_ptr_depress_multi(method, raw, count, n);
    free(raw);

    return out;
}
*/
void *slow5_pread_depress(struct slow5_press *comp, int fd, size_t count, off_t offset, size_t *n) {
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

void *slow5_pread_depress_multi(slow5_press_method_t method, int fd, size_t count, off_t offset, size_t *n) {
    void *raw = (void *) malloc(count);
    SLOW5_MALLOC_CHK(raw);

    if (pread(fd, raw, count, offset) == -1) {
        SLOW5_WARNING("pread could not read %ld bytes as expected.",(long)count);
        free(raw);
        return NULL;
    }

    void *out = slow5_ptr_depress_multi(method, raw, count, n);
    free(raw);

    return out;
}


/* --- Compress with printf to some file --- */

int slow5_fprintf_compress(struct slow5_press *comp, FILE *fp, const char *format, ...) {
    int ret = -1;

    va_list ap;
    va_start(ap, format);

    ret = vfprintf_compress(comp, fp, format, ap);

    va_end(ap);

    return ret;
}

int slow5_printf_compress(struct slow5_press *comp, const char *format, ...) {
    int ret = -1;

    va_list ap;
    va_start(ap, format);

    ret = vfprintf_compress(comp, stdout, format, ap);

    va_end(ap);

    return ret;
}

static int vfprintf_compress(struct slow5_press *comp, FILE *fp, const char *format, va_list ap) {
    int ret = -1;

    if (comp != NULL) {

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

void slow5_compress_footer_next(struct slow5_press *comp) {

    if (comp != NULL && comp->stream != NULL) {

        switch (comp->method) {

            case SLOW5_COMPRESS_NONE:
                break;

            case SLOW5_COMPRESS_GZIP: {

                struct slow5_gzip_stream *gzip = comp->stream->gzip;

                if (gzip != NULL) {
                    gzip->flush = Z_FINISH;
                }

            } break;
        }
    }
}




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

    int ret = inflateInit2(&strm, GZIP_WBITS);

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

    uLong out_sz = SLOW5_Z_OUT_CHUNK;
    unsigned char *out = (unsigned char *) malloc(sizeof *out * out_sz);

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
