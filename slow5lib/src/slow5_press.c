#define _XOPEN_SOURCE 700
#include <zlib.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <slow5/slow5.h>
#include <slow5/slow5_error.h>
#include <slow5/slow5_press.h>
#include <slow5/slow5_defs.h>
#include "slow5_misc.h"

extern enum slow5_log_level_opt  slow5_log_level;
extern enum slow5_exit_condition_opt  slow5_exit_condition;

static int zlib_init_deflate(z_stream *strm);
static int zlib_init_inflate(z_stream *strm);

static void *ptr_compress_zlib(struct slow5_zlib_stream *zlib, const void *ptr, size_t count, size_t *n);
static void *ptr_depress_zlib(struct slow5_zlib_stream *zlib, const void *ptr, size_t count, size_t *n);
static void *ptr_depress_zlib_solo(const void *ptr, size_t count, size_t *n);
static size_t fwrite_compress_zlib(struct slow5_zlib_stream *zlib, const void *ptr, size_t size, size_t nmemb, FILE *fp);
static int vfprintf_compress(struct slow5_press *comp, FILE *fp, const char *format, va_list ap);

/* --- Init / free slow5_press structure --- */

/*
 * init slow5_press struct
 * returns NULL on error and sets errno
 * errors
 * SLOW5_ERR_MEM
 * SLOW5_ERR_ARG    method is bad
 * SLOW5_ERR_PRESS  (de)compression failure
 */
struct slow5_press *slow5_press_init(slow5_press_method_t method) {

    struct slow5_press *comp = NULL;

    comp = (struct slow5_press *) calloc(1, sizeof *comp);
    if (!comp) {
        SLOW5_MALLOC_ERROR();
        slow5_errno = SLOW5_ERR_MEM;
        return NULL;
    }
    comp->method = method;

    if (method == SLOW5_COMPRESS_NONE) {

    } else if (method == SLOW5_COMPRESS_ZLIB) {
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

    } else {
        SLOW5_ERROR("Invalid (de)compression method '%d'.", method);
        free(comp);
        slow5_errno = SLOW5_ERR_ARG;
        return NULL;
    }

    return comp;
}

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

void slow5_press_free(struct slow5_press *comp) {

    if (comp) {

        switch (comp->method) {

            case SLOW5_COMPRESS_NONE:
                break;

            case SLOW5_COMPRESS_ZLIB: {
                (void) deflateEnd(&(comp->stream->zlib->strm_deflate));
                (void) inflateEnd(&(comp->stream->zlib->strm_inflate));
                free(comp->stream->zlib);
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
        }
    }

    if (n) {
        *n = n_tmp;
    }

    return out;
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

/*
 * ptr cannot be NULL
 */
void *slow5_ptr_depress_solo(slow5_press_method_t method, const void *ptr, size_t count, size_t *n) {
    void *out = NULL;
    size_t n_tmp = 0;

    if (!ptr) {
        SLOW5_ERROR("Argument '%s' cannot be NULL.", SLOW5_TO_STR(ptr))
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
void *slow5_ptr_depress(struct slow5_press *comp, const void *ptr, size_t count, size_t *n) {
    void *out = NULL;

    if (!comp || !ptr) {
        if (!comp) {
            SLOW5_ERROR("Argument '%s' cannot be NULL.", SLOW5_TO_STR(comp))
        }
        if (!ptr) {
            SLOW5_ERROR("Argument '%s' cannot be NULL.", SLOW5_TO_STR(ptr))
        }
        if (n) {
            *n = 0;
        }
        return NULL;
    } else {

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
                    out = ptr_depress_zlib(comp->stream->zlib, ptr, count, n);
                    if (!out) {
                        SLOW5_ERROR("%s", "zlib decompression failed.")
                    }
                }
                break;
        }
    }

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

    if (n) {
        *n = n_cur;
    }
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
        if (ret == Z_STREAM_ERROR || ret == Z_DATA_ERROR) {
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


/* --- Compress / decompress a ptr to some file --- */

size_t slow5_fwrite_compress(struct slow5_press *comp, const void *ptr, size_t size, size_t nmemb, FILE *fp) {
    size_t bytes = -1;

    if (comp) {
        switch (comp->method) {

            case SLOW5_COMPRESS_NONE:
                bytes = fwrite(ptr, size, nmemb, fp);
                break;

            case SLOW5_COMPRESS_ZLIB:
                if (comp->stream && comp->stream->zlib) {
                    bytes = fwrite_compress_zlib(comp->stream->zlib, ptr, size, nmemb, fp);
                }
                break;
        }
    }

    return bytes;
}

static size_t fwrite_compress_zlib(struct slow5_zlib_stream *zlib, const void *ptr, size_t size, size_t nmemb, FILE *fp) {

    size_t bytes = 0;
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


/*
 * decompress count bytes from some file
 * returns pointer to decompressed memory of size *n bytes to be later freed
 * returns NULL on error and *n set to 0
 * DONE
 */
void *slow5_fread_depress(struct slow5_press *comp, size_t count, FILE *fp, size_t *n) {
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

/*
 * pread size bytes from offset then decompress the data according to the method
 */
void *slow5_pread_depress_solo(slow5_press_method_t method, int fd, size_t count, off_t offset, size_t *n) {
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

void slow5_compress_footer_next(struct slow5_press *comp) {

    if (comp && comp->stream) {
        switch (comp->method) {
            case SLOW5_COMPRESS_NONE:
                break;
            case SLOW5_COMPRESS_ZLIB: {
                struct slow5_zlib_stream *zlib = comp->stream->zlib;

                if (zlib) {
                    zlib->flush = Z_FINISH;
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
