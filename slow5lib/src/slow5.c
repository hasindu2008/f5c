/**
 * @file slow5.c
 * @brief SLOW5 implementation
 * @author Sasha Jenner (jenner.sasha@gmail.com), Hasindu Gamaarachchi (hasindu@garvan.org.au), Hiruna Samarakoon
 * @date 27/02/2021
 */

/*
MIT License

Copyright (c) 2020 Hasindu Gamaarachchi
Copyright (c) 2020 Sasha Jenner
Copyright (c) 2020 Hiruna Samarakoon

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

#define _XOPEN_SOURCE 700
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <slow5/slow5.h>
#include "slow5_extra.h"
#include "slow5_idx.h"
#include "slow5_misc.h"
#include "klib/ksort.h"


/* IMPORTANT: The comments in this are NOT the API documentation
The API documentation is available at https://hasindu2008.github.io/slow5tools/
The comments here are for internal use and do not rely on them. Open a GitHub issue for any questions.
*/

KSORT_INIT(str_slow5, ksstr_t, ks_lt_str)


// TODO fail with getline if end of file occurs on a non-empty line
// TODO (void) cast if ignoring return value
// TODO sizeof of macros at compile time rather than strlen
// TODO Make sure all mallocs are checked for success
// TODO check klib mallocs

#define INT16_MAX_LENGTH (6) /* Max length is 6 (−32768) for a int16_t */
#define UINT32_MAX_LENGTH (10) /* Max length is 10 (4294967295) for a uint32_t */

/* TODO is this too much? Or put to a page length */
#define SLOW5_HDR_DATA_BUF_INIT_CAP (1024) /* Initial string buffer capacity for parsing the data header: 2^10 */
#define SLOW5_SIGNAL_BUF_FIXED_CAP (8) /* Fixed string buffer capacity for storing signal: 2^3 since INT16_MAX_LENGTH=6 */
/* TODO is this good? Or put to a page length */
#define SLOW5_HDR_STR_INIT_CAP (1024) /* Initial capacity for converting the header to a string: 2^10 */
#define SLOW5_AUX_META_CAP_INIT (32) /* Initial capacity for the number of auxiliary fields: 2^5 */
#define SLOW5_AUX_ENUM_LABELS_CAP_INIT (32) /* Initial capacity for the number of enum labels: 2^5 */
#define SLOW5_AUX_ARRAY_CAP_INIT (256) /* Initial capacity for parsing auxiliary array: 2^8 */
#define SLOW5_AUX_ARRAY_STR_CAP_INIT (1024) /* Initial capacity for storing auxiliary array string: 2^10 */

#define SLOW5_FSTREAM_BUFF_SIZE (131072)  /* buffer size for freads and fwrites */

static inline void slow5_free(struct slow5_file *s5p);
static int slow5_rec_aux_parse(char *tok, char *read_mem, uint64_t offset, size_t read_size, struct slow5_rec *read, enum slow5_fmt format, struct slow5_aux_meta *aux_meta);
static inline khash_t(slow5_s2a) *slow5_rec_aux_init(void);
static inline void slow5_rec_set_aux_map(khash_t(slow5_s2a) *aux_map, const char *field, const uint8_t *data, size_t len, uint64_t bytes, enum slow5_aux_type type);
static char *get_missing_str(size_t *len);
static int slow5_version_sanity(struct slow5_hdr *hdr);
static struct slow5_version slow5_press_version_bump(struct slow5_version current, slow5_press_method_t method);

static inline slow5_file_t *slow5_open_write(const char *filename);
static inline slow5_file_t *slow5_open_append(const char *filename,  enum slow5_fmt format);

enum slow5_log_level_opt slow5_log_level = SLOW5_LOG_INFO;
enum slow5_exit_condition_opt slow5_exit_condition = SLOW5_EXIT_OFF;

__thread int slow5_errno_intern = 0;

inline int *slow5_errno_location(void) {
    return &slow5_errno_intern;
}

/* Definitions */

// slow5 file


/*
 * core function for opening a SLOW5 file
 * if slow5_fmt is SLOW5_FORMAT_UNKNOWN determines format from pathname extension
 * parses and populates the header
 * initialises decompression buffers
 *
 * errors
 * SLOW5_ERR_ARG    fp is NULL
 * SLOW5_ERR_UNK    format could not be determined from extension
 * SLOW5_ERR_MEM    memory allocation failed
 * slow5_hdr_init errors
 * slow5_press_init errors
 * SLOW5_ERR_IO     ftello, fileno failed
 */
struct slow5_file *slow5_init(FILE *fp, const char *pathname, enum slow5_fmt format) {
    /* pathname allowed to be NULL at this point */
    if (!fp) {
        SLOW5_ERROR("Argument '%s' cannot be NULL.", SLOW5_TO_STR(fp));
        slow5_errno = SLOW5_ERR_ARG;
        return NULL;
    }

    // Attempt to determine format from pathname
    if (format == SLOW5_FORMAT_UNKNOWN &&
            (format = slow5_path_get_fmt(pathname)) == SLOW5_FORMAT_UNKNOWN) {
        SLOW5_ERROR("Unknown slow5 format for file '%s'. Extension must be '%s' or '%s'.",
                pathname, SLOW5_ASCII_EXTENSION, SLOW5_BINARY_EXTENSION);
        slow5_errno = SLOW5_ERR_UNK;
        return NULL;
    }

    char *fread_buff = (char *)calloc(SLOW5_FSTREAM_BUFF_SIZE, sizeof(char));
    if (!fread_buff) {
        SLOW5_MALLOC_ERROR();
        slow5_errno = SLOW5_ERR_MEM;
        return NULL;
    }

    if (setvbuf(fp, fread_buff, _IOFBF, SLOW5_FSTREAM_BUFF_SIZE) != 0){
        SLOW5_WARNING("Could not set a large buffer for file stream of '%s': %s.", pathname, strerror(errno));;
        free(fread_buff);
        fread_buff = NULL;
    }
    else {
        SLOW5_LOG_DEBUG("Buffer for file stream of '%s' was set to %d.", pathname, SLOW5_FSTREAM_BUFF_SIZE);
    }

    // TODO Attempt to determine from magic number

    slow5_press_method_t method;
    struct slow5_hdr *header = slow5_hdr_init(fp, format, &method);
    if (!header) {
        free(fread_buff);
        SLOW5_ERROR("Parsing slow5 header of file '%s' failed.", pathname);
        return NULL;
    }

    struct slow5_file *s5p = (struct slow5_file *) calloc(1, sizeof *s5p);
    if (!s5p) {
        SLOW5_MALLOC_ERROR();
        slow5_hdr_free(header);
        free(fread_buff);
        slow5_errno = SLOW5_ERR_MEM;
        return NULL;
    }

    s5p->fp = fp;
    s5p->format = format;
    s5p->header = header;
    s5p->meta.fread_buffer = fread_buff;

    s5p->compress = slow5_press_init(method);
    if (!s5p->compress) {
        free(fread_buff);
        slow5_hdr_free(header);
        free(s5p);
        return NULL;
    }

    if ((s5p->meta.fd = fileno(fp)) == -1) {
        SLOW5_ERROR("Obtaining file descriptor with fileno() failed: %s.", strerror(errno));
        free(fread_buff);
        slow5_press_free(s5p->compress);
        slow5_hdr_free(header);
        free(s5p);
        slow5_errno = SLOW5_ERR_IO;
        return NULL;
    }

    s5p->meta.pathname = pathname;
    if ((s5p->meta.start_rec_offset = ftello(fp)) == -1) {
        SLOW5_ERROR("Obtaining file offset with ftello() failed: %s.", strerror(errno));
        free(fread_buff);
        slow5_press_free(s5p->compress);
        slow5_hdr_free(header);
        free(s5p);
        slow5_errno = SLOW5_ERR_IO;
        return NULL;
    }

    return s5p;
}

/*
 * initialise an empty SLOW5 file structure
 * if slow5_fmt is SLOW5_FORMAT_UNKNOWN determines format from pathname extension
 * allocate memory for header, SLOW5 file structure and populate file number and offset
 *
 * errors
 * SLOW5_ERR_OTH    big endian machine
 * SLOW5_ERR_ARG    fp is NULL
 * SLOW5_ERR_UNK    format could not be determined from extension
 * SLOW5_ERR_MEM    memory allocation failed
 * SLOW5_ERR_IO     ftello, fileno failed
 */
struct slow5_file *slow5_init_empty(FILE *fp, const char *pathname, enum slow5_fmt format) {

    if (slow5_is_big_endian()) {
        SLOW5_ERROR("%s","Big endian machine detected. slow5lib only support little endian at this time. Please open a github issue stating your machine spec <https://github.com/hasindu2008/slow5lib/issues>.");
        slow5_errno = SLOW5_ERR_OTH;
        return NULL;
    }
    // pathname allowed to be NULL at this point
    if (!fp) {
        SLOW5_ERROR("Argument '%s' cannot be NULL.", SLOW5_TO_STR(fp));
        slow5_errno = SLOW5_ERR_ARG;
        return NULL;
    }

    // Attempt to determine format from pathname
    if (format == SLOW5_FORMAT_UNKNOWN &&
            (format = slow5_path_get_fmt(pathname)) == SLOW5_FORMAT_UNKNOWN) {
        SLOW5_ERROR("Unknown slow5 format for file '%s'. Extension must be '%s' or '%s'.",
                pathname, SLOW5_ASCII_EXTENSION, SLOW5_BINARY_EXTENSION);
        slow5_errno = SLOW5_ERR_UNK; //slow5_path_get_fmt does not set slow5_errno, so set here
        return NULL;
    }

    struct slow5_file *s5p;
    struct slow5_hdr *header = slow5_hdr_init_empty();
    if (!header) {
        SLOW5_ERROR("%s","Initiallising an empty slow5 header failed.");
        return NULL;
    }

    header->version = SLOW5_VERSION_STRUCT;
    s5p = (struct slow5_file *) calloc(1, sizeof *s5p);
    if (!s5p) {
        SLOW5_MALLOC_ERROR();
        slow5_errno = SLOW5_ERR_MEM;
        return NULL;
    }

    s5p->fp = fp;
    s5p->format = format;
    s5p->header = header;

    if ((s5p->meta.fd = fileno(fp)) == -1) {
        SLOW5_ERROR("Obtaining file descriptor with fileno() failed: %s.", strerror(errno));
        slow5_errno = SLOW5_ERR_IO;
        slow5_close(s5p);
        s5p = NULL;
    }
    s5p->meta.pathname = pathname;
    if ((s5p->meta.start_rec_offset = ftello(fp)) == -1) {
        SLOW5_ERROR("Obtaining file offset with ftello() failed: %s.", strerror(errno));
        slow5_errno = SLOW5_ERR_IO;
        slow5_close(s5p);
        s5p = NULL;
    }

    return s5p;
}

/**
 * Open a slow5 file with a specific mode given it's pathname.
 *
 * Attempt to guess the file's slow5 format from the pathname's extension.
 * Return NULL on error.
 *
 * If successful, return a slow5 file structure with the header parsed.
 * slow5_close() should be called when finished with the structure.
 *
 *
 * @param   pathname    relative or absolute path to slow5 file
 * @param   mode        "r" for reading, "w" for writing a new file, "a" for appending to an existing file
 * @return              slow5 file structure
 */
struct slow5_file *slow5_open(const char *pathname, const char *mode) {
    return slow5_open_with(pathname, mode, SLOW5_FORMAT_UNKNOWN);
}

/**
 * Open a slow5 file of a specific format with a mode given it's pathname.
 *
 * Return NULL if pathname or mode is NULL, or if the format specified doesn't match the file.
 * slow5_open_with(pathname, mode, SLOW5_FORMAT_UNKNOWN) is equivalent to slow5_open(pathname, mode).
 *
 * Otherwise, return a slow5 file structure with the header parsed.
 * slow5_close() should be called when finished with the structure.
 *
 * On error, NULL is returned and slow5_errno is set to indicate the error.
 * SLOW5_ERR_OTH    Big endian is not supported.
 * SLOW5_ERR_ARG    The pathname or mode provided was NULL.
 * SLOW5_ERR_IO     The file could not be opened. See errno for details.
 * slow5_init errors
 * slow5_open_write errors
 * slow5_open_append errors
 *
 * @param   pathname    path to slow5 file
 * @param   mode        "r" for reading, "w" for writing a new file, "a" for appending to an existing file
 * @param   format      format of the slow5 file
 * @return              slow5 file structure
 */
struct slow5_file *slow5_open_with(const char *pathname, const char *mode, enum slow5_fmt format) {
    if (slow5_is_big_endian()) {
        SLOW5_ERROR_EXIT("%s", "Big endian machine detected. slow5lib only supports little endian at this time. Please open a github issue stating your machine spec <https://github.com/hasindu2008/slow5lib/issues>.");
        slow5_errno = SLOW5_ERR_OTH;
        return NULL;
    }
    if (!pathname || !mode) {
        if (!pathname) {
            SLOW5_ERROR_EXIT("Argument '%s' cannot be NULL.", SLOW5_TO_STR(pathname));
        }
        if (!mode) {
            SLOW5_ERROR_EXIT("Argument '%s' cannot be NULL.", SLOW5_TO_STR(mode));
        }
        slow5_errno = SLOW5_ERR_ARG;
        return NULL;
    }

    else if (strcmp(mode, "w") == 0){
        struct slow5_file *s5p = slow5_open_write(pathname);
        if(s5p) {
            s5p->meta.mode = mode;
        } else{
            SLOW5_EXIT_IF_ON_ERR();
        }
        return s5p;
    }
    else if (strcmp(mode, "a") == 0){
        struct slow5_file *s5p = slow5_open_append(pathname, format);
        if(s5p) {
            s5p->meta.mode = mode;
        } else {
            SLOW5_EXIT_IF_ON_ERR();
        }
        return s5p;
    }
    else if (strcmp(mode, "r") != 0){
        SLOW5_WARNING("Currently, the only supported modes are 'r', 'w' and 'a'. You entered '%s'.", mode);
    }

    FILE *fp = fopen(pathname, mode);
    if (!fp) {
        SLOW5_ERROR_EXIT("Error opening file '%s': %s.", pathname, strerror(errno));
        slow5_errno = SLOW5_ERR_IO;
        return NULL;
    }

    struct slow5_file *s5p = slow5_init(fp, pathname, format);
    if (!s5p) {
        if (fclose(fp) == EOF) {
            SLOW5_ERROR("Error closing file '%s': %s.", pathname, strerror(errno));
        }
        SLOW5_EXIT_IF_ON_ERR();
    } else {
        s5p->meta.mode = mode;
    }

    return s5p;
}

/* Creates an empty SLOW5 file with the given pathname,
*  initialise SLOW5 file structure with space for a single read group
*  errors
* SLOW5_ERR_IO     The file could not be opened. See errno for details.
* slow5_init_empty errors
* slow5_aux_meta_init_empty errors
* slow5_press_init errors
*/
static inline slow5_file_t *slow5_open_write(const char *filename){

    FILE *fp = fopen(filename, "w");
    if(fp==NULL){
        SLOW5_ERROR("Error opening file '%s': %s.", filename, strerror(errno));
        slow5_errno = SLOW5_ERR_IO;
        return NULL;
    }

    slow5_file_t *s5p = slow5_init_empty(fp, filename, SLOW5_FORMAT_UNKNOWN);
    if(!s5p){
        SLOW5_ERROR("Error initialising an empty SLOW5 file '%s'",filename);
        fclose(fp);
        return NULL;
    }

    slow5_hdr_t *header=s5p->header;
    if (slow5_hdr_add_rg(header) < 0){ //todo error handling down the chain
        SLOW5_ERROR("Error adding read group 0 for %s",filename);
        slow5_close(s5p);
        return NULL;
    }
    header->num_read_groups = 1;

    struct slow5_aux_meta *aux_meta = slow5_aux_meta_init_empty();
    if(!aux_meta){
        SLOW5_ERROR("Error initializing aux meta for %s",filename);
        slow5_close(s5p);
        return NULL;
    }
    header->aux_meta = aux_meta;

    //this structure is only to be used in single threaded writes
    if(s5p->format == SLOW5_FORMAT_BINARY){
        slow5_press_method_t press_out = {SLOW5_COMPRESS_ZLIB, SLOW5_COMPRESS_SVB_ZD};
        s5p->compress = slow5_press_init(press_out);
        if(!s5p->compress){
            SLOW5_ERROR("Could not initialise the slow5 compression method. %s","");
            slow5_close(s5p);
            return NULL;
        }
    }

    return s5p;
}

/* Opens a SLOW5 file with the given pathname for appending,
*  populates the SLOW5 file structure and set file pointer to the end of the file
*  errors
* SLOW5_ERR_IO     The file could not be opened or fseek failed. See errno for details.
* SLOW5_ERR_TRUNC  EOF marker not present
* SLOW5_ERR_UNK    Unknown file format (not .slow5 or .blow5 extension)
* slow5_init errors
* slow5_is_eof errors
*/
static inline slow5_file_t *slow5_open_append(const char *filename, enum slow5_fmt format){

    FILE *fp = fopen(filename, "r+");
    if (!fp) {
        SLOW5_ERROR("Error opening file '%s': %s.", filename, strerror(errno));
        slow5_errno = SLOW5_ERR_IO;
        return NULL;
    }

    struct slow5_file *s5p = slow5_init(fp, filename, format);
    if (!s5p) {
        if (fclose(fp) == EOF) {
            SLOW5_ERROR("Error closing file '%s': %s.", filename, strerror(errno));
        }
        return NULL;
    }

    if(s5p->format==SLOW5_FORMAT_BINARY){
        if(fseek(s5p->fp , 0, SEEK_END) !=0 ){
            SLOW5_ERROR("Fseek to the end of file (SEEK_END) failed '%s': %s.", filename, strerror(errno));
            slow5_errno = SLOW5_ERR_IO;
            slow5_close(s5p);
            return NULL;
        }
        const char eof[] = SLOW5_BINARY_EOF;
        if(slow5_is_eof(s5p->fp, eof, sizeof eof)!=1){
            SLOW5_ERROR("No valid slow5 EOF marker at the end of the SLOW5 file %s.",filename);
            slow5_close(s5p);
            return NULL;
        }
    }

    if(s5p->format==SLOW5_FORMAT_BINARY){
        const char eof[] = SLOW5_BINARY_EOF;
        if(fseek(s5p->fp, - (sizeof *eof) * (sizeof eof) , SEEK_END) != 0){
            SLOW5_ERROR("Fseek to the end of file (SEEK_END-eof_marker_size) failed '%s': %s.", filename, strerror(errno));
            slow5_errno = SLOW5_ERR_IO;
            slow5_close(s5p);
            return NULL;
        }
    } else if (s5p->format==SLOW5_FORMAT_ASCII){
        if(fseek(s5p->fp, 0 , SEEK_END) != 0){
            SLOW5_ERROR("Fseek to the end of file failed '%s': %s.", filename, strerror(errno));
            slow5_errno = SLOW5_ERR_IO;
            slow5_close(s5p);
            return NULL;
        }
    } else {
        SLOW5_ERROR("Unknown slow5 format for file '%s'. Extension must be '%s' or '%s'.",
                filename, SLOW5_ASCII_EXTENSION, SLOW5_BINARY_EXTENSION);
        slow5_errno = SLOW5_ERR_UNK;
        slow5_close(s5p);
        return NULL;
    }

    return s5p;

}

/* Close a slow5 file and free its memory.
 *
 * Return:
 *  0   file successfully closed and memory freed
 *  EOF error occured
 *
 * Errors:
 * SLOW5_ERR_IO     file closing error, see errno for details
 * slow5_idx_write errors
 *
 * @param   s5p     slow5 file
 * @return  see above
 */
int slow5_close(struct slow5_file *s5p) {
    int ret = 0;

    if (!s5p) {
        ret = EOF;
    } else {

        if(s5p->meta.mode && (strcmp(s5p->meta.mode, "w") == 0 || strcmp(s5p->meta.mode, "a") == 0)){
            if(s5p->format == SLOW5_FORMAT_BINARY){
                SLOW5_LOG_DEBUG("Writing EOF marker to file '%s'", s5p->meta.pathname);
                if(slow5_eof_fwrite(s5p->fp) < 0){
                    SLOW5_ERROR_EXIT("%s","Error writing EOF!\n");
                    slow5_errno = SLOW5_ERR_IO;
                    ret = EOF;
                }
            }
        }

        if (fclose(s5p->fp) == EOF) {
            SLOW5_ERROR("Error closing slow5 file '%s': %s.", s5p->meta.pathname, strerror(errno)); //not critical - so do not call exit
            slow5_errno = SLOW5_ERR_IO;
            ret = EOF;
        }

        // TODO slow5_index_close
        if (s5p->index && s5p->index->fp && s5p->index->dirty) { // if the index has been changed, write it back
            if (fseek(s5p->index->fp, 0L, SEEK_SET) != 0) {
                SLOW5_ERROR("Failed to fseek() to start of index file '%s': %s.", s5p->index->pathname, strerror(errno));
                slow5_errno = SLOW5_ERR_IO;
                ret = EOF;
            } else {
                int err = slow5_idx_write(s5p->index, s5p->header->version);
                if (err != 0) {
                    SLOW5_ERROR("Writing index file to '%s' failed.", s5p->index->pathname);
                    slow5_errno = err;
                    ret = EOF;
                }
            }
        }

        slow5_free(s5p);
    }

    return ret;
}

static inline void slow5_free(struct slow5_file *s5p) {
    if (s5p) {
        slow5_press_free(s5p->compress);
        slow5_hdr_free(s5p->header);
        slow5_idx_free(s5p->index);
        free(s5p->meta.fread_buffer);
        free(s5p);
    }
}

int slow5_set_press(slow5_file_t *s5p, enum slow5_press_method rec_press, enum slow5_press_method sig_press){

    if(s5p==NULL){
        SLOW5_ERROR_EXIT("Argument '%s' cannot be NULL.", SLOW5_TO_STR(s5p));
        slow5_errno = SLOW5_ERR_ARG;
        return -1;
    }
    //only works in binary mode and opened in 'w' mode
    if(!(s5p->meta.mode && strcmp(s5p->meta.mode,"w")==0)){
        SLOW5_ERROR_EXIT("%s","File must have been opened for writing.");
        slow5_errno = SLOW5_ERR_ARG;
        return -1;
    }

    if(!(s5p->format == SLOW5_FORMAT_BINARY)){
        SLOW5_ERROR_EXIT("%s","File should be in binary format (blow5).");
        slow5_errno = SLOW5_ERR_ARG;
        return -1;
    }


    //free the existing press if any
    slow5_press_free(s5p->compress);

    //this structure is only to be used in single threaded writes
    if(s5p->format == SLOW5_FORMAT_BINARY){
        slow5_press_method_t press_out = {rec_press,sig_press};
        s5p->compress = slow5_press_init(press_out);
        if(!s5p->compress){
            slow5_errno = SLOW5_ERR_PRESS;
            SLOW5_ERROR_EXIT("Could not initialise the slow5 compression method. %s","");
            return -1;
        }
    }

    return 0;

}

// slow5 header

struct slow5_hdr *slow5_hdr_init_empty(void) {

    struct slow5_hdr *header = (struct slow5_hdr *) calloc(1, sizeof *(header));
    if (!header) {
        SLOW5_MALLOC_ERROR();
        slow5_errno = SLOW5_ERR_MEM;
        return NULL;
    }

    return header;
}


/*
 * parses a slow5 header
 * errors
 * SLOW5_ERR_ARG        fp or method is NULL
 * SLOW5_ERR_MEM        memory allocation failed
 * SLOW5_ERR_HDRPARSE   header parsing error
 * SLOW5_ERR_TRUNC      eof reached prematurely
 * SLOW5_ERR_IO
 * SLOW5_ERR_MAGIC
 * SLOW5_ERR_VERSION
 * slow5_hdr_data_init errors
 * slow5_aux_meta_init errors
 */
struct slow5_hdr *slow5_hdr_init(FILE *fp, enum slow5_fmt format, slow5_press_method_t *method) {

    if (!fp || !method) {
        if (!fp) {
            SLOW5_ERROR("Argument '%s' cannot be NULL.", SLOW5_TO_STR(fp));
        }
        if (!method) {
            SLOW5_ERROR("Argument '%s' cannot be NULL.", SLOW5_TO_STR(method));
        }
        slow5_errno = SLOW5_ERR_ARG;
        return NULL;
    }

    struct slow5_hdr *header = (struct slow5_hdr *) calloc(1, sizeof *(header));
    if (!header) {
        SLOW5_MALLOC_ERROR();
        slow5_errno = SLOW5_ERR_MEM;
        return NULL;
    }

    method->signal_method = SLOW5_COMPRESS_NONE;
    struct slow5_version max_supported = SLOW5_VERSION_ARRAY;

    char *buf = NULL;

    // Parse slow5 header
    if (format == SLOW5_FORMAT_ASCII) {
        method->record_method = SLOW5_COMPRESS_NONE;

        // Buffer for file parsing
        size_t cap = SLOW5_HDR_DATA_BUF_INIT_CAP;
        buf = (char *) malloc(cap * sizeof *buf);
        if (!buf) {
            SLOW5_MALLOC_ERROR();
            free(header);
            return NULL;
        }
        char *bufp;
        ssize_t buf_len;
        int err;

        // 1st line - slow5_version
        if ((buf_len = getline(&buf, &cap, fp)) == -1) {
            SLOW5_ERROR("%s", "Malformed slow5 header. No newline character in whole file.");
            slow5_errno = SLOW5_ERR_HDRPARSE;
            goto err;
        }
        buf[buf_len - 1] = '\0'; // Remove newline for later parsing
        // "#slow5_version"
        bufp = buf;
        char *tok = slow5_strsep(&bufp, SLOW5_SEP_COL);
        if (strcmp(tok, SLOW5_HDR_FILE_VERSION_ID) != 0) {
            SLOW5_ERROR("Malformed slow5 header. Expected '%s', instead found '%s'.",
                    SLOW5_HDR_FILE_VERSION_ID, tok);
            slow5_errno = SLOW5_ERR_HDRPARSE;
            goto err;
        }
        if (!bufp) {
            SLOW5_ERROR("Malformed slow5 header. Missing %s separator after '%s'.", SLOW5_SEP_COL_NAME, tok);
            slow5_errno = SLOW5_ERR_HDRPARSE;
            goto err;
        }

        // Parse file version
        tok = bufp;
        char *toksub;
        toksub = slow5_strsep(&tok, SLOW5_HDR_FILE_VERSION_SEP); // Major version
        header->version.major = slow5_ato_uint8(toksub, &err);
        if (err == -1) {
            SLOW5_ERROR("Malformed slow5 header. Bad major version '%s'.", toksub);
            slow5_errno = SLOW5_ERR_HDRPARSE;
            goto err;
        }

        if (!tok) {
            SLOW5_ERROR("Malformed slow5 header. Missing '%s' separator after major version '%s'.",
                    SLOW5_HDR_FILE_VERSION_SEP, toksub);
            slow5_errno = SLOW5_ERR_HDRPARSE;
            goto err;
        }
        toksub = slow5_strsep(&tok, SLOW5_HDR_FILE_VERSION_SEP); // Minor version
        header->version.minor = slow5_ato_uint8(toksub, &err);
        if (err == -1) {
            SLOW5_ERROR("Malformed slow5 header. Bad minor version '%s'.", toksub);
            slow5_errno = SLOW5_ERR_HDRPARSE;
            goto err;
        }

        if (!tok) {
            SLOW5_ERROR("Malformed slow5 header. Missing '%s' separator after minor version '%s'.",
                    SLOW5_HDR_FILE_VERSION_SEP, toksub);
            slow5_errno = SLOW5_ERR_HDRPARSE;
            goto err;
        }
        toksub = slow5_strsep(&tok, SLOW5_HDR_FILE_VERSION_SEP); // Patch version
        header->version.patch = slow5_ato_uint8(toksub, &err);
        if (err == -1) {
            SLOW5_ERROR("Malformed slow5 header. Bad minor version '%s'.", toksub);
            slow5_errno = SLOW5_ERR_HDRPARSE;
            goto err;
        }

        if (tok) { // Still more data
            SLOW5_ERROR("Malformed slow5 header. More data exists in '%s' after patch version.", tok);
            slow5_errno = SLOW5_ERR_HDRPARSE;
            goto err;
        }

        if (slow5_is_version_compatible(header->version, max_supported) == 0) {
            SLOW5_ERROR("File version '" SLOW5_VERSION_STRING_FORMAT "' is higher than the max slow5 version '" SLOW5_VERSION_STRING "' supported by this slow5lib! Please use a newer version of slow5lib.",
                    header->version.major, header->version.minor, header->version.patch);
            slow5_errno = SLOW5_ERR_VERSION;
            goto err;
        }

        // 2nd line - num_read_groups
        if ((buf_len = getline(&buf, &cap, fp)) == -1) {
            SLOW5_ERROR("%s", "Malformed slow5 header. No newline character after slow5 version.");
            slow5_errno = SLOW5_ERR_HDRPARSE;
            goto err;
        }
        buf[buf_len - 1] = '\0'; // Remove newline for later parsing
        // "#num_read_groups"
        bufp = buf;
        tok = slow5_strsep(&bufp, SLOW5_SEP_COL);
        if (strcmp(tok, SLOW5_HDR_NUM_GROUPS_ID) != 0) {
            SLOW5_ERROR("Malformed slow5 header. Expected '%s', instead found '%s'.", SLOW5_HDR_NUM_GROUPS_ID, tok);
            slow5_errno = SLOW5_ERR_HDRPARSE;
            goto err;
        }
        if (!bufp) {
            SLOW5_ERROR("Malformed slow5 header. Missing %s separator after '%s'.", SLOW5_SEP_COL_NAME, tok);
            slow5_errno = SLOW5_ERR_HDRPARSE;
            goto err;
        }
        // Parse num read groups
        tok = slow5_strsep(&bufp, SLOW5_SEP_COL);
        header->num_read_groups = slow5_ato_uint32(tok, &err);
        if (err == -1 || header->num_read_groups == 0) {
            SLOW5_ERROR("Malformed slow5 header. Invalid number of read groups '%s'.", tok);
            slow5_errno = SLOW5_ERR_HDRPARSE;
            goto err;
        }

        // Parse the data header
        if (slow5_hdr_data_init(fp, &buf, &cap, header, NULL) != 0) {
            goto err;
        }

        // Parse column datatypes and names
        header->aux_meta = slow5_aux_meta_init(fp, &buf, &cap, NULL, &err);
        if (err == -1) {
            slow5_hdr_data_free(header);
            goto err;
        }

    } else if (format == SLOW5_FORMAT_BINARY) {
        const char magic[] = SLOW5_BINARY_MAGIC_NUMBER;

        char buf_magic[sizeof magic + 1];
        uint32_t header_size;

        uint8_t record_method = 0;
        uint8_t signal_method = 0;

        if (fread(buf_magic, sizeof *magic, sizeof magic, fp) != sizeof magic) {
            SLOW5_ERROR("Malformed blow5 header. Failed to read the magic number.%s", feof(fp) ? " EOF reached." : "");
            goto err_fread;
        } else if (memcmp(magic, buf_magic, sizeof *magic * sizeof magic) != 0) {
            SLOW5_ERROR("%s", "Malformed blow5 header. Invalid magic number.");
            free(header);
            slow5_errno = SLOW5_ERR_MAGIC;
            return NULL;
        } else if (fread(&header->version.major, sizeof header->version.major, 1, fp) != 1) {
            SLOW5_ERROR("Malformed blow5 header. Failed to read the major version.%s", feof(fp) ? " EOF reached." : "");
            goto err_fread;
        } else if (fread(&header->version.minor, sizeof header->version.minor, 1, fp) != 1) {
            SLOW5_ERROR("Malformed blow5 header. Failed to read the minor version.%s", feof(fp) ? " EOF reached." : "");
            goto err_fread;
        } else if (fread(&header->version.patch, sizeof header->version.patch, 1, fp) != 1) {
            SLOW5_ERROR("Malformed blow5 header. Failed to read the patch version.%s", feof(fp) ? " EOF reached." : "");
            goto err_fread;
        } else if (fread(&record_method, sizeof record_method, 1, fp) != 1) {
            SLOW5_ERROR("Malformed blow5 header. Failed to read the record compression method.%s", feof(fp) ? " EOF reached." : "");
            goto err_fread;
        } else if (fread(&header->num_read_groups, sizeof header->num_read_groups, 1, fp) != 1) {
            SLOW5_ERROR("Malformed blow5 header. Failed to read the number of read groups.%s", feof(fp) ? " EOF reached." : "");
            goto err_fread;
        } else if (slow5_signal_press_version_cmp(header->version) >= 0 && fread(&signal_method, sizeof signal_method, 1, fp) != 1) {
            SLOW5_ERROR("Malformed blow5 header. Failed to read the signal compression method.%s", feof(fp) ? " EOF reached." : "");
            goto err_fread;
        } else if (fseek(fp, SLOW5_BINARY_HDR_SIZE_OFFSET, SEEK_SET) == -1) {
            SLOW5_ERROR("Failed to fseek() to offset %ld: %s.", SLOW5_BINARY_HDR_SIZE_OFFSET, strerror(errno));
            free(header);
            slow5_errno = SLOW5_ERR_IO;
            return NULL;
        } else if (fread(&header_size, sizeof header_size, 1, fp) != 1) {
            SLOW5_ERROR("Malformed blow5 header. Failed to read the ascii header size.%s", feof(fp) ? " EOF reached." : "");
            goto err_fread;
        }

        if (slow5_is_version_compatible(header->version, max_supported) == 0) {
            SLOW5_ERROR("File version '" SLOW5_VERSION_STRING_FORMAT "' is higher than the max slow5 version '" SLOW5_VERSION_STRING "' supported by this slow5lib! Please use a newer version of slow5lib.",
                    header->version.major, header->version.minor, header->version.patch);
            free(header);
            slow5_errno = SLOW5_ERR_VERSION;
            return NULL;
        }

        method->record_method = slow5_decode_record_press(record_method);
        method->signal_method = slow5_decode_signal_press(signal_method);

        size_t cap = SLOW5_HDR_DATA_BUF_INIT_CAP;
        buf = (char *) malloc(cap * sizeof *buf);
        if (!buf) {
            SLOW5_MALLOC_ERROR();
            free(header);
            slow5_errno = SLOW5_ERR_MEM;
            return NULL;
        }

        // Header data
        uint32_t header_act_size;

        if (slow5_hdr_data_init(fp, &buf, &cap, header, &header_act_size) != 0) {
            goto err;
        }

        int err;
        header->aux_meta = slow5_aux_meta_init(fp, &buf, &cap, &header_act_size, &err);
        if (err == -1) {
            slow5_hdr_data_free(header);
            goto err;
        }

        if (header_act_size != header_size) {
            SLOW5_ERROR("Expected a slow5 header of size '%" PRIu32 "' bytes, but instead '%" PRIu32 "' bytes were read.",
                    header_size, header_act_size);
            free(buf);
            slow5_hdr_free(header);
            slow5_errno = SLOW5_ERR_HDRPARSE;
            return NULL;
        }
    }

    //TODO some sanity checks - could be made errors later - hasindu
    if(slow5_version_sanity(header) != 0){
        SLOW5_WARNING("%s","Version sanity check of the SLOW5 file failed, which means that it does not conform to specification");
    }

    free(buf);
    return header;

    err:
        free(buf);
        free(header);
        return NULL;

    err_fread:
        free(header);
        slow5_errno = SLOW5_ERR_TRUNC ? feof(fp) : SLOW5_ERR_IO;
        return NULL;
}


//get the list of hdr data keys in sorted order (only the returned pointer must be freed, not the ones inside - subject to change)
//currently used in pyslow5
//len is the number of elements in the return char** list
//returns null if there are no data header attributes
const char **slow5_get_hdr_keys(const slow5_hdr_t *header, uint64_t *len) {
    if (len) {
        *len = header->data.num_attrs;
    }
    if (header->data.num_attrs == 0) {
        return NULL;
    }
    const char **data_attrs = (const char **) malloc(header->data.num_attrs * sizeof *data_attrs);
    SLOW5_MALLOC_CHK(data_attrs);
    uint32_t i = 0;
    for (khint_t j = kh_begin(header->data.attrs); j != kh_end(header->data.attrs); ++ j) {
        if (kh_exist(header->data.attrs, j)) {
            data_attrs[i] = kh_key(header->data.attrs, j);
            ++ i;
        }
    }

    // Sort header data attributes alphabetically
    ks_mergesort(str_slow5, header->data.num_attrs, data_attrs, 0);

    return data_attrs;
}

/**
 * Get the header as a string in the specified format.
 *
 * Returns NULL if s5p is NULL
 * or format is SLOW5_FORMAT_UNKNOWN
 * or an internal error occurs.
 *
 * @param   header  slow5 header
 * @param   format  slow5 format to write the entry in
 * @param   comp    compression method
 * @param   n       number of bytes written to the returned buffer
 * @return  malloced memory storing the slow5 header representation,
 *          to use free() on afterwards
 */
// TODO don't allow comp of SLOW5_COMPRESS_ZLIB for SLOW5_FORMAT_ASCII

//  flattened header returned as a void * (incase of BLOW5 magic number is also included)

void *slow5_hdr_to_mem(struct slow5_hdr *header, enum slow5_fmt format, slow5_press_method_t comp, size_t *n) {
    char *mem = NULL;

    if (header == NULL || format == SLOW5_FORMAT_UNKNOWN) {
        return mem;
    }

    if(slow5_version_sanity(header)!=0){ //TODO is this proper? -hasindu
        SLOW5_ERROR("%s","Version sanity check of the SLOW5 file failed, which means that it does not conform to specification");
        return mem;
    }


    size_t len = 0;
    size_t cap = SLOW5_HDR_STR_INIT_CAP;
    mem = (char *) malloc(cap * sizeof *mem);
    SLOW5_MALLOC_CHK(mem);
    uint32_t header_size;

    if (format == SLOW5_FORMAT_ASCII) {

        struct slow5_version *version = &header->version;

        // Relies on SLOW5_HDR_DATA_BUF_INIT_CAP being bigger than
        // strlen(ASCII_SLOW5_HDR) + UINT32_MAX_LENGTH + strlen("\0")
        int len_ret = sprintf(mem, SLOW5_ASCII_HDR_FORMAT,
                              version->major,
                              version->minor,
                              version->patch,
                              header->num_read_groups);
        if (len_ret <= 0) {
            free(mem);
            return NULL;
        }
        len = len_ret;

    } else if (format == SLOW5_FORMAT_BINARY) {

        struct slow5_version version_tmp = slow5_press_version_bump(header->version,comp);
        struct slow5_version *version = &version_tmp;

        uint8_t record_comp = slow5_encode_record_press(comp.record_method);
        uint8_t signal_comp = slow5_encode_signal_press(comp.signal_method);

        // Relies on SLOW5_HDR_DATA_BUF_INIT_CAP
        // being at least 68 + 1 (for '\0') bytes
        const char magic[] = SLOW5_BINARY_MAGIC_NUMBER;
        memcpy(mem, magic, sizeof magic * sizeof *magic);
        len += sizeof magic * sizeof *magic;
        memcpy(mem + len, &version->major, sizeof version->major);
        len += sizeof version->major;
        memcpy(mem + len, &version->minor, sizeof version->minor);
        len += sizeof version->minor;
        memcpy(mem + len, &version->patch, sizeof version->patch);
        len += sizeof version->patch;
        memcpy(mem + len, &record_comp, sizeof record_comp);
        len += sizeof record_comp;
        memcpy(mem + len, &header->num_read_groups, sizeof header->num_read_groups);
        len += sizeof header->num_read_groups;
        memcpy(mem + len, &signal_comp, sizeof signal_comp);
        len += sizeof signal_comp;

        memset(mem + len, '\0', SLOW5_BINARY_HDR_SIZE_OFFSET - len);
        len = SLOW5_BINARY_HDR_SIZE_OFFSET;

        // Skip header size for later
        len += sizeof header_size;
    }

    size_t len_to_cp;
    // Get unsorted list of header data attributes.
    if (header->data.num_attrs != 0) {
        const char **data_attrs = slow5_get_hdr_keys(header,NULL);

        // Write header data attributes to string
        for (uint32_t i = 0; i < header->data.num_attrs; ++ i) {
            const char *attr = data_attrs[i];

            // Realloc if necessary
            if (len + 1 + strlen(attr) >= cap) { // + 1 for SLOW5_HDR_DATA_PREFIX_CHAR
                cap *= 2;
                mem = (char *) realloc(mem, cap * sizeof *mem);
                SLOW5_MALLOC_CHK(mem);
            }

            mem[len] = SLOW5_HDR_DATA_PREFIX_CHAR;
            ++ len;
            memcpy(mem + len, attr, strlen(attr));
            len += strlen(attr);

            for (uint64_t j = 0; j < (uint64_t) header->num_read_groups; ++ j) {
                const khash_t(slow5_s2s) *hdr_data = header->data.maps.a[j];
                khint_t pos = kh_get(slow5_s2s, hdr_data, attr);

                if (pos != kh_end(hdr_data)) {
                    const char *value = kh_value(hdr_data, pos);

                    // Realloc if necessary
                    if (len + 1 >= cap) { // +1 for SLOW5_SEP_COL_CHAR
                        cap *= 2;
                        mem = (char *) realloc(mem, cap * sizeof *mem);
                        SLOW5_MALLOC_CHK(mem);
                    }

                    mem[len] = SLOW5_SEP_COL_CHAR;
                    ++ len;

                    if (value != NULL) {
                        len_to_cp = strlen(value);

                        // special case for SLOW5_ASCII_MISSING_CHAR
                        bool is_empty = (len_to_cp == 0);
                        if (is_empty) {
                            len_to_cp ++;
                        }

                        // Realloc if necessary
                        if (len + len_to_cp >= cap) {
                            cap *= 2;
                            mem = (char *) realloc(mem, cap * sizeof *mem);
                            SLOW5_MALLOC_CHK(mem);
                        }

                        if (is_empty) { // special case for SLOW5_ASCII_MISSING_CHAR
                            mem[len] = SLOW5_ASCII_MISSING_CHAR;
                        } else {
                            memcpy(mem + len, value, len_to_cp);
                        }
                        len += len_to_cp;
                    } else {
                        // Realloc if necessary
                        if (len + 1 >= cap) { // +1 for SLOW5_ASCII_MISSING_CHAR
                            cap *= 2;
                            mem = (char *) realloc(mem, cap * sizeof *mem);
                            SLOW5_MALLOC_CHK(mem);
                        }
                        mem[len] = SLOW5_ASCII_MISSING_CHAR;
                        ++ len;
                    }
                } else {
                    // Realloc if necessary
                    if (len + 2 >= cap) { // +2 for SLOW5_ASCII_MISSING_CHAR and SLOW5_SEP_COL_CHAR
                        cap *= 2;
                        mem = (char *) realloc(mem, cap * sizeof *mem);
                        SLOW5_MALLOC_CHK(mem);
                    }
                    mem[len] = SLOW5_SEP_COL_CHAR;
                    mem[len + 1] = SLOW5_ASCII_MISSING_CHAR;
                    len += 2;
                }
            }

            // Realloc if necessary
            if (len + 1 >= cap) { // +1 for '\n'
                cap *= 2;
                mem = (char *) realloc(mem, cap * sizeof *mem);
                SLOW5_MALLOC_CHK(mem);
            }

            mem[len] = '\n';
            ++ len;
        }

        free(data_attrs);
    }

    //  data type header
    char *str_to_cp = slow5_hdr_types_to_str(header->aux_meta, &len_to_cp);
    // Realloc if necessary
    if (len + len_to_cp >= cap) {
        cap *= 2;
        mem = (char *) realloc(mem, cap * sizeof *mem);
        SLOW5_MALLOC_CHK(mem);
    }
    memcpy(mem + len, str_to_cp, len_to_cp);
    len += len_to_cp;
    free(str_to_cp);

    // Column header
    str_to_cp = slow5_hdr_attrs_to_str(header->aux_meta, &len_to_cp);
    // Realloc if necessary
    if (len + len_to_cp >= cap) {
        cap *= 2;
        mem = (char *) realloc(mem, cap * sizeof *mem);
        SLOW5_MALLOC_CHK(mem);
    }
    memcpy(mem + len, str_to_cp, len_to_cp);
    len += len_to_cp;
    free(str_to_cp);

    if (format == SLOW5_FORMAT_ASCII) {
        // Realloc if necessary
        if (len + 1 >= cap) { // +1 for '\0'
            cap *= 2;
            mem = (char *) realloc(mem, cap * sizeof *mem);
            SLOW5_MALLOC_CHK(mem);
        }

        mem[len] = '\0';
    } else if (format == SLOW5_FORMAT_BINARY) { //write the header size in bytes (which was skipped previously)
        header_size = len - (SLOW5_BINARY_HDR_SIZE_OFFSET + sizeof header_size);
        memcpy(mem + SLOW5_BINARY_HDR_SIZE_OFFSET, &header_size, sizeof header_size);
    }

    if (n != NULL) {
        *n = len;
    }

    return (void *) mem;
}

char *slow5_hdr_types_to_str(struct slow5_aux_meta *aux_meta, size_t *len) {
    char *types = NULL;
    size_t types_len = 0;

    if (aux_meta != NULL) {
        size_t types_cap = SLOW5_HDR_STR_INIT_CAP;
        types = (char *) malloc(types_cap);
        SLOW5_MALLOC_CHK(types);

        // Assumption that SLOW5_HDR_STR_INIT_CAP > strlen(SLOW5_ASCII_TYPE_HDR_MIN)
        const char *str_to_cp = SLOW5_ASCII_TYPE_HDR_MIN;
        size_t len_to_cp = strlen(str_to_cp);
        memcpy(types, str_to_cp, len_to_cp);
        types_len += len_to_cp;

        for (uint64_t i = 0; i < aux_meta->num; ++ i) {
            str_to_cp = SLOW5_AUX_TYPE_META[aux_meta->types[i]].type_str;
            len_to_cp = strlen(str_to_cp);

            while (types_len + len_to_cp + 1 >= types_cap) { // +1 for SLOW5_SEP_COL_CHAR
                types_cap *= 2;
                types = (char *) realloc(types, types_cap);
                SLOW5_MALLOC_CHK(types);
            }

            types[types_len ++] = SLOW5_SEP_COL_CHAR;
            memcpy(types + types_len, str_to_cp, len_to_cp);
            types_len += len_to_cp;

            /* add enum labels */
            if (aux_meta->types[i] == SLOW5_ENUM ||
                    aux_meta->types[i] == SLOW5_ENUM_ARRAY) {

                if (!aux_meta->enum_num_labels || aux_meta->enum_num_labels[i] == 0) {
                    /* something went wrong */
                    return NULL;
                }

                if (types_len + 1 >= types_cap) { /* +1 for SLOW5_HDR_ENUM_LABELS_BEGIN */
                    types_cap *= 2;
                    types = (char *) realloc(types, types_cap);
                    SLOW5_MALLOC_CHK(types);
                }
                types[types_len ++] = SLOW5_HDR_ENUM_LABELS_BEGIN;

                /* copy over first label without comma */
                str_to_cp = aux_meta->enum_labels[i][0];
                len_to_cp = strlen(str_to_cp);

                while (types_len + len_to_cp >= types_cap) {
                    types_cap *= 2;
                    types = (char *) realloc(types, types_cap);
                    SLOW5_MALLOC_CHK(types);
                }

                memcpy(types + types_len, str_to_cp, len_to_cp);
                types_len += len_to_cp;

                /* copy over rest of the labels with preceding commas */
                for (uint16_t j = 1; j < aux_meta->enum_num_labels[i]; ++ j) {
                    str_to_cp = aux_meta->enum_labels[i][j];
                    len_to_cp = strlen(str_to_cp);

                    while (types_len + len_to_cp + 1 >= types_cap) { // +1 for SLOW5_SEP_ARRAY_CHAR
                        types_cap *= 2;
                        types = (char *) realloc(types, types_cap);
                        SLOW5_MALLOC_CHK(types);
                    }

                    types[types_len ++] = SLOW5_SEP_ARRAY_CHAR;
                    memcpy(types + types_len, str_to_cp, len_to_cp);
                    types_len += len_to_cp;
                }

                if (types_len + 1 >= types_cap) { /* +1 for SLOW5_HDR_ENUM_LABELS_END */
                    types_cap *= 2;
                    types = (char *) realloc(types, types_cap);
                    SLOW5_MALLOC_CHK(types);
                }
                types[types_len ++] = SLOW5_HDR_ENUM_LABELS_END;
            }
        }

        if (types_len + 2 >= types_cap) { // +2 for '\n' and '\0'
            types_cap *= 2;
            types = (char *) realloc(types, types_cap);
            SLOW5_MALLOC_CHK(types);
        }
        types[types_len ++] = '\n';
        types[types_len] = '\0';

    } else {
        types = strdup(SLOW5_ASCII_TYPE_HDR_MIN "\n");
        types_len = strlen(types);
    }

    *len = types_len;

    return types;
}

char *slow5_hdr_attrs_to_str(struct slow5_aux_meta *aux_meta, size_t *len) {
    char *attrs = NULL;
    size_t attrs_len = 0;

    if (aux_meta != NULL) {
        size_t attrs_cap = SLOW5_HDR_STR_INIT_CAP;
        attrs = (char *) malloc(attrs_cap);
        SLOW5_MALLOC_CHK(attrs);

        // Assumption that SLOW5_HDR_STR_INIT_CAP > strlen(SLOW5_ASCII_TYPE_HDR_MIN)
        const char *str_to_cp = SLOW5_ASCII_COLUMN_HDR_MIN;
        size_t len_to_cp = strlen(str_to_cp);
        memcpy(attrs, str_to_cp, len_to_cp);
        attrs_len += len_to_cp;

        for (uint64_t i = 0; i < aux_meta->num; ++ i) {
            str_to_cp = aux_meta->attrs[i];
            len_to_cp = strlen(str_to_cp);

            while (attrs_len + len_to_cp + 1 >= attrs_cap) { // +1 for SLOW5_SEP_COL_CHAR
                attrs_cap *= 2;
                attrs = (char *) realloc(attrs, attrs_cap);
                SLOW5_MALLOC_CHK(attrs);
            }

            attrs[attrs_len ++] = SLOW5_SEP_COL_CHAR;
            memcpy(attrs + attrs_len, str_to_cp, len_to_cp);
            attrs_len += len_to_cp;
        }

        if (attrs_len + 2 >= attrs_cap) { // +2 for '\n' and '\0'
            attrs_cap *= 2;
            attrs = (char *) realloc(attrs, attrs_cap);
            SLOW5_MALLOC_CHK(attrs);
        }
        attrs[attrs_len ++] = '\n';
        attrs[attrs_len] = '\0';

    } else {
        attrs = strdup(SLOW5_ASCII_COLUMN_HDR_MIN "\n");
        attrs_len = strlen(attrs);
    }

    *len = attrs_len;

    return attrs;
}


/**
 * Print the header in the specified format to a file pointer.
 *
 * On success, the number of bytes written is returned.
 * On error, -1 is returned.
 *
 * @param   fp              output file pointer
 * @param   header          slow5 header
 * @param   format          slow5 format to write the entry in
 * @param   comp            compression method
 * @return  number of bytes written, -1 on error
 */
int slow5_hdr_fwrite(FILE *fp, struct slow5_hdr *header, enum slow5_fmt format, slow5_press_method_t comp) {
    int ret;
    void *hdr;
    size_t hdr_size;

    if (fp == NULL || header == NULL || (hdr = slow5_hdr_to_mem(header, format, comp, &hdr_size)) == NULL) {
        return -1;
    }

    size_t n = fwrite(hdr, hdr_size, 1, fp);
    if (n != 1) {
        ret = -1;
    } else {
        ret = hdr_size; // TODO is this okay ie. size_t -> int
    }

    free(hdr);
    return ret;
}


int slow5_hdr_write(slow5_file_t *s5p){

    slow5_press_method_t method={SLOW5_COMPRESS_NONE,SLOW5_COMPRESS_NONE};
    if (s5p->format == SLOW5_FORMAT_BINARY){
        method.record_method = s5p->compress->record_press->method;
        method.signal_method = s5p->compress->signal_press->method;
    }
    int ret = slow5_hdr_fwrite(s5p->fp, s5p->header, s5p->format, method);
    return ret;
}


/**
 * Get a header data map.
 *
 * Returns NULL if an input parameter is NULL.
 *
 * @param   read_group     read group number
 * @param   header     pointer to the header
 * @return  the header data map for that read_group, or NULL on error
 */
khash_t(slow5_s2s) *slow5_hdr_get_data(uint32_t read_group, const struct slow5_hdr *header) {
    if (header == NULL || read_group >= header->num_read_groups) {
        return NULL;
    }

    return header->data.maps.a[read_group];
}


/**
 * Get a header data attribute.
 *
 * Returns NULL if the attribute name doesn't exist
 * or an input parameter is NULL.
 *
 * @param   attr    attribute name
 * @param   read_group     read group number
 * @param   header     pointer to the header
 * @return  the attribute's value, or NULL on error
 */
char *slow5_hdr_get(const char *attr, uint32_t read_group, const struct slow5_hdr *header) {
    char *value;

    if (attr == NULL || header == NULL || read_group >= header->num_read_groups) {
        return NULL;
    }

    khash_t(slow5_s2s) *hdr_data = header->data.maps.a[read_group];

    khint_t pos = kh_get(slow5_s2s, hdr_data, attr);
    if (pos == kh_end(hdr_data)) {
        return NULL;
    } else {
        value = kh_value(hdr_data, pos);
    }

    return value;
}

char **slow5_get_aux_names(const slow5_hdr_t *header, uint64_t *len) {

    uint64_t my_len = header->aux_meta ? header->aux_meta->num : 0;

    if (len) {
        *len = my_len;
    }

    if (my_len == 0) {
        return NULL;
    } else {
        return (header->aux_meta->attrs);
    }
}


enum slow5_aux_type *slow5_get_aux_types(const slow5_hdr_t *header, uint64_t *len) {
    uint64_t my_len = (header->aux_meta == NULL) ? 0 : header->aux_meta->num;

    if (len) {
        *len = my_len;
    }

    if (my_len == 0) {
        return NULL;
    } else {
        return (header->aux_meta->types);
    }
}


/**
 * get the enum labels for a specific auxiliary field and set the number of labels in *n
 * return NULL on error and slow5_errno set to
 * SLOW5_ERR_ARG    if header, field NULL, n can be NULL
 * SLOW5_ERR_NOAUX  if auxiliary header is NULL
 * SLOW5_ERR_TYPE   if the enum labels or num_labels array is NULL, or the field type is not an enum type
 * SLOW5_ERR_NOFLD  if the auxiliary field was not found
 * SLOW5_ERR_MEM    memory allocation error
 */
char **slow5_get_aux_enum_labels(const slow5_hdr_t *header, const char *field, uint8_t *n) {
    if (!header || !field) {
        if (!header) {
            SLOW5_ERROR_EXIT("Argument '%s' cannot be NULL.", SLOW5_TO_STR(header));
        }
        if (!field) {
            SLOW5_ERROR_EXIT("Argument '%s' cannot be NULL.", SLOW5_TO_STR(field));
        }
        slow5_errno = SLOW5_ERR_ARG;
        return NULL;
    }

    if (!header->aux_meta) {
        SLOW5_ERROR_EXIT("%s", "Missing auxiliary hash map.");
        slow5_errno = SLOW5_ERR_NOAUX;
        return NULL;
    }

    if (!header->aux_meta->enum_labels || !header->aux_meta->enum_num_labels) {
        SLOW5_ERROR_EXIT("%s", "No enum auxiliary type exists.");
        slow5_errno = SLOW5_ERR_TYPE;
        return NULL;
    }

    uint32_t i;
    khint_t pos = kh_get(slow5_s2ui32, header->aux_meta->attr_to_pos, field);
    if (pos == kh_end(header->aux_meta->attr_to_pos)) {
        SLOW5_ERROR_EXIT("Field '%s' not found.", field);
        slow5_errno = SLOW5_ERR_NOFLD;
        return NULL;
    } else {
        i = kh_val(header->aux_meta->attr_to_pos, pos);
    }

    enum slow5_aux_type type = header->aux_meta->types[i];
    if (type != SLOW5_ENUM && type != SLOW5_ENUM_ARRAY) {
        SLOW5_ERROR_EXIT("Field '%s' does not have type enum.", field);
        slow5_errno = SLOW5_ERR_TYPE;
        return NULL;
    }

    if (n) {
        *n = header->aux_meta->enum_num_labels[i];
    }

    return header->aux_meta->enum_labels[i];
}

/**
 * Add a new header data attribute.
 *
 * Returns -1 if an input parameter is NULL.
 * Returns -2 if the attribute already exists.
 * Returns -3 if internal error.
 * Returns 0 other.
 *
 * @param   attr        attribute name
 * @param   header      pointer to the header
 * @return  0 on success, <0 on error as described above
 */
int slow5_hdr_add_attr(const char *attr, struct slow5_hdr *header) {
    if (attr == NULL || header == NULL) {
        return -1;
    }

    if (header->data.attrs == NULL) {
        header->data.attrs = kh_init(slow5_s);
    }

    // See if attr already there
    if (kh_get(slow5_s, header->data.attrs, attr) == kh_end(header->data.attrs)) {
        // Add attr
        int ret;
        char *attr_cp = strdup(attr);
        kh_put(slow5_s, header->data.attrs, attr_cp, &ret);
        if (ret == -1) {
            free(attr_cp);
            return -3;
        }
        ++ header->data.num_attrs;
    } else {
        return -2;
    }

    return 0;
}


/*
 * Add a new header data attribute.
 *
 * slow5_hdr_add is the wrapper to slow5_hdr_add_attr which exposed to public
 */
int slow5_hdr_add(const char *attr, slow5_hdr_t *header){
    return slow5_hdr_add_attr(attr, header);
}

/**
 * Add a new header read group.
 *
 * All values are set to NULL for the new read group.
 *
 * Returns -1 if an input parameter is NULL.
 * Returns the new read group number otherwise.
 *
 * @param   header  slow5 header
 * @return  < 0 on error as described above
 */
int64_t slow5_hdr_add_rg(struct slow5_hdr *header) {
    int64_t rg_num = -1;

    if (header != NULL) {
        rg_num = header->num_read_groups ++;
        kv_push(khash_t(slow5_s2s) *, header->data.maps, kh_init(slow5_s2s));
    }

    return rg_num;
}

/**
 * Add a new header read group with its data.
 *
 * Returns -1 if an input parameter is NULL.
 * Returns the new read group number otherwise.
 *
 * @param   header      slow5 header
 * @param   new_data    khash map of the new read group
 * @return  < 0 on error as described above
 */
int64_t slow5_hdr_add_rg_data(struct slow5_hdr *header, khash_t(slow5_s2s) *new_data) {
    if (header == NULL || new_data == NULL) {
        return -1;
    }

    int64_t rg_num = slow5_hdr_add_rg(header);

    for (khint_t i = kh_begin(new_data); i != kh_end(new_data); ++ i) {
        if (kh_exist(new_data, i)) {
            const char *attr = kh_key(new_data, i);
            char *value = kh_value(new_data, i);

            if (slow5_hdr_add_attr(attr, header) == -3) {
                SLOW5_ERROR("%s", "Internal klib error.");
                /* TODO remove rg */
                return -1;
            }
            (void) slow5_hdr_set(attr, value, rg_num, header);
        }
    }

    return rg_num;
}

/**
 * Set a header data attribute for a particular read_group.
 *
 * Returns -1 if the attribute name doesn't exist
 * or the read group is out of range
 * or an input parameter is NULL.
 * Returns 0 other.
 *
 * @param   attr        attribute name
 * @param   value       new attribute value
 * @param   read_group  the read group
 * @param   s5p         slow5 file
 * @return  0 on success, -1 on error
 */
int slow5_hdr_set(const char *attr, const char *value, uint32_t read_group, struct slow5_hdr *header) {

    if (attr == NULL || value == NULL || header == NULL || read_group >= header->num_read_groups) {
        return -1;
    }

    khint_t pos = kh_get(slow5_s, header->data.attrs, attr);
    if (pos != kh_end(header->data.attrs)) {
        const char *attr_lib = kh_key(header->data.attrs, pos);
        khash_t(slow5_s2s) *map = header->data.maps.a[read_group];

        khint_t k = kh_get(slow5_s2s, map, attr);
        if (k != kh_end(map)) {
            free(kh_value(map, k));
            kh_value(map, k) = strdup(value);
        } else {
            int ret;
            k = kh_put(slow5_s2s, map, attr, &ret);
            kh_value(map, k) = strdup(value);
            kh_key(map, k) = attr_lib;
        }
    } else { // Attribute doesn't exist
        return -1;
    }

    return 0;
}

void slow5_hdr_free(struct slow5_hdr *header) {
    if (header) {
        slow5_hdr_data_free(header);
        slow5_aux_meta_free(header->aux_meta);

        free(header);
    }
}


/************************************* slow5 header data *************************************/

/*
 * read slow5 data header from file
 * return 0 on success, -1 on error
 * errors
 * SLOW5_ERR_HDRPARSE
 * SLOW5_ERR_MEM
 * SLOW5_ERR_OTH
 */
int slow5_hdr_data_init(FILE *fp, char **bufp, size_t *cap, struct slow5_hdr *header, uint32_t *hdr_len) {

    uint32_t hdr_len_tmp = 0;

    // Parse slow5 header data
    ssize_t buf_len;
    // Get first line of header data
    if ((buf_len = getline(bufp, cap, fp)) == -1) {
        SLOW5_ERROR("%s", "Malformed slow5 header data. No newline characters after number of read groups.");
        slow5_errno = SLOW5_ERR_HDRPARSE;
        goto err;
    }
    char *buf = *bufp;
    buf[buf_len - 1] = '\0'; // Remove newline for later parsing
    hdr_len_tmp += buf_len;

    /* init header data vector of hash maps; one map per read group */
    kv_init(header->data.maps);
    kv_resize(khash_t(slow5_s2s) *, header->data.maps, header->num_read_groups);
    if (!header->data.maps.a) { /* malloc error */
        SLOW5_MALLOC_ERROR()
        slow5_errno = SLOW5_ERR_MEM;
        goto err;
    }
    for (uint64_t i = 0; i < (uint64_t) header->num_read_groups; ++ i) {
        kv_A(header->data.maps, i) = kh_init(slow5_s2s);

        if (!kv_A(header->data.maps, i)) { /* calloc error */
            SLOW5_MALLOC_ERROR()
            slow5_errno = SLOW5_ERR_MEM;
            header->data.maps.n = i;
            goto err;
        }
    }
    header->data.maps.n = header->num_read_groups;

    /* init header data set of @[NAME] attributes */
    header->data.attrs = kh_init(slow5_s);
    if (!header->data.attrs) { /* calloc error */
        SLOW5_MALLOC_ERROR()
        slow5_errno = SLOW5_ERR_MEM;
        goto err;
    }


    // While the column header hasn't been reached
    while (strncmp(buf, SLOW5_ASCII_TYPE_HDR_MIN, strlen(SLOW5_ASCII_TYPE_HDR_MIN)) != 0) {

        // Ensure prefix is there
        if (buf[0] != SLOW5_HDR_DATA_PREFIX_CHAR) {
            SLOW5_ERROR("Malformed slow5 header data. Expected '%s...' or '%s...', instead found '%s'.",
                    SLOW5_HDR_DATA_PREFIX, SLOW5_ASCII_TYPE_HDR_MIN, buf); /* TODO truncate buf size in message */
            slow5_errno = SLOW5_ERR_HDRPARSE;
            goto err;
        }
        char *shift = buf + strlen(SLOW5_HDR_DATA_PREFIX); // Remove prefix

        // Get the attribute name
        char *attr = strdup(slow5_strsep(&shift, SLOW5_SEP_COL));

        int absent;
        kh_put(slow5_s, header->data.attrs, attr, &absent);
        if (absent == -1 || absent == 0) {
            if (absent == -1) {
                SLOW5_ERROR("Header data attribute '%s' failed to be inserted into attribute set.", attr);
                slow5_errno = SLOW5_ERR_OTH;
            }
            if (absent == 0) {
                SLOW5_ERROR("Malformed slow5 header data. Duplicate data attribute '%s'.", attr);
                slow5_errno = SLOW5_ERR_HDRPARSE;
            }
            free(attr);
            goto err;
        }
        ++ header->data.num_attrs;

        // Iterate through the values
        char *val;
        uint32_t i = 0;
        while ((val = slow5_strsep(&shift, SLOW5_SEP_COL)) != NULL && i < header->num_read_groups) {

            // Set key
            int absent;
            khint_t pos = kh_put(slow5_s2s, header->data.maps.a[i], attr, &absent);
            if (absent == -1) { /* ignoring 0 since duplicated attr should have been caught above */
                SLOW5_ERROR("Header data attribute '%s' failed to be inserted into attribute hash map.", attr);
                slow5_errno = SLOW5_ERR_OTH;
                free(attr);
                goto err;
            }

            // if the value is ".", we store an empty string
            if (strcmp(val, SLOW5_ASCII_MISSING) == 0) {
                val[0] = '\0';
            }

            char *val_dup = strdup(val);

            // Set value
            kh_val(header->data.maps.a[i], pos) = val_dup;

            ++ i;
        }
        // Ensure that read group number of entries are read
        if (i != header->num_read_groups || val) {
            SLOW5_ERROR("Malformed slow5 header data. Mismatch between number of entries (%s%" PRIu32 ") for attribute '%s' and number of read groups (%" PRIu32 ").",
                    val ? ">" : "", i, attr, header->num_read_groups);
            slow5_errno = SLOW5_ERR_HDRPARSE;
            goto err;
        }

        // Get next line
        if ((buf_len = getline(bufp, cap, fp)) == -1) {
            SLOW5_ERROR("Malformed slow5 header data. No more lines after '%s...'.", *bufp);
            slow5_errno = SLOW5_ERR_HDRPARSE;
            goto err;
        }
        buf = *bufp;
        buf[buf_len - 1] = '\0'; // Remove newline for later parsing
        hdr_len_tmp += buf_len;
    }

    if (hdr_len) {
        *hdr_len = hdr_len_tmp;
    }
    return 0;

    err:
        slow5_hdr_data_free(header);
        return -1;
}

struct slow5_aux_meta *slow5_aux_meta_init_empty(void) {
    struct slow5_aux_meta *aux_meta = (struct slow5_aux_meta *) calloc(1, sizeof *aux_meta);
    if(!aux_meta){
        SLOW5_MALLOC_ERROR()
        slow5_errno = SLOW5_ERR_MEM;
        return NULL;
    }

    aux_meta->cap = SLOW5_AUX_META_CAP_INIT;
    aux_meta->attrs = (char **) malloc(aux_meta->cap * sizeof *(aux_meta->attrs));
    aux_meta->types = (enum slow5_aux_type *) malloc(aux_meta->cap * sizeof *(aux_meta->types));
    aux_meta->sizes = (uint8_t *) malloc(aux_meta->cap * sizeof *(aux_meta->sizes));

    if(!aux_meta->attrs || !aux_meta->types || !aux_meta->sizes) {
        SLOW5_MALLOC_ERROR()
        slow5_errno = SLOW5_ERR_MEM;
        free(aux_meta->attrs);
        free(aux_meta->types);
        free(aux_meta->sizes);
        free(aux_meta);
        return NULL;
    }

    return aux_meta;
}

/*
 * parses the auxiliary data type and attribute row in the header (from file)
 * returns populated slow5_aux_meta
 * BE CAREFUL about using this (see assumptions below), it's intended for internal use in slow5_hdr_init
 *
 * assumptions:
 * only hdr_len can be NULL
 * *cap has the capacity of *bufp
 * already true that strcmp(*bufp, SLOW5_ASCII_TYPE_HDR_MIN) == 0
 * *bufp stores the header data type line
 *
 * returns -1 on error, 0 on success
 * errors:
 * SLOW5_ERR_MEM
 * SLOW5_ERR_HDRPARSE
 * SLOW5_ERR_OTH
 */
struct slow5_aux_meta *slow5_aux_meta_init(FILE *fp, char **bufp, size_t *cap, uint32_t *hdr_len, int *err) {

    struct slow5_aux_meta *aux_meta = NULL;
    char *buf = *bufp;
    ssize_t buf_len = strlen(buf);

    // Parse data types and deduce the sizes
    if (buf_len != strlen(SLOW5_ASCII_TYPE_HDR_MIN)) {
        char *shift = buf + strlen(SLOW5_ASCII_TYPE_HDR_MIN);

        char *tok = slow5_strsep(&shift, SLOW5_SEP_COL);
        if (strcmp(tok, "") != 0) {
            SLOW5_ERROR("Malformed slow5 header. Expected '%s...', instead found '%s%s'.",
                    SLOW5_ASCII_TYPE_HDR_MIN, buf, "..." ? shift : "");
            slow5_errno = SLOW5_ERR_HDRPARSE;
            goto err;
        }

        aux_meta = (struct slow5_aux_meta *) calloc(1, sizeof *aux_meta);
        if (!aux_meta) {
            SLOW5_MALLOC_ERROR();
            slow5_errno = SLOW5_ERR_MEM;
            goto err;
        }
        aux_meta->cap = SLOW5_AUX_META_CAP_INIT;
        aux_meta->types = (enum slow5_aux_type *) malloc(aux_meta->cap * sizeof *(aux_meta->types));
        aux_meta->sizes = (uint8_t *) malloc(aux_meta->cap * sizeof *(aux_meta->sizes));
        if (!aux_meta->sizes || !aux_meta->types) {
            SLOW5_MALLOC_ERROR();
            slow5_aux_meta_free(aux_meta);
            slow5_errno = SLOW5_ERR_MEM;
            goto err;
        }

        while ((tok = slow5_strsep(&shift, SLOW5_SEP_COL)) != NULL) {
            int err;
            enum slow5_aux_type type = slow5_str_to_aux_type(tok, &err);

            if (err == -1) {
                SLOW5_ERROR("Malformed slow5 header. Unknown type '%s' for auxiliary column %" PRIu32 ".",
                        tok, aux_meta->num + 1);
                slow5_aux_meta_free(aux_meta);
                slow5_errno = SLOW5_ERR_HDRPARSE;
                goto err;
            }

            aux_meta->types[aux_meta->num] = type;
            aux_meta->sizes[aux_meta->num] = SLOW5_AUX_TYPE_META[type].size;

            /* parse enum labels */
            if (type == SLOW5_ENUM || type == SLOW5_ENUM_ARRAY) {

                /* if hasn't been allocated yet */
                if (!aux_meta->enum_labels) {
                    aux_meta->enum_labels = (char ***) malloc(aux_meta->cap * sizeof *(aux_meta->enum_labels));
                    aux_meta->enum_num_labels = (uint8_t *) calloc(aux_meta->cap, sizeof *(aux_meta->enum_num_labels));
                    if (!aux_meta->enum_labels || !aux_meta->enum_num_labels) {
                        SLOW5_MALLOC_ERROR();
                        free(aux_meta->enum_labels);
                        free(aux_meta->enum_num_labels);
                        slow5_aux_meta_free(aux_meta);
                        slow5_errno = SLOW5_ERR_MEM;
                        goto err;
                    }
                }

                uint8_t enum_num_labels;
                char **enum_labels = slow5_aux_meta_enum_parse(tok, type, &enum_num_labels);
                if (!enum_labels) {
                    SLOW5_ERROR("Malformed slow5 header. Invalid enum labels '%s', expected 'enum%s{LABEL_0,LABEL_1,...}'.",
                            tok, SLOW5_IS_PTR(type) ? "*" : "");
                    slow5_aux_meta_free(aux_meta);
                    slow5_errno = SLOW5_ERR_HDRPARSE;
                    goto err;
                } else {
                    aux_meta->enum_num_labels[aux_meta->num] = enum_num_labels;
                    aux_meta->enum_labels[aux_meta->num] = enum_labels;
                }
            }

            ++ aux_meta->num;
            if (aux_meta->num > aux_meta->cap) { /* resize buffers capacity */
                size_t old_cap = aux_meta->cap;
                aux_meta->cap = aux_meta->cap << 1;

                enum slow5_aux_type *types_tmp = (enum slow5_aux_type *) realloc(aux_meta->types, aux_meta->cap * sizeof *(aux_meta->types));
                uint8_t *sizes_tmp = (uint8_t *) realloc(aux_meta->sizes, aux_meta->cap * sizeof *(aux_meta->sizes));
                if (!types_tmp || !sizes_tmp) { /* realloc error */
                    SLOW5_MALLOC_ERROR();
                    free(types_tmp);
                    free(sizes_tmp);
                    slow5_aux_meta_free(aux_meta);
                    slow5_errno = SLOW5_ERR_MEM;
                    goto err;
                }
                aux_meta->types = types_tmp;
                aux_meta->sizes = sizes_tmp;

                /* if enum labels allocated reallocate them as well */
                if (aux_meta->enum_labels) {
                    char ***enum_labels_tmp = (char ***) realloc(aux_meta->enum_labels, aux_meta->cap * sizeof *(aux_meta->enum_labels));
                    uint8_t *enum_num_labels_tmp = (uint8_t *) realloc(aux_meta->enum_num_labels, aux_meta->cap * sizeof *(aux_meta->enum_num_labels));
                    memset(enum_num_labels_tmp + old_cap, 0, aux_meta->cap - old_cap);
                    if (!enum_labels_tmp || !enum_num_labels_tmp) { /* realloc error */
                        SLOW5_MALLOC_ERROR();
                        free(enum_labels_tmp);
                        free(enum_num_labels_tmp);
                        slow5_aux_meta_free(aux_meta);
                        slow5_errno = SLOW5_ERR_MEM;
                        goto err;
                    }
                    aux_meta->enum_labels = enum_labels_tmp;
                    aux_meta->enum_num_labels = enum_num_labels_tmp;
                }
            }
        }
    }

    // Get column names
    if ((buf_len = getline(bufp, cap, fp)) == -1) {
        SLOW5_ERROR("Malformed slow5 header. No more lines after '%s...'.", *bufp);
        slow5_aux_meta_free(aux_meta);
        slow5_errno = SLOW5_ERR_HDRPARSE;
        goto err;
    }
    buf = *bufp;
    if (hdr_len) {
        *hdr_len += buf_len;
    }
    buf[-- buf_len] = '\0'; // Remove newline for later parsing

    /* Check mandatory columns are there */
    if (strncmp(buf, SLOW5_ASCII_COLUMN_HDR_MIN, strlen(SLOW5_ASCII_COLUMN_HDR_MIN)) != 0) {
        SLOW5_ERROR("Malformed slow5 header. Expected '%s...', instead found '%s...'.",
                SLOW5_ASCII_COLUMN_HDR_MIN, buf);
        slow5_aux_meta_free(aux_meta);
        slow5_errno = SLOW5_ERR_HDRPARSE;
        goto err;
    }

    // Parse auxiliary attributes
    if (buf_len != strlen(SLOW5_ASCII_COLUMN_HDR_MIN)) {
        char *shift = buf + strlen(SLOW5_ASCII_COLUMN_HDR_MIN);

        char *tok = slow5_strsep(&shift, SLOW5_SEP_COL);
        if (strcmp(tok, "") != 0) {
            SLOW5_ERROR("Malformed slow5 header. Expected '%s...', instead found '%s...'.",
                    SLOW5_ASCII_COLUMN_HDR_MIN, buf);
            slow5_errno = SLOW5_ERR_HDRPARSE;
            slow5_aux_meta_free(aux_meta);
            goto err;
        }
        if (!aux_meta) {
            SLOW5_ERROR("%s", "Malformed slow5 header. Found auxiliary columns but missing auxiliary types.");
            slow5_errno = SLOW5_ERR_HDRPARSE;
            goto err;
        }

        aux_meta->attr_to_pos = kh_init(slow5_s2ui32);
        aux_meta->attrs = (char **) malloc(aux_meta->cap * sizeof *(aux_meta->attrs));
        if (!aux_meta->attr_to_pos || !aux_meta->attrs) {
            SLOW5_MALLOC_ERROR();
            aux_meta->num = 0;
            slow5_aux_meta_free(aux_meta);
            slow5_errno = SLOW5_ERR_MEM;
            goto err;
        }

        for (uint64_t i = 0; i < aux_meta->num; ++ i) {
            if ((tok = slow5_strsep(&shift, SLOW5_SEP_COL)) == NULL) {
                SLOW5_ERROR("Malformed slow5 header. Mismatch between the number of auxiliary types (%" PRIu32 ") and auxiliary attributes (%" PRIu64 ").",
                        aux_meta->num, i);
                aux_meta->num = i;
                slow5_aux_meta_free(aux_meta);
                slow5_errno = SLOW5_ERR_HDRPARSE;
                goto err;
            }

            aux_meta->attrs[i] = strdup(tok);
            if (!aux_meta->attrs[i]) {
                SLOW5_MALLOC_ERROR();
                aux_meta->num = i;
                slow5_aux_meta_free(aux_meta);
                slow5_errno = SLOW5_ERR_MEM;
                goto err;
            }

            int absent;
            khint_t pos = kh_put(slow5_s2ui32, aux_meta->attr_to_pos, aux_meta->attrs[i], &absent);
            if (absent == -1 || absent == 0) {
                if (absent == -1) {
                    SLOW5_ERROR("Header auxiliary attribute '%s' failed to be inserted into auxiliary attribute hash map.",
                            aux_meta->attrs[i]);
                    slow5_errno = SLOW5_ERR_OTH;
                }
                if (absent == 0) {
                    SLOW5_ERROR("Malformed slow5 header data. Duplicate auxiliary attribute '%s'.",
                            aux_meta->attrs[i]);
                    slow5_errno = SLOW5_ERR_HDRPARSE;
                }
                aux_meta->num = i + 1;
                slow5_aux_meta_free(aux_meta);
                slow5_errno = SLOW5_ERR_HDRPARSE;
                goto err;
            }
            kh_value(aux_meta->attr_to_pos, pos) = i;
        }
        if ((tok = slow5_strsep(&shift, SLOW5_SEP_COL)) != NULL) {
            SLOW5_ERROR("Malformed slow5 header. Mismatch between the number of auxiliary types (%" PRIu32 ") and auxiliary attributes (%s%" PRIu32 ").",
                    aux_meta->num, shift ? ">" : "", aux_meta->num + 1);
            slow5_aux_meta_free(aux_meta);
            goto err;
        }
    } else if (aux_meta) {
        SLOW5_ERROR("%s", "Malformed slow5 header. Found auxiliary types but missing auxiliary columns.");
        slow5_aux_meta_free(aux_meta);
        slow5_errno = SLOW5_ERR_HDRPARSE;
        goto err;
    }

    *err = 0;
    return aux_meta;

    err:
        *err = -1;
        return NULL;
}

/**
 * parses "enum[*]{LABEL_0,LABEL_1,...}"
 * expects the first enum[*] part to be verified
 * tok, n cannot be NULL
 * rules: no spaces anywhere, no empty labels, no duplicate labels
 * returns the array of labels as strings
 * returns NULL on parsing error
 * free the labels and the array after use
 */
char **slow5_aux_meta_enum_parse(char *tok, enum slow5_aux_type type, uint8_t *n) {

    const char *type_str = SLOW5_AUX_TYPE_META[type].type_str;
    const size_t tok_len = strlen(tok);

    if (tok_len == strlen(type_str)) {
        SLOW5_ERROR("Expected '%c' after '%s'.",
                SLOW5_HDR_ENUM_LABELS_BEGIN, type_str);
        return NULL;
    }
    if (tok[strlen(type_str)] != SLOW5_HDR_ENUM_LABELS_BEGIN) {
        SLOW5_ERROR("Expected '%c' after '%s', instead of '%c'.",
                SLOW5_HDR_ENUM_LABELS_BEGIN, type_str, tok[strlen(type_str)]);
        return NULL;
    }
    if (tok[tok_len - 1] != SLOW5_HDR_ENUM_LABELS_END) {
        SLOW5_ERROR("Expected '%c' at the end of '%s'.",
                SLOW5_HDR_ENUM_LABELS_END, tok);
        return NULL;
    }

    tok[tok_len - 1] = '\0'; /* ignore trailing } */
    tok += strlen(type_str) + 1; /* skip "enum[*]{" */

    uint8_t num = 0;
    uint16_t cap = SLOW5_AUX_ENUM_LABELS_CAP_INIT;
    char **labels = malloc(cap * sizeof *labels);
    if (!labels) {
        SLOW5_MALLOC_ERROR();
        return NULL;
    }

    char *label = slow5_strsep(&tok, SLOW5_SEP_ARRAY);
    do {
        /* check if it's a std c label */
        int ret;
        if ((ret = slow5_is_c_label(label)) != 0) {

            if (ret == -1) {
                SLOW5_ERROR("Enum label for %d is empty.", num);
            } else if (ret == -2) {
                SLOW5_ERROR("Enum label '%s' for %d has a character which is not alpha-numeric nor an underscore.", label, num);
            } else { /* ret == -3 */
                SLOW5_ERROR("Enum label '%s' for %d starts with a number.", label, num);
            }

            for (uint16_t i = 0; i < num; ++ i) {
                free(labels[i]);
            }
            free(labels);
            return NULL;
        }

        /* very naive algorithm TODO could use hashsets */
        for (uint16_t i = 0; i < num; ++ i) {
            if (strcmp(label, labels[i]) == 0) { /* duplicate label */
                SLOW5_ERROR("Enum label '%s' for %d is a duplicate of label %d.",
                        label, num, i);

                for (uint16_t i = 0; i < num; ++ i) {
                    free(labels[i]);
                }
                free(labels);
                return NULL;
            }
        }

        char *label_cp = strdup(label);
        if (!label_cp) {
            SLOW5_MALLOC_ERROR();
            for (uint16_t i = 0; i < num; ++ i) {
                free(labels[i]);
            }
            free(labels);
            return NULL;
        }

        if (num >= cap) {
            cap = cap << 1;
            char **labels_tmp = realloc(labels, cap * sizeof *labels);
            if (!labels_tmp) { /* memory allocation error */
                SLOW5_MALLOC_ERROR();
                for (uint16_t i = 0; i < num; ++ i) {
                    free(labels[i]);
                }
                free(labels);
                return NULL;
            }
            labels = labels_tmp;
        }

        labels[num] = label_cp;
        ++ num;
    } while ((label = slow5_strsep(&tok, SLOW5_SEP_ARRAY)) != NULL);

    *n = num; /* n better not be NULL */

    return labels;
}

// Return
// 0    success
// -1   null input
// -2   other failure
// -3   use slow5_aux_meta_add_enum instead if type is SLOW5_ENUM or SLOW5_ENUM_ARRAY
/* TODO this error checking is bad */
int slow5_aux_meta_add(struct slow5_aux_meta *aux_meta, const char *attr, enum slow5_aux_type type) {
    if (aux_meta == NULL || attr == NULL) {
        return -1;
    }

    if (type == SLOW5_ENUM || type == SLOW5_ENUM_ARRAY) {
        return -3;
    }

    if (aux_meta->attr_to_pos == NULL) {
        aux_meta->attr_to_pos = kh_init(slow5_s2ui32);
    }

    if (aux_meta->num == aux_meta->cap) {
        aux_meta->cap = aux_meta->cap << 1; // TODO is this ok?
        aux_meta->attrs = (char **) realloc(aux_meta->attrs, aux_meta->cap * sizeof *(aux_meta->attrs));
        SLOW5_MALLOC_CHK(aux_meta->attrs);
        aux_meta->types = (enum slow5_aux_type *) realloc(aux_meta->types, aux_meta->cap * sizeof *(aux_meta->types));
        SLOW5_MALLOC_CHK(aux_meta->types);
        aux_meta->sizes = (uint8_t *) realloc(aux_meta->sizes, aux_meta->cap * sizeof *(aux_meta->sizes));
        SLOW5_MALLOC_CHK(aux_meta->sizes);
    }

    aux_meta->attrs[aux_meta->num] = strdup(attr);

    int absent;
    khint_t pos = kh_put(slow5_s2ui32, aux_meta->attr_to_pos, aux_meta->attrs[aux_meta->num], &absent);
    if (absent == -1 || absent == 0) {
        free(aux_meta->attrs[aux_meta->num]);
        return -2;
    }
    kh_value(aux_meta->attr_to_pos, pos) = aux_meta->num;

    aux_meta->types[aux_meta->num] = type;
    aux_meta->sizes[aux_meta->num] = SLOW5_AUX_TYPE_META[type].size;

    ++ aux_meta->num;

    return 0;
}

int slow5_aux_add(const char *field, enum slow5_aux_type type, slow5_hdr_t *header){
    int ret = slow5_aux_meta_add(header->aux_meta, field, type);
    return ret;
}


// Return
// 0    success
// -1   null input
// -2   other failure
// -3   use slow5_aux_meta_add instead if type is not SLOW5_ENUM or SLOW5_ENUM_ARRAY
// -4   bad enum_labels not good c labels
/* TODO this error checking is bad */
int slow5_aux_meta_add_enum(struct slow5_aux_meta *aux_meta, const char *attr, enum slow5_aux_type type, const char **enum_labels, uint8_t enum_num_labels) {
    if (!aux_meta || !attr || !enum_labels) {
        return -1;
    }

    if (type != SLOW5_ENUM && type != SLOW5_ENUM_ARRAY) {
        return -3;
    }

    if (aux_meta->attr_to_pos == NULL) {
        aux_meta->attr_to_pos = kh_init(slow5_s2ui32);
    }

    /* hard copy enum_labels */
    char **enum_labels_cp = malloc(enum_num_labels * sizeof *enum_labels_cp);
    if (!enum_labels_cp) { /* mem allocation error */
        SLOW5_MALLOC_ERROR();
        return -2;
    }
    for (uint16_t i = 0; i < enum_num_labels; ++ i) {
        /* check enum label is valid */
        if (slow5_is_c_label(enum_labels[i]) != 0) {
            SLOW5_ERROR("Bad enum label '%s' for %d", enum_labels[i], i);
            for (uint16_t j = 0; j < i; ++ j) {
                free(enum_labels_cp[j]);
            }
            free(enum_labels_cp);
            return -4;
        }

        enum_labels_cp[i] = strdup(enum_labels[i]);

        if (!enum_labels_cp[i]) { /* mem allocation error */
            SLOW5_MALLOC_ERROR();
            for (uint16_t j = 0; j < i; ++ j) {
                free(enum_labels_cp[j]);
            }
            free(enum_labels_cp);
            return -2;
        }
    }

    if (aux_meta->num == aux_meta->cap) {
        size_t old_cap = aux_meta->cap;
        aux_meta->cap = old_cap << 1;

        /* TODO this reallocing is bad */
        aux_meta->attrs = (char **) realloc(aux_meta->attrs, aux_meta->cap * sizeof *(aux_meta->attrs));
        SLOW5_MALLOC_CHK(aux_meta->attrs);
        aux_meta->types = (enum slow5_aux_type *) realloc(aux_meta->types, aux_meta->cap * sizeof *(aux_meta->types));
        SLOW5_MALLOC_CHK(aux_meta->types);
        aux_meta->sizes = (uint8_t *) realloc(aux_meta->sizes, aux_meta->cap * sizeof *(aux_meta->sizes));
        SLOW5_MALLOC_CHK(aux_meta->sizes);
        if (aux_meta->enum_labels) {
            aux_meta->enum_labels = (char ***) realloc(aux_meta->enum_labels, aux_meta->cap * sizeof *(aux_meta->enum_labels));
            SLOW5_MALLOC_CHK(aux_meta->enum_labels);
            aux_meta->enum_num_labels = (uint8_t *) realloc(aux_meta->enum_num_labels, aux_meta->cap * sizeof *(aux_meta->enum_num_labels));
            SLOW5_MALLOC_CHK(aux_meta->enum_num_labels);
            memset(aux_meta->enum_num_labels + old_cap, 0, aux_meta->cap - old_cap);
        }
    }

    /* enum labels not allocated yet */
    if (!aux_meta->enum_labels) {
        aux_meta->enum_labels = (char ***) malloc(aux_meta->cap * sizeof *(aux_meta->enum_labels));
        aux_meta->enum_num_labels = (uint8_t *) calloc(aux_meta->cap, sizeof *(aux_meta->enum_num_labels));
        if (!aux_meta->enum_labels || !aux_meta->enum_num_labels) {
            SLOW5_MALLOC_ERROR();
            free(aux_meta->enum_labels);
            free(aux_meta->enum_num_labels);
            for (uint16_t i = 0; i < enum_num_labels; ++ i) {
                free(enum_labels_cp[i]);
            }
            free(enum_labels_cp);
            return -2;
        }
    }

    aux_meta->attrs[aux_meta->num] = strdup(attr);
    aux_meta->enum_labels[aux_meta->num] = (char **) enum_labels_cp;
    aux_meta->enum_num_labels[aux_meta->num] = enum_num_labels;

    int absent;
    khint_t pos = kh_put(slow5_s2ui32, aux_meta->attr_to_pos, aux_meta->attrs[aux_meta->num], &absent);
    if (absent == -1 || absent == 0) {
        free(aux_meta->attrs[aux_meta->num]);
        for (uint16_t i = 0; i < enum_num_labels; ++ i) {
            free(enum_labels_cp[i]);
        }
        free(enum_labels_cp);
        return -2;
    }
    kh_value(aux_meta->attr_to_pos, pos) = aux_meta->num;

    aux_meta->types[aux_meta->num] = type;
    aux_meta->sizes[aux_meta->num] = SLOW5_AUX_TYPE_META[type].size;

    ++ aux_meta->num;

    return 0;
}

void slow5_aux_meta_free(struct slow5_aux_meta *aux_meta) {
    if (aux_meta) {
        if (aux_meta->attrs) {
            for (uint64_t i = 0; i < aux_meta->num; ++ i) {
                free(aux_meta->attrs[i]);
            }
            free(aux_meta->attrs);
        }
        kh_destroy(slow5_s2ui32, aux_meta->attr_to_pos);
        free(aux_meta->types);
        free(aux_meta->sizes);
        if (aux_meta->enum_labels) {
            for (uint64_t i = 0; i < aux_meta->num; ++ i) {
                for (uint16_t j = 0; j < aux_meta->enum_num_labels[i]; ++ j) {
                    free((void *) aux_meta->enum_labels[i][j]);
                }

                if (aux_meta->enum_num_labels[i] > 0) {
                    free(aux_meta->enum_labels[i]);
                }
            }

            free(aux_meta->enum_labels);
            free(aux_meta->enum_num_labels);
        }
        free(aux_meta);
    }
}

void slow5_hdr_data_free(struct slow5_hdr *header) {

    if (header) {
        if (header->data.attrs && header->data.maps.a) {

            for (khint_t i = kh_begin(header->data.attrs); i < kh_end(header->data.attrs); ++ i) {
                if (kh_exist(header->data.attrs, i)) {
                    char *attr = (char *) kh_key(header->data.attrs, i);

                    // Free header data map
                    for (size_t j = 0; j < kv_size(header->data.maps); ++ j) {
                        khash_t(slow5_s2s) *map = header->data.maps.a[j];

                        khint_t pos = kh_get(slow5_s2s, map, attr);
                        if (pos != kh_end(map)) {
                            free(kh_value(map, pos));
                            kh_del(slow5_s2s, map, pos);
                        }
                    }

                    free(attr);
                }
            }
        }

        // Free header data map
        for (size_t j = 0; j < kv_size(header->data.maps); ++ j) {
            kh_destroy(slow5_s2s, header->data.maps.a[j]);
        }
        kv_destroy(header->data.maps);
        kh_destroy(slow5_s, header->data.attrs);
    }
}


// slow5 record

/*
 * get slow5 record with read_id as it is stored in the file s5p with length *n
 * on error returns NULL, sets *n=0 if possible, and sets slow5_errno
 * slow5_errno errors:
 * SLOW5_ERR_ARG
 * SLOW5_ERR_NOIDX
 * SLOW5_ERR_NOTFOUND
 * SLOW5_ERR_UNK
 * SLOW5_ERR_MEM
 * SLOW5_ERR_IO
 */
void *slow5_get_mem(const char *read_id, size_t *n, const struct slow5_file *s5p) {
    if (!read_id || !s5p) {
        if (!read_id) {
            SLOW5_ERROR("Argument '%s' cannot be NULL.", SLOW5_TO_STR(read_id));
        }
        if (!s5p) {
            SLOW5_ERROR("Argument '%s' cannot be NULL.", SLOW5_TO_STR(s5p));
        }
        slow5_errno = SLOW5_ERR_ARG;
        goto err;
    }

    if (!s5p->index) {
        /* index not loaded */
        SLOW5_ERROR("%s", "No slow5 index has been loaded.");
        slow5_errno = SLOW5_ERR_NOIDX;
        goto err;
    }

    /* get index record */
    struct slow5_rec_idx read_index;
    if (slow5_idx_get(s5p->index, read_id, &read_index) == -1) {
        /* read_id not found in index */
        slow5_errno = SLOW5_ERR_NOTFOUND;
        goto err;
    }

    size_t bytes;
    off_t offset;
    if (s5p->format == SLOW5_FORMAT_BINARY) {
        bytes = read_index.size - sizeof (slow5_rec_size_t);
        offset = read_index.offset + sizeof (slow5_rec_size_t);
    } else if (s5p->format == SLOW5_FORMAT_ASCII) {
        bytes = read_index.size;
        offset = read_index.offset;
    } else {
        SLOW5_ERROR("Unknown slow5 format '%d'.", s5p->format);
        slow5_errno = SLOW5_ERR_UNK;
        goto err;
    }

    uint8_t *mem = (uint8_t *) malloc(bytes);
    if (!mem) {
        SLOW5_MALLOC_ERROR();
        slow5_errno = SLOW5_ERR_MEM;
        goto err;
    }

    if (s5p->format == SLOW5_FORMAT_ASCII) {
        bytes -= 1;
        /* null terminate */
        mem[bytes] = '\0';
    }

    if (pread(s5p->meta.fd, mem, bytes, offset) != bytes) {
        SLOW5_ERROR("Failed to pread '%zu' bytes at offset '%zu' from slow5 file '%s'.",
                bytes, offset, s5p->meta.pathname);
        free(mem);
        slow5_errno = SLOW5_ERR_IO;
        goto err;
    }

    if (n) {
        *n = bytes;
    }
    return mem;

    err:
        if (n) {
            *n = 0;
        }
        return NULL;
}


/**
 * slow5_get
 * Get a read entry from a slow5 file corresponding to a read_id.
 *
 * Allocates memory for *read if it is NULL.
 * Otherwise, the data in *read is freed and overwritten.
 * slow5_rec_free() should always be called when finished with the structure.
 *
 * Require the slow5 index to be loaded beforehand using slow5_idx_load()
 *
 * Return:
 *  =0   the read was successfully found and stored
 *  <0   error code
 *
 * Errors:
 * SLOW5_ERR_ARG        read_id, read or s5p is NULL
 * SLOW5_ERR_NOTFOUND   read_id was not found in the index
 * SLOW5_ERR_IO         other error when reading the slow5 file
 * SLOW5_ERR_RECPARSE   record parsing error
 * SLOW5_ERR_NOIDX      the index has not been loaded
 * SLOW5_ERR_MEM        memory allocation error
 * SLOW5_ERR_UNK        slow5 file format is unknown
 * SLOW5_ERR_PRESS      decompression error
 *
 * @param   read_id the read identifier
 * @param   read    address of a slow5_rec pointer
 * @param   s5p     slow5 file
 * @return  error code described above
 */
int slow5_get(const char *read_id, struct slow5_rec **read, struct slow5_file *s5p) {

    if (!read) {
        SLOW5_ERROR_EXIT("Argument '%s' cannot be NULL.", SLOW5_TO_STR(read));
        return slow5_errno = SLOW5_ERR_ARG;
    }

    size_t bytes;
    char *mem;
    if (!(mem = slow5_get_mem(read_id, &bytes, s5p))) {
        SLOW5_EXIT_IF_ON_ERR();
        return slow5_errno;
    }

    if (slow5_rec_depress_parse(&mem, &bytes, read_id, read, s5p) != 0) {
        SLOW5_EXIT_IF_ON_ERR();
        free(mem);
        return slow5_errno;
    }

    free(mem);
    return 0;
}


//gets the list of read ids from the SLOW5 index
//the list of read is is a pointer and must not be freed by user
//*len will have the number of read ids
//NULL will be returned in ase of error
char **slow5_get_rids(const slow5_file_t *s5p, uint64_t *len) {

    if (!s5p->index) {
        /* index not loaded */
        SLOW5_ERROR("%s", "No slow5 index has been loaded.");
        slow5_errno = SLOW5_ERR_NOIDX;
        return NULL;
        *len=0;
    }

    if(!s5p->index->ids){
        SLOW5_ERROR("%s", "No read ID list in the index.");
        slow5_errno = SLOW5_ERR_OTH;
        return NULL;
        *len=0;
    }

    *len = s5p->index->num_ids;
    return s5p->index->ids;

}


/*
 * decompress record if s5p has a compression method then parse to read
 * set mem to decompressed mem and bytes to new bytes
 * mem, bytes, s5p cannot be null
 * return 0 on success
 * return a SLOW5_ERR_* and set slow5_errno on failure
 * doesn't free *mem
 */
int slow5_rec_depress_parse(char **mem, size_t *bytes, const char *read_id, struct slow5_rec **read, struct slow5_file *s5p) {

    /* assuming that if compress is initialised so is record_press */
    if (s5p->compress && s5p->compress->record_press->method != SLOW5_COMPRESS_NONE) {
        size_t new_bytes;
        char *new_mem;
        if (!(new_mem = slow5_ptr_depress_solo(s5p->compress->record_press->method, *mem, *bytes, &new_bytes)) || new_bytes == 0) {
            if (read_id) {
                SLOW5_ERROR("Failed to decompress read with ID '%s' from slow5 file '%s'.",
                        read_id, s5p->meta.pathname);
            } else {
                SLOW5_ERROR("Failed to decompress read from slow5 file '%s'.", s5p->meta.pathname);
            }
            return slow5_errno = SLOW5_ERR_PRESS;
        }
        free(*mem);
        *bytes = new_bytes;
        *mem = new_mem;
    }

    enum slow5_press_method signal_comp = SLOW5_COMPRESS_NONE;
    if (s5p->compress) {
        signal_comp = s5p->compress->signal_press->method;
    }

    if (slow5_rec_parse(*mem, *bytes, read_id, read, s5p->format, s5p->header->aux_meta, signal_comp) == -1) {
        SLOW5_ERROR("%s", "Record parsing failed.");
        return slow5_errno = SLOW5_ERR_RECPARSE;
    }

    return 0;
}

int slow_decode(void **mem, size_t *bytes, slow5_rec_t **read, slow5_file_t *s5p){
    return slow5_rec_depress_parse((char **)mem, bytes, NULL, read, s5p);
}

/*
 * parse read_mem with read_size bytes, intended read ID read_id, the given format and auxiliary meta data aux_meta
 * if read compression is used read_mem should be decompressed beforehand, signal decompression is handled here though
 * populate the read
 * read_mem, read cannot be NULL
 * returns -1 on error, 0 on success
 */
int slow5_rec_parse(char *read_mem, size_t read_size, const char *read_id, struct slow5_rec **readp, enum slow5_fmt format, struct slow5_aux_meta *aux_meta, enum slow5_press_method signal_method) {

    if (!*readp) {
        /* allocate memory for read */
        *readp = (struct slow5_rec *) calloc(1, sizeof **readp);
        if (!*readp) {
            SLOW5_MALLOC_ERROR();
            return -1;
        }
    } else {
        /* free previously allocated memory except for raw_signal */
        free((*readp)->read_id);
        (*readp)->read_id = NULL;
        slow5_rec_aux_free((*readp)->aux_map);
        (*readp)->aux_map = NULL;
    }

    struct slow5_rec *read = *readp;

    int ret = 0;
    uint8_t i = 0;
    const uint64_t prev_len_raw_signal = read->len_raw_signal;

    if (format == SLOW5_FORMAT_ASCII) {

        char *tok = slow5_strsep(&read_mem, SLOW5_SEP_COL);
        if (read_mem == NULL) {
            SLOW5_ERROR("Read ID could not be parsed from '%.10s%s'. Tab separator was not found.",
                    tok, strnlen(tok, 11) == 11 ? "..." : "");
            return -1;
        }

        bool more_to_parse = false;
        int err;
        do {
            switch (i) {
                case COL_read_id:
                    /* Ensure line matches requested id */
                    if (read_id != NULL) {
                        if (strcmp(tok, read_id) != 0) {
                            SLOW5_ERROR("Requested read ID '%s' does not match the read ID '%s' in fetched record.", read_id, tok);
                            ret = -1;
                            break;
                        }
                    }
                    read->read_id_len = strlen(tok);
                    read->read_id = strndup(tok, read->read_id_len);
                    if (read->read_id == NULL) {
                        if (errno == ENOMEM) {
                            SLOW5_ERROR("%s", "Insufficient memory available to allocate read ID.")
                        } else {
                            SLOW5_ERROR("Duplicating read ID failed: %s", strerror(errno));
                        }
                        ret = -1;
                    }
                    break;

                case COL_read_group:
                    /* TODO check slow5_ato_uint32 and the other ones below */
                    read->read_group = slow5_ato_uint32(tok, &err);
                    if (err == -1) {
                        SLOW5_ERROR("Read group could not be parsed from '%s'.", tok);
                        ret = -1;
                    }
                    break;

                case COL_digitisation:
                    read->digitisation = slow5_strtod_check(tok, &err);
                    if (err == -1) {
                        SLOW5_ERROR("Digitisation could not be parsed from '%s'.", tok);
                        ret = -1;
                    }
                    break;

                case COL_offset:
                    read->offset = slow5_strtod_check(tok, &err);
                    if (err == -1) {
                        SLOW5_ERROR("Offset could not be parsed from '%s'.", tok);
                        ret = -1;
                    }
                    break;

                case COL_range:
                    read->range = slow5_strtod_check(tok, &err);
                    if (err == -1) {
                        SLOW5_ERROR("Range could not be parsed from '%s'.", tok);
                        ret = -1;
                    }
                    break;

                case COL_sampling_rate:
                    read->sampling_rate = slow5_strtod_check(tok, &err);
                    if (err == -1) {
                        SLOW5_ERROR("Sampling rate could not be parsed from '%s'.", tok);
                        ret = -1;
                    }
                    break;

                case COL_len_raw_signal:
                    read->len_raw_signal = slow5_ato_uint64(tok, &err);
                    if (err == -1) {
                        SLOW5_ERROR("Raw signal length could not be parsed from '%s'.", tok);
                        ret = -1;
                    }
                    break;

                case COL_raw_signal: {
                    if (read->len_raw_signal == 0) {
                        break;
                    }

                    if (read->raw_signal == NULL) {
                        read->raw_signal = (int16_t *) malloc(read->len_raw_signal * sizeof *(read->raw_signal));
                        SLOW5_MALLOC_CHK(read->raw_signal);
                        if (read->raw_signal == NULL) {
                            read->len_raw_signal = 0;
                            ret = -1;
                            break;
                        }
                    } else if (prev_len_raw_signal < read->len_raw_signal) {
                        int16_t *raw_signal_tmp = (int16_t *) realloc(read->raw_signal, read->len_raw_signal * sizeof *(read->raw_signal));
                        SLOW5_MALLOC_CHK(raw_signal_tmp);
                        if (raw_signal_tmp == NULL) {
                            read->len_raw_signal = prev_len_raw_signal;
                            ret = -1;
                            break;
                        }
                        read->raw_signal = raw_signal_tmp;
                    }

                    char *signal_tok = slow5_strsep(&tok, SLOW5_SEP_ARRAY);
                    if (tok == NULL && read->len_raw_signal > 1) {
                        SLOW5_ERROR("Raw signals could not be parsed from '%.10s%s'. Comma separator was not found.",
                                signal_tok, strnlen(signal_tok, 11) == 11 ? "..." : "");
                        ret = -1;
                        break;
                    }

                    uint64_t j = 0;

                    /* Parse raw signal */
                    do {
                        (read->raw_signal)[j] = slow5_ato_int16(signal_tok, &err);
                        if (err == -1) {
                            SLOW5_ERROR("Raw signal '%" PRIu64 "' could not be parsed from '%s'.", j, signal_tok);
                            ret = -1;
                            break;
                        }
                        ++ j;
                    } while ((signal_tok = slow5_strsep(&tok, SLOW5_SEP_ARRAY)) != NULL);

                    if (ret != -1 && j != read->len_raw_signal) {
                        SLOW5_ERROR("Raw signal length '%" PRIu64 "' does not match the number of signals parsed '%" PRIu64 "'. Raw signal is potentially truncated.", read->len_raw_signal, j);
                        ret = -1;
                    }

                } break;

                /* All columns parsed */
                default:
                    more_to_parse = true;
                    break;

            }
            ++ i;

        } while (ret != -1 &&
                 !more_to_parse &&
                 (tok = slow5_strsep(&read_mem, SLOW5_SEP_COL)) != NULL);

        if (i < SLOW5_COLS_NUM) {
            /* Not all main columns parsed */
            SLOW5_ERROR("%" PRIu8 "/%" PRIu32 " slow5 main columns were parsed.", i, SLOW5_COLS_NUM);
            ret = -1;
        } else if (!more_to_parse && aux_meta != NULL) {
            SLOW5_ERROR("%s", "Missing auxiliary fields in record, but present in header.");
            ret = -1;
        } else if (more_to_parse) {
            /* Main columns parsed and more to parse */
            if (aux_meta == NULL) {
                SLOW5_ERROR("%s", "Missing auxiliary fields in header, but present in record.");
                ret = -1;
            } else {
                ret = slow5_rec_aux_parse(tok, read_mem, 0, read_size, read, format, aux_meta);
            }
        }

    } else if (format == SLOW5_FORMAT_BINARY) {

        bool main_cols_parsed = false;

        size_t size = 0;
        uint64_t offset = 0;

        while (!main_cols_parsed && ret != -1) {

            switch (i) {

                case COL_read_id:
                    size = sizeof read->read_id_len;
                    memcpy(&read->read_id_len, read_mem + offset, size);
                    offset += size;

                    size = read->read_id_len * sizeof *read->read_id;
                    read->read_id = strndup(read_mem + offset, size);
                    if (read->read_id == NULL) {
                        if (errno == ENOMEM) {
                            SLOW5_ERROR("%s", "Insufficient memory available to allocate read ID.")
                        } else {
                            SLOW5_ERROR("Duplicating read ID failed: %s", strerror(errno));
                        }
                        ret = -1;
                        break;
                    }

                    /* Ensure line matches requested id */
                    if (read_id != NULL && strcmp(read->read_id, read_id) != 0) {
                        SLOW5_ERROR("Requested read ID '%s' does not match the read ID '%s' in fetched record.", read_id, read->read_id);
                        ret = -1;
                        break;
                    }
                    offset += size;
                    break;

                case COL_read_group:
                    size = sizeof read->read_group;
                    memcpy(&read->read_group, read_mem + offset, size);
                    offset += size;
                    break;

                case COL_digitisation:
                    size = sizeof read->digitisation;
                    memcpy(&read->digitisation, read_mem + offset, size);
                    offset += size;
                    break;

                case COL_offset:
                    size = sizeof read->offset;
                    memcpy(&read->offset, read_mem + offset, size);
                    offset += size;
                    break;

                case COL_range:
                    size = sizeof read->range;
                    memcpy(&read->range, read_mem + offset, size);
                    offset += size;
                    break;

                case COL_sampling_rate:
                    size = sizeof read->sampling_rate;
                    memcpy(&read->sampling_rate, read_mem + offset, size);
                    offset += size;
                    break;

                case COL_len_raw_signal:
                    size = sizeof read->len_raw_signal;
                    memcpy(&read->len_raw_signal, read_mem + offset, size);
                    offset += size;
                    break;

                case COL_raw_signal: {
                    if (signal_method == SLOW5_COMPRESS_NONE) {
                        size = read->len_raw_signal * sizeof *(read->raw_signal);
                    } else { /* integer encoding */
                        size = read->len_raw_signal;
                    }
                    if (read->raw_signal == NULL) {
                        read->raw_signal = (int16_t *) malloc(size+16);
                        SLOW5_MALLOC_CHK(read->raw_signal);
                        if (read->raw_signal == NULL) {
                            read->len_raw_signal = 0;
                            ret = -1;
                            break;
                        }
                    } else if (prev_len_raw_signal < read->len_raw_signal) {
                        int16_t *raw_signal_tmp = (int16_t *) realloc(read->raw_signal, size+16);
                        SLOW5_MALLOC_CHK(raw_signal_tmp);
                        if (raw_signal_tmp == NULL) {
                            read->len_raw_signal = prev_len_raw_signal;
                            ret = -1;
                            break;
                        }
                        read->raw_signal = raw_signal_tmp;
                    }

                    memcpy(read->raw_signal, read_mem + offset, size);
                    offset += size;

                    if (signal_method != SLOW5_COMPRESS_NONE) {
                        size_t raw_sig_depress_bytes;
                        int16_t *raw_sig_depress = slow5_ptr_depress_solo(signal_method, read->raw_signal, read->len_raw_signal, &raw_sig_depress_bytes);
                        if (!raw_sig_depress) {
                            SLOW5_ERROR("%s", "Decompressing raw signal failed.");
                            ret = -1;
                            break;
                        }

                        free(read->raw_signal);
                        read->raw_signal = raw_sig_depress;
                        read->len_raw_signal = raw_sig_depress_bytes / sizeof *read->raw_signal;
                    }

                 } break;

                default: /* All columns parsed */
                    main_cols_parsed = true;
                    break;
            }

            ++ i;
        }

        if (offset > read_size) { /* Read too much */
            SLOW5_ERROR("Read more bytes than record size. At offset '%" PRIu64 "' but record size is '%zu'.", offset, read_size);
            ret = -1;
        } else if (aux_meta == NULL && offset < read_size) { /* More to read but no auxiliary meta */
            SLOW5_ERROR("No auxiliary meta data but more data in record. At offset '%" PRIu64 "' but record size is '%zu'.", offset, read_size);
            ret = -1;
        } else if (aux_meta != NULL && offset == read_size) { /* No more to read but auxiliary meta expected */
            SLOW5_ERROR("Missing auxiliary data in record. At offset '%" PRIu64 "' which equals the record size.", offset);
            ret = -1;
        } else if (aux_meta != NULL) {
            /* Parse auxiliary data */
            ret = slow5_rec_aux_parse(NULL, read_mem, offset, read_size, read, format, aux_meta);
        }
    }

    return ret;
}

/*
 * parse auxiliary data from read_mem with read_size bytes, the given format and auxiliary meta data aux_meta
 * setup the read's auxiliary map
 * read_mem, read, aux_meta cannot be NULL
 * returns -1 on error, 0 on success
 */
static int slow5_rec_aux_parse(char *tok, char *read_mem, uint64_t offset, size_t read_size, struct slow5_rec *read, enum slow5_fmt format, struct slow5_aux_meta *aux_meta) {

    int ret = 0;

    khash_t(slow5_s2a) *aux_map = slow5_rec_aux_init();
    if (aux_map == NULL) {
        return -1;
    }

    if (format == SLOW5_FORMAT_ASCII) {

        for (int64_t i = 0; i < aux_meta->num; ++ i) {
            if (tok == NULL) {
                SLOW5_ERROR("Auxiliary field '%s' is missing in record.", aux_meta->attrs[i]);
                slow5_rec_aux_free(aux_map);
                return -1;
            }

            uint64_t bytes = 0;
            uint64_t len = 1;
            uint8_t *data = NULL;
            if (SLOW5_IS_PTR(aux_meta->types[i])) {
                /* Type is an array */
                if (strcmp(tok, SLOW5_ASCII_MISSING) == 0) {
                    len = 0;
                } else if (aux_meta->types[i] == SLOW5_STRING) {
                    len = strlen(tok);
                    data = (uint8_t *) malloc((len + 1) * aux_meta->sizes[i]);
                    SLOW5_MALLOC_CHK(data);
                    if (data == NULL) {
                        slow5_rec_aux_free(aux_map);
                        return -1;
                    }
                    memcpy(data, tok, (len + 1) * aux_meta->sizes[i]);
                } else {
                    /* Split tok by SLOW5_SEP_ARRAY */
                    char *tok_sep;
                    uint64_t array_cap = SLOW5_AUX_ARRAY_CAP_INIT;
                    uint64_t array_i = 0;
                    data = (uint8_t *) malloc(array_cap * aux_meta->sizes[i]);
                    SLOW5_MALLOC_CHK(data);
                    if (data == NULL) {
                        slow5_rec_aux_free(aux_map);
                        return -1;
                    }

                    while ((tok_sep = slow5_strsep(&tok, SLOW5_SEP_ARRAY)) != NULL) {
                        /*
                         * Storing comma-separated array
                         * Dynamic array creation
                         */

                        /* Memcpy giving the primitive type not the array type */
                        if (slow5_memcpy_type_from_str(data + (aux_meta->sizes[i] * array_i), tok_sep, SLOW5_TO_PRIM_TYPE(aux_meta->types[i])) == -1) {
                            SLOW5_ERROR("Array element '%" PRIu64 "' of auxiliary field '%s' could not be parsed from '%s'.", array_i, aux_meta->attrs[i], tok);
                            free(data);
                            slow5_rec_aux_free(aux_map);
                            return -1;
                        }

                        ++ array_i;

                        if (array_i >= array_cap) {
                            array_cap = array_cap << 1;
                            uint8_t *data_tmp = (uint8_t *) realloc(data, array_cap * aux_meta->sizes[i]);
                            SLOW5_MALLOC_CHK(data_tmp);
                            if (data_tmp == NULL) {
                                free(data);
                                slow5_rec_aux_free(aux_map);
                                return -1;
                            }
                            data = data_tmp;
                        }
                    }

                    len = array_i;
                }
                bytes = len * aux_meta->sizes[i];

            } else { /* primitive type */
                data = (uint8_t *) malloc(aux_meta->sizes[i]);
                SLOW5_MALLOC_CHK(data);
                if (data == NULL) {
                    slow5_rec_aux_free(aux_map);
                    return -1;
                }

                if (strcmp(tok, SLOW5_ASCII_MISSING) == 0) {
                    slow5_memcpy_null_type(data, aux_meta->types[i]);
                } else if (slow5_memcpy_type_from_str(data, tok, aux_meta->types[i]) == -1) {
                    SLOW5_ERROR("Auxiliary field '%s' could not be parsed from '%s'.", aux_meta->attrs[i], tok);
                    free(data);
                    slow5_rec_aux_free(aux_map);
                    return -1;
                }
                bytes = aux_meta->sizes[i];
            }

            int absent;
            khint_t pos = kh_put(slow5_s2a, aux_map, aux_meta->attrs[i], &absent);
            if (absent == -1) {
                SLOW5_ERROR("Inserting key '%s' into the auxiliary hash table failed.", aux_meta->attrs[i]);
                free(data);
                slow5_rec_aux_free(aux_map);
                return -1;
            } else if (absent == 0) {
                SLOW5_ERROR("Auxiliary field '%s' is duplicated.", aux_meta->attrs[i]);
                free(data);
                slow5_rec_aux_free(aux_map);
                return -1;
            }
            struct slow5_rec_aux_data *aux_data = &kh_value(aux_map, pos);
            aux_data->len = len;
            aux_data->bytes = bytes;
            aux_data->data = data;
            aux_data->type = aux_meta->types[i];

            tok = slow5_strsep(&read_mem, SLOW5_SEP_COL);
        }
        /* Ensure line ends */
        if (tok != NULL) {
            SLOW5_ERROR("%s", "More auxiliary record data exists but all header auxiliary columns were parsed.");
            slow5_rec_aux_free(aux_map);
            return -1;
        }
        read->aux_map = aux_map;

    } else if (format == SLOW5_FORMAT_BINARY) {

        size_t size = 0;

        for (int64_t i = 0; i < aux_meta->num; ++ i) {
            if (offset >= read_size) {
                SLOW5_ERROR("Auxiliary field '%s' is missing in record. At offset '%" PRIu64 "' and record size is '%" PRIu64 "'.", aux_meta->attrs[i], offset, read_size);
                slow5_rec_aux_free(aux_map);
                return -1;
            }

            uint64_t len = 1;
            if (SLOW5_IS_PTR(aux_meta->types[i])) {
                /* Type is an array */
                size = sizeof len;
                memcpy(&len, read_mem + offset, size);
                offset += size;
            }

            uint8_t *data;
            uint64_t bytes = len * aux_meta->sizes[i];

            if (len == 0) {
                data = NULL;
            } else if (aux_meta->types[i] == SLOW5_STRING) {
                data = (uint8_t *) malloc(bytes + 1);
                SLOW5_MALLOC_CHK(data);
                if (data == NULL) {
                    slow5_rec_aux_free(aux_map);
                    return -1;
                }
                memcpy(data, read_mem + offset, bytes);
                offset += bytes;
                data[bytes] = '\0';
            } else {
                data = (uint8_t *) malloc(bytes);
                SLOW5_MALLOC_CHK(data);
                if (data == NULL) {
                    slow5_rec_aux_free(aux_map);
                    return -1;
                }
                memcpy(data, read_mem + offset, bytes);
                offset += bytes;
            }

            int absent;
            khint_t pos = kh_put(slow5_s2a, aux_map, aux_meta->attrs[i], &absent);
            if (absent == -1) {
                SLOW5_ERROR("Inserting key '%s' into the auxiliary hash table failed.", aux_meta->attrs[i]);
                slow5_rec_aux_free(aux_map);
                return -1;
            } else if (absent == 0) {
                SLOW5_ERROR("Auxiliary field '%s' is duplicated.", aux_meta->attrs[i]);
                slow5_rec_aux_free(aux_map);
                return -1;
            }
            struct slow5_rec_aux_data *aux_data = &kh_value(aux_map, pos);
            aux_data->len = len;
            aux_data->bytes = bytes;
            aux_data->data = data;
            aux_data->type = aux_meta->types[i];
        }

        if (offset != read_size) {
            SLOW5_ERROR("More auxiliary record data exists but all header auxiliary columns were parsed. At offset '%" PRIu64 "' but record size is '%zu'.", offset, read_size);
            slow5_rec_aux_free(aux_map);
            return -1;
        }
        read->aux_map = aux_map;
    }

    return ret;
}

/* init a slow5 record auxiliary hash map DONE */
static inline khash_t(slow5_s2a) *slow5_rec_aux_init(void) {
    khash_t(slow5_s2a) *aux_map = kh_init(slow5_s2a);
    SLOW5_MALLOC_CHK(aux_map)
    return aux_map;
}

/* free a slow5 record auxiliary hash map DONE */
void slow5_rec_aux_free(khash_t(slow5_s2a) *aux_map) {
    if (aux_map) {
        for (khint_t i = kh_begin(aux_map); i != kh_end(aux_map); ++ i) {
            if (kh_exist(aux_map, i)) {
                kh_del(slow5_s2a, aux_map, i);
                struct slow5_rec_aux_data *aux_data = &kh_value(aux_map, i);
                /* TODO aux_data better not be null here because I'm not checking */
                free(aux_data->data);
            }
        }

        kh_destroy(slow5_s2a, aux_map);
    }
}

/*
 * get next slow5 record from current file pointer as it is stored in s5p with length *n
 * on error returns NULL, sets *n=0 if possible, and sets slow5_errno
 * slow5_errno errors:
 * SLOW5_ERR_ARG
 * SLOW5_ERR_IO
 * SLOW5_ERR_EOF    end of file reached (and blow5 eof marker found)
 * slow5_is_eof errors
 * SLOW5_ERR_TRUNC
 * SLOW5_ERR_MEM
 * SLOW5_ERR_UNK
 */
void *slow5_get_next_mem(size_t *n, const struct slow5_file *s5p) {
    if (!s5p) {
        SLOW5_ERROR("Argument '%s' cannot be NULL.", SLOW5_TO_STR(s5p));
        slow5_errno = SLOW5_ERR_ARG;
        goto err;
    }

    char *mem = NULL;
    size_t bytes;

    if (s5p->format == SLOW5_FORMAT_ASCII) {
        size_t cap = 0;
        ssize_t bytes_tmp = getline(&mem, &cap, s5p->fp);
        if (bytes_tmp == -1) { /* getline failed */
            if (feof(s5p->fp)) {
                free(mem);
                slow5_errno = SLOW5_ERR_EOF;
                goto err;
            }
            SLOW5_ERROR("Reading the next slow5 record until newline failed: %s", strerror(errno));
            free(mem);
            slow5_errno = SLOW5_ERR_IO;
            goto err;
        }
        bytes = bytes_tmp;
        mem[-- bytes] = '\0'; /* remove newline for parsing */

    } else if (s5p->format == SLOW5_FORMAT_BINARY) {
        slow5_rec_size_t bytes_tmp;
        size_t bytes_read;
        /* try to read the record size */
        if ((bytes_read = fread(&bytes_tmp, 1, sizeof bytes_tmp, s5p->fp)) != sizeof bytes_tmp) {
            const char eof[] = SLOW5_BINARY_EOF;
            int is_eof = 0;
            if (bytes_read == sizeof eof) {
                /* check if eof */
                is_eof = slow5_is_eof(s5p->fp, eof, sizeof eof);
                if (is_eof == -1) { /* io/mem error */
                    SLOW5_ERROR("%s", "Internal error while checking for blow5 eof marker.");
                } else if (is_eof == -2) {
                    SLOW5_ERROR("%s", "Malformed blow5. End of file marker found, but end of file not reached.");
                } else if (is_eof == 1) {
                    slow5_errno = SLOW5_ERR_EOF;
                }
            }
            if (!is_eof) { /* catches the case when bytes_read != sizeof eof */
                SLOW5_ERROR("Malformed blow5 record. Failed to read the record size.%s", feof(s5p->fp) ? " Missing blow5 end of file marker." : "");
                if (feof(s5p->fp)) {
                    slow5_errno = SLOW5_ERR_TRUNC;
                } else {
                    slow5_errno = SLOW5_ERR_IO;
                }
            }
            goto err;
        }
        bytes = bytes_tmp;

        mem = (char *) malloc(bytes);
        if (!mem) {
            SLOW5_MALLOC_ERROR();
            slow5_errno = SLOW5_ERR_MEM;
            goto err;
        }

        if (fread(mem, bytes, 1, s5p->fp) != 1) {
            SLOW5_ERROR("Malformed blow5 record. Failed to read '%zu' bytes from blow5 file '%s'.%s",
                    bytes, s5p->meta.pathname, feof(s5p->fp) ? " EOF reached unexpectedly." : "");
            if (feof(s5p->fp)) {
                slow5_errno = SLOW5_ERR_TRUNC;
            } else {
                slow5_errno = SLOW5_ERR_IO;
            }
            free(mem);
            goto err;
        }

    } else {
        SLOW5_ERROR("Invalid slow5 format '%d'.", s5p->format);
        slow5_errno = SLOW5_ERR_UNK;
        goto err;
    }

    if (n) {
        *n = bytes;
    }
    return mem;

    err:
        if (n) {
            *n = 0;
        }
        return NULL;
}


int slow5_get_next_bytes(void **mem, size_t *bytes, slow5_file_t *s5p){
    *mem = slow5_get_next_mem(bytes, s5p);
    if (*mem == NULL) {
        return -1;
    } else {
        return 0;
    }
}

/**
 * Get the read entry under the current file pointer of a slow5 file.
 *
 * Allocates memory for *read if it is NULL.
 * Otherwise, the data in *read is freed and overwritten.
 *
 * slow5_rec_free() should be called when finished with the read, even on failure.
 *
 * Return
 * >=0  the read was successfully found and stored
 * <0   error code described below
 *
 * Not thread safe
 *
 * Error
 * SLOW5_ERR_ARG        read or s5p is NULL
 * SLOW5_ERR_EOF        EOF reached (and eof marker found if blow5)
 * SLOW5_ERR_IO         other reading error when reading the slow5 file
 * SLOW5_ERR_MEM        memory allocation error
 * SLOW5_ERR_PRESS      record decompression error
 * SLOW5_ERR_RECPARSE   record parsing error
 * SLOW5_ERR_TRUNC      record truncated
 * SLOW5_ERR_UNK        slow5 format is unknown
 *
 * @param   read    address of a slow5_rec pointer
 * @param   s5p     slow5 file
 * @return  error code described above
 */
int slow5_get_next(struct slow5_rec **read, struct slow5_file *s5p) {
    if (!read) {
        SLOW5_ERROR_EXIT("Argument '%s' cannot be NULL.", SLOW5_TO_STR(read));
        return slow5_errno = SLOW5_ERR_ARG;
    }

    size_t bytes;
    char *mem;
    if (!(mem = slow5_get_next_mem(&bytes, s5p))) {
        if (slow5_errno != SLOW5_ERR_EOF) {
            SLOW5_EXIT_IF_ON_ERR();
        }
        return slow5_errno;
    }

    if (slow5_rec_depress_parse(&mem, &bytes, NULL, read, s5p) != 0) {
        SLOW5_EXIT_IF_ON_ERR();
        free(mem);
        return slow5_errno;
    }

    free(mem);
    return 0;
}

// For non-array types
// Return
// -1   input invalid
// -2   attr not found
// -3   type is an array type
// -4   data is invalid (eg enum out of range)
// TODO add setting a non-aux?
int slow5_rec_set(struct slow5_rec *read, struct slow5_aux_meta *aux_meta, const char *attr, const void *data) {
    if (read == NULL || aux_meta == NULL || aux_meta->num == 0 || attr == NULL || data == NULL) {
        return -1;
    }

    khint_t pos = kh_get(slow5_s2ui32, aux_meta->attr_to_pos, attr);
    if (pos == kh_end(aux_meta->attr_to_pos)) {
        return -2;
    }

    uint32_t i = kh_value(aux_meta->attr_to_pos, pos);

    if (SLOW5_IS_PTR(aux_meta->types[i])) {
        return -3;
    }

    // TODO warning if setting null attribute here

    uint8_t *data_cast = (uint8_t *) data;
    if (aux_meta->types[i] == SLOW5_ENUM) {
        if (!aux_meta->enum_labels) {
            return -1;
        } else if (*data_cast >= aux_meta->enum_num_labels[i]) {
            /* check if enum number is out of bounds */
            return -4;
        }
    }

    if (read->aux_map == NULL) {
        read->aux_map = kh_init(slow5_s2a);
    }
    slow5_rec_set_aux_map(read->aux_map, attr, data_cast, 1, aux_meta->sizes[i], aux_meta->types[i]);

    return 0;
}


int slow5_aux_set(slow5_rec_t *read, const char *field, const void *data, slow5_hdr_t *header){
    int ret = slow5_rec_set(read, header->aux_meta, field, data);
    return ret;
}


// For array types
// Return
// 0    success
// -1   input invalid
// -2   attr not found
// -3   type is not an array type
// -4   a data value is invalid (eg enum out of range)
int slow5_rec_set_array(struct slow5_rec *read, struct slow5_aux_meta *aux_meta, const char *attr, const void *data, size_t len) {
    if (read == NULL || aux_meta == NULL || aux_meta->num == 0 || attr == NULL || data == NULL) {
        return -1;
    }

    khint_t pos = kh_get(slow5_s2ui32, aux_meta->attr_to_pos, attr);
    if (pos == kh_end(aux_meta->attr_to_pos)) {
        return -2;
    }

    uint32_t i = kh_value(aux_meta->attr_to_pos, pos);

    if (!SLOW5_IS_PTR(aux_meta->types[i])) {
        return -3;
    }

    // TODO warning if setting null attribute here

    uint8_t *data_cast = (uint8_t *) data;
    if (aux_meta->types[i] == SLOW5_ENUM_ARRAY) {
        if (!aux_meta->enum_labels) {
            return -1;
        } else {
            /* check if an enum number is out of bounds */
            for (uint16_t j = 0; j < len; ++ j) {
                if (data_cast[j] >= aux_meta->enum_num_labels[i]) {
                    return -4;
                }
            }
        }
    }

    if (read->aux_map == NULL) {
        read->aux_map = kh_init(slow5_s2a);
    }
    slow5_rec_set_aux_map(read->aux_map, attr, data_cast, len, aux_meta->sizes[i] * len, aux_meta->types[i]);

    return 0;
}

int slow5_aux_set_string(slow5_rec_t *read, const char *field, const char *data, slow5_hdr_t *header) {
    int ret = slow5_rec_set_string(read, header->aux_meta, field, data);
    return ret;
}

static inline void slow5_rec_set_aux_map(khash_t(slow5_s2a) *aux_map, const char *attr, const uint8_t *data, size_t len, uint64_t bytes, enum slow5_aux_type type) {
    khint_t pos = kh_get(slow5_s2a, aux_map, attr);
    struct slow5_rec_aux_data *aux_data;
    if (pos != kh_end(aux_map)) {
        aux_data = &kh_value(aux_map, pos);
    } else {
        int ret;
        pos = kh_put(slow5_s2a, aux_map, attr, &ret);
        if (ret == -1) {
            SLOW5_ERROR("%s", "Internal klib error.");
            return;
        }
        aux_data = &kh_value(aux_map, pos);
    }
    aux_data->len = len;
    aux_data->bytes = bytes;
    aux_data->type = type;
    if (type == SLOW5_STRING) {
        aux_data->data = (uint8_t *) malloc(bytes + 1);
        aux_data->data[bytes] = '\0';
    } else {
        aux_data->data = (uint8_t *) malloc(bytes);
    }
    SLOW5_MALLOC_CHK(aux_data->data);
    memcpy(aux_data->data, data, bytes);
}

#define SLOW5_AUX_GET(read, field, raw_type, null, prim_type) \
    int tmp_err = 0; \
    raw_type val = null; \
    if (!read || !field) { \
        if (!read) { \
            SLOW5_ERROR_EXIT("Argument '%s' cannot be NULL.", SLOW5_TO_STR(read)) \
        } \
        if (!field) { \
            SLOW5_ERROR_EXIT("Argument '%s' cannot be NULL.", SLOW5_TO_STR(field)) \
        } \
        tmp_err = slow5_errno = SLOW5_ERR_ARG; \
    } else if (read->aux_map == NULL) { \
        SLOW5_ERROR_EXIT("%s", "Missing auxiliary hash map.") \
        tmp_err = slow5_errno = SLOW5_ERR_NOAUX; \
    } else { \
        khint_t pos = kh_get(slow5_s2a, read->aux_map, field); \
        if (pos == kh_end(read->aux_map)) { \
            SLOW5_ERROR_EXIT("Field '%s' not found.", field) \
            tmp_err = slow5_errno = SLOW5_ERR_NOFLD; \
        } else  { \
            struct slow5_rec_aux_data aux_data = kh_value(read->aux_map, pos); \
            if (aux_data.type == prim_type) { \
                val = *((raw_type*) aux_data.data); \
            } else { \
                SLOW5_ERROR_EXIT("Desired type '%s' is different to actual type '%s' of field '%s'.", \
                        SLOW5_TO_STR(raw_type), SLOW5_AUX_TYPE_META[prim_type].type_str, field); \
                tmp_err = slow5_errno = SLOW5_ERR_TYPE; \
            } \
        } \
    } \
    if (err != NULL) { /* this refers to err in the prototype slow5_aux_get_*_array */ \
        *err = tmp_err; \
    } \
    return val;
/*
 * get the auxiliary entry of a primitive type field from a slow5 record
 * err is deprecated, check slow5_errno instead
 * SLOW5_ERR_ARG if read or field is NULL
 * SLOW5_ERR_NOAUX if no auxiliary hash map for the record
 * SLOW5_ERR_NOFLD if the field was not found
 * SLOW5_ERR_TYPE if the desired return type does not match the field's type
 * returns NULL equivalent of return type on error
 * DONE
 */
int8_t slow5_aux_get_int8(const struct slow5_rec *read, const char *field, int *err) {
    SLOW5_AUX_GET(read, field, int8_t, SLOW5_INT8_T_NULL, SLOW5_INT8_T)
}
int16_t slow5_aux_get_int16(const struct slow5_rec *read, const char *field, int *err) {
    SLOW5_AUX_GET(read, field, int16_t, SLOW5_INT16_T_NULL, SLOW5_INT16_T)
}
int32_t slow5_aux_get_int32(const struct slow5_rec *read, const char *field, int *err) {
    SLOW5_AUX_GET(read, field, int32_t, SLOW5_INT32_T_NULL, SLOW5_INT32_T)
}
int64_t slow5_aux_get_int64(const struct slow5_rec *read, const char *field, int *err) {
    SLOW5_AUX_GET(read, field, int64_t, SLOW5_INT64_T_NULL, SLOW5_INT64_T)
}
uint8_t slow5_aux_get_uint8(const struct slow5_rec *read, const char *field, int *err) {
    SLOW5_AUX_GET(read, field, uint8_t, SLOW5_UINT8_T_NULL, SLOW5_UINT8_T)
}
uint16_t slow5_aux_get_uint16(const struct slow5_rec *read, const char *field, int *err) {
    SLOW5_AUX_GET(read, field, uint16_t, SLOW5_UINT16_T_NULL, SLOW5_UINT16_T)
}
uint32_t slow5_aux_get_uint32(const struct slow5_rec *read, const char *field, int *err) {
    SLOW5_AUX_GET(read, field, uint32_t, SLOW5_UINT32_T_NULL, SLOW5_UINT32_T)
}
uint64_t slow5_aux_get_uint64(const struct slow5_rec *read, const char *field, int *err) {
    SLOW5_AUX_GET(read, field, uint64_t, SLOW5_UINT64_T_NULL, SLOW5_UINT64_T)
}
float slow5_aux_get_float(const struct slow5_rec *read, const char *field, int *err) {
    SLOW5_AUX_GET(read, field, float, SLOW5_FLOAT_NULL, SLOW5_FLOAT)
}
double slow5_aux_get_double(const struct slow5_rec *read, const char *field, int *err) {
    SLOW5_AUX_GET(read, field, double, SLOW5_DOUBLE_NULL, SLOW5_DOUBLE)
}
char slow5_aux_get_char(const struct slow5_rec *read, const char *field, int *err) {
    SLOW5_AUX_GET(read, field, char, SLOW5_CHAR_NULL, SLOW5_CHAR)
}
uint8_t slow5_aux_get_enum(const struct slow5_rec *read, const char *field, int *err) {
    SLOW5_AUX_GET(read, field, uint8_t, SLOW5_ENUM_NULL, SLOW5_ENUM)
}

#define SLOW5_AUX_GET_ARRAY(read, field, len, raw_type, ptr_type) \
    int tmp_err = 0; \
    raw_type val = NULL; \
    if (!read || !field) { \
        if (!read) { \
            SLOW5_ERROR_EXIT("Argument '%s' cannot be NULL.", SLOW5_TO_STR(read)) \
        } \
        if (!field) { \
            SLOW5_ERROR_EXIT("Argument '%s' cannot be NULL.", SLOW5_TO_STR(field)) \
        } \
        tmp_err = slow5_errno = SLOW5_ERR_ARG; \
    } else if (read->aux_map == NULL) { \
        SLOW5_ERROR_EXIT("%s", "Missing auxiliary hash map.") \
        tmp_err = slow5_errno = SLOW5_ERR_NOAUX; \
    } else { \
        khint_t pos = kh_get(slow5_s2a, read->aux_map, field); \
        if (pos == kh_end(read->aux_map)) { \
            SLOW5_ERROR_EXIT("Field '%s' not found.", field) \
            tmp_err = slow5_errno = SLOW5_ERR_NOFLD; \
        } else { \
            struct slow5_rec_aux_data aux_data = kh_value(read->aux_map, pos); \
            if (aux_data.type == ptr_type) { \
                val = (raw_type) aux_data.data; \
                if (len != NULL) { \
                    *len = aux_data.len; \
                } \
            } else { \
                SLOW5_ERROR_EXIT("Desired type '%s' is different to actual type '%s' of field '%s'.", \
                        SLOW5_TO_STR(raw_type), SLOW5_AUX_TYPE_META[SLOW5_TO_PRIM_TYPE(ptr_type)].type_str, field); \
                tmp_err = slow5_errno = SLOW5_ERR_TYPE; \
            } \
        } \
    } \
    if (err != NULL) { /* this refers to err in the prototype slow5_aux_get_*_array */ \
        *err = tmp_err; \
    } \
    return val;
/*
 * get the auxiliary entry of an array type field from a slow5 record
 * err is deprecated, check slow5_errno instead
 * SLOW5_ERR_ARG if read or field is NULL
 * SLOW5_ERR_NOAUX if no auxiliary hash map for the record
 * SLOW5_ERR_NOFLD if the field was not found
 * SLOW5_ERR_TYPE if the desired return type does not match the field's type
 * return NULL on error
 * DONE
 */
int8_t *slow5_aux_get_int8_array(const struct slow5_rec *read, const char *field, uint64_t *len, int *err) {
    SLOW5_AUX_GET_ARRAY(read, field, len, int8_t*, SLOW5_INT8_T_ARRAY)
}
int16_t *slow5_aux_get_int16_array(const struct slow5_rec *read, const char *field, uint64_t *len, int *err) {
    SLOW5_AUX_GET_ARRAY(read, field, len, int16_t*, SLOW5_INT16_T_ARRAY)
}
int32_t *slow5_aux_get_int32_array(const struct slow5_rec *read, const char *field, uint64_t *len, int *err) {
    SLOW5_AUX_GET_ARRAY(read, field, len, int32_t*, SLOW5_INT32_T_ARRAY)
}
int64_t *slow5_aux_get_int64_array(const struct slow5_rec *read, const char *field, uint64_t *len, int *err) {
    SLOW5_AUX_GET_ARRAY(read, field, len, int64_t*, SLOW5_INT64_T_ARRAY)
}
uint8_t *slow5_aux_get_uint8_array(const struct slow5_rec *read, const char *field, uint64_t *len, int *err) {
    SLOW5_AUX_GET_ARRAY(read, field, len, uint8_t*, SLOW5_UINT8_T_ARRAY)
}
uint16_t *slow5_aux_get_uint16_array(const struct slow5_rec *read, const char *field, uint64_t *len, int *err) {
    SLOW5_AUX_GET_ARRAY(read, field, len, uint16_t*, SLOW5_UINT16_T_ARRAY)
}
uint32_t *slow5_aux_get_uint32_array(const struct slow5_rec *read, const char *field, uint64_t *len, int *err) {
    SLOW5_AUX_GET_ARRAY(read, field, len, uint32_t*, SLOW5_UINT32_T_ARRAY)
}
uint64_t *slow5_aux_get_uint64_array(const struct slow5_rec *read, const char *field, uint64_t *len, int *err) {
    SLOW5_AUX_GET_ARRAY(read, field, len, uint64_t*, SLOW5_UINT64_T_ARRAY)
}
float *slow5_aux_get_float_array(const struct slow5_rec *read, const char *field, uint64_t *len, int *err) {
    SLOW5_AUX_GET_ARRAY(read, field, len, float*, SLOW5_FLOAT_ARRAY)
}
double *slow5_aux_get_double_array(const struct slow5_rec *read, const char *field, uint64_t *len, int *err) {
    SLOW5_AUX_GET_ARRAY(read, field, len, double*, SLOW5_DOUBLE_ARRAY)
}
char *slow5_aux_get_string(const struct slow5_rec *read, const char *field, uint64_t *len, int *err) {
    SLOW5_AUX_GET_ARRAY(read, field, len, char*, SLOW5_STRING)
}
uint8_t *slow5_aux_get_enum_array(const struct slow5_rec *read, const char *field, uint64_t *len, int *err) {
    SLOW5_AUX_GET_ARRAY(read, field, len, uint8_t*, SLOW5_ENUM_ARRAY)
}

/**
 * Add a read entry to the slow5 file.
 *
 * Return
 *  0   the read was successfully stored
 * -1   read or s5p or read->read_id is NULL
 * -2   the index was not previously init and failed to init
 * -3   duplicate read id
 * -4   writing failure
 *
 * @param   read    slow5_rec ptr
 * @param   s5p     slow5 file
 * @return  error code described above
 */
int slow5_add_rec(struct slow5_rec *read, struct slow5_file *s5p) {
    if (read == NULL || read->read_id == NULL || s5p == NULL) {
        return -1;
    }

    // Create index if NULL
    if (s5p->index == NULL && (s5p->index = slow5_idx_init(s5p)) == NULL) {
        // index failed to init
        return -2;
    }

    // Duplicate read id
    if (slow5_idx_get(s5p->index, read->read_id, NULL) == 0) {
        return -3;
    }

    // Append record to file
    void *mem = NULL;
    size_t bytes;
    if ((mem = slow5_rec_to_mem(read, s5p->header->aux_meta, s5p->format, s5p->compress, &bytes)) == NULL) {
        return -4;
    }
    if (fseek(s5p->fp, 0L, SEEK_END) != 0) {
        free(mem);
        return -4;
    }
    uint64_t offset = ftello(s5p->fp);
    if (fwrite(mem, bytes, 1, s5p->fp) != 1) {
        free(mem);
        return -4;
    }
    free(mem);

    // Update index
    slow5_idx_insert(s5p->index, strdup(read->read_id), offset, bytes);

    //after updating mark dirty
    s5p->index->dirty = 1;

    return 0;
}

/**
 * Remove a read entry at a read_id in a slow5 file.
 *
 * Return
 *  0   the read was successfully stored
 * -1   an input parameter is NULL
 * -2   the index was not previously init and failed to init
 * -3   read_id was not found in the index
 *
 * @param   read_id the read identifier
 * @param   s5p     slow5 file
 * @return  error code described above
 */
int slow5_rec_rm(const char *read_id, struct slow5_file *s5p) {
    if (read_id == NULL || s5p == NULL) {
        return -1;
    }

    // Create index if NULL
    if (s5p->index == NULL && (s5p->index = slow5_idx_init(s5p)) == NULL) {
        // index failed to init
        return -2;
    }

    // Get index record
    struct slow5_rec_idx read_index;
    if (slow5_idx_get(s5p->index, read_id, &read_index) == -1) {
        // read_id not found in index
        return -3;
    }

    // TODO
    // remove record from file
    // update index

    return 0;
}

/**
 * Print a read entry in the correct format with newline character to a file pointer.
 *
 * Error if fp or read is NULL,
 * of if the format is SLOW5_FORMAT_UNKNOWN.
 *
 * On success, the number of bytes written is returned.
 * On error, -1 is returned.
 *
 * @param   fp      output file pointer
 * @param   read    slow5_rec pointer
 * @return  number of bytes written, -1 on error
 */
int slow5_rec_fwrite(FILE *fp, struct slow5_rec *read, struct slow5_aux_meta *aux_meta, enum slow5_fmt format, struct slow5_press *compress) {
    int ret;
    void *read_mem;
    size_t read_size;

    if (fp == NULL || read == NULL || (read_mem = slow5_rec_to_mem(read, aux_meta, format, compress, &read_size)) == NULL) {
        return -1;
    }

    size_t n = fwrite(read_mem, read_size, 1, fp);
    if (n != 1) {
        ret = -1;
    } else {
        ret = read_size; // TODO is this okay
    }

    free(read_mem);
    return ret;
}

int slow5_write_bytes(void *mem, size_t bytes, slow5_file_t *s5p){
    size_t n = fwrite(mem, bytes, 1, s5p->fp);
    int ret;
    if (n != 1) {
        ret = -1;
    } else {
        ret = 0; // TODO is this okay
    }
    return ret;
}

int slow5_write(slow5_rec_t *rec, slow5_file_t *s5p){
    int ret = slow5_rec_fwrite(s5p->fp, rec, s5p->header->aux_meta, s5p->format, s5p->compress);
    return ret;
}


/**
 * Get the read entry in the specified format.
 *
 * Returns NULL if read is NULL,
 * or format is SLOW5_FORMAT_UNKNOWN,
 * or the read attribute values are invalid
 *
 * @param   read        slow5_rec pointer
 * @param   format      slow5 format to write the entry in
 * @param   compress    compress structure
 * @param   n           number of bytes written to the returned buffer
 * @return  malloced string to use free() on, NULL on error
 */
void *slow5_rec_to_mem(struct slow5_rec *read, struct slow5_aux_meta *aux_meta, enum slow5_fmt format, struct slow5_press *compress, size_t *n) {
    char *mem = NULL;

    if (read == NULL || format == SLOW5_FORMAT_UNKNOWN) {
        return NULL;
    }

    size_t curr_len = 0;

    if (format == SLOW5_FORMAT_ASCII) {

        char *digitisation_str = slow5_double_to_str(read->digitisation, NULL);
        char *offset_str = slow5_double_to_str(read->offset, NULL);
        char *range_str = slow5_double_to_str(read->range, NULL);
        char *sampling_rate_str = slow5_double_to_str(read->sampling_rate, NULL);

        // Set read id to "" if NULL
        const char *read_id = read->read_id;
        if (read->read_id == NULL) {
            read_id = "";
        }

        int curr_len_tmp = slow5_asprintf(&mem,
                SLOW5_COLS(SLOW5_GENERATE_FORMAT_STRING_SEP, SLOW5_GENERATE_NULL),
                read_id,
                read->read_group,
                digitisation_str,
                offset_str,
                range_str,
                sampling_rate_str,
                read->len_raw_signal);
        free(digitisation_str);
        free(offset_str);
        free(range_str);
        free(sampling_rate_str);
        if (curr_len_tmp > 0) {
            curr_len = curr_len_tmp;
        } else {
            free(mem);
            return NULL;
        }

        // TODO memory optimise
        // <max length> = <current length> + (<max signal length> + ','/'\n') * <number of signals> + '\0'
        // <max length> = <current length> + '\n' + '\0'
        const size_t max_len = read->len_raw_signal != 0 ? curr_len + (INT16_MAX_LENGTH + 1) * read->len_raw_signal + 1 : curr_len + 1 + 1;
        mem = (char *) realloc(mem, max_len * sizeof *mem);
        SLOW5_MALLOC_CHK(mem);

        char sig_buf[SLOW5_SIGNAL_BUF_FIXED_CAP];

        uint64_t i;
        for (i = 1; i < read->len_raw_signal; ++ i) {
            int sig_len = sprintf(sig_buf, SLOW5_FORMAT_STRING_RAW_SIGNAL SLOW5_SEP_ARRAY, read->raw_signal[i-1]);

            memcpy(mem + curr_len, sig_buf, sig_len);
            curr_len += sig_len;
        }
        if (read->len_raw_signal > 0) {
            // Trailing signal
            int len_to_cp = sprintf(sig_buf, SLOW5_FORMAT_STRING_RAW_SIGNAL, read->raw_signal[i-1]);
            memcpy(mem + curr_len, sig_buf, len_to_cp);
            curr_len += len_to_cp;
        }

        // Auxiliary fields
        size_t cap = max_len;
        if (read->aux_map != NULL && aux_meta != NULL) {
            for (uint64_t i = 0; i < aux_meta->num; ++ i) {

                // Realloc if necessary
                if (curr_len + 1 >= cap) { // +1 for '\t'
                    cap *= 2;
                    mem = (char *) realloc(mem, cap);
                    SLOW5_MALLOC_CHK(mem);
                }
                mem[curr_len ++] = '\t';

                size_t type_len;
                char *type_str;

                khint_t pos = kh_get(slow5_s2a, read->aux_map, aux_meta->attrs[i]);
                if (pos != kh_end(read->aux_map)) {
                    struct slow5_rec_aux_data aux_data = kh_value(read->aux_map, pos);
                    type_str = slow5_data_to_str(aux_data.data, aux_data.type, aux_data.len, &type_len);
                } else {
                    type_str = get_missing_str(&type_len);
                }

                // Realloc if necessary
                if (curr_len + type_len >= cap) {
                    cap *= 2;
                    mem = (char *) realloc(mem, cap);
                    SLOW5_MALLOC_CHK(mem);
                }

                memcpy(mem + curr_len, type_str, type_len);
                curr_len += type_len;

                free(type_str);
            }
        }

        // Trailing newline
        // Realloc if necessary
        if (curr_len + 2 >= cap) { // +2 for '\n' and '\0'
            cap *= 2;
            mem = (char *) realloc(mem, cap);
            SLOW5_MALLOC_CHK(mem);
        }
        strcpy(mem + curr_len, "\n"); // Copies null byte as well
        curr_len += 1;

    } else if (format == SLOW5_FORMAT_BINARY) {

        bool compress_to_free = false;
        if (compress == NULL) {
            slow5_press_method_t method = {SLOW5_COMPRESS_NONE,SLOW5_COMPRESS_NONE};
            compress = slow5_press_init(method);
            /* TODO error if this fails */
            // I added a quick fix below to prevent a dangerous situation, but error codes and cleaning up must be propoerly done - hasindu
            if(compress==NULL){
                return NULL;
            }

            compress_to_free = true;
        }

        size_t cap = sizeof read->read_id_len +
                     read->read_id_len * sizeof *read->read_id +
                     sizeof read->read_group +
                     sizeof read->digitisation +
                     sizeof read->offset +
                     sizeof read->range +
                     sizeof read->sampling_rate +
                     sizeof read->len_raw_signal +
                     read->len_raw_signal * sizeof read->raw_signal;
        mem = (char *) malloc(cap * sizeof *mem);
        SLOW5_MALLOC_CHK(mem);

        memcpy(mem + curr_len, &read->read_id_len, sizeof read->read_id_len);
        curr_len += sizeof read->read_id_len;
        memcpy(mem + curr_len, read->read_id, read->read_id_len * sizeof *read->read_id);
        curr_len += read->read_id_len * sizeof *read->read_id;
        memcpy(mem + curr_len, &read->read_group, sizeof read->read_group);
        curr_len += sizeof read->read_group;
        memcpy(mem + curr_len, &read->digitisation, sizeof read->digitisation);
        curr_len += sizeof read->digitisation;
        memcpy(mem + curr_len, &read->offset, sizeof read->offset);
        curr_len += sizeof read->offset;
        memcpy(mem + curr_len, &read->range, sizeof read->range);
        curr_len += sizeof read->range;
        memcpy(mem + curr_len, &read->sampling_rate, sizeof read->sampling_rate);
        curr_len += sizeof read->sampling_rate;

        size_t bytes_raw_sig = read->len_raw_signal * sizeof *read->raw_signal;

        if (compress->signal_press->method != SLOW5_COMPRESS_NONE) { /* signal compression */
            uint8_t *raw_sig_svb = slow5_ptr_compress_solo(compress->signal_press->method, read->raw_signal, read->len_raw_signal * sizeof *read->raw_signal, &bytes_raw_sig);
            if (!raw_sig_svb) {
                SLOW5_ERROR("%s", "Signal compression failed.");
                free(mem);
                if (compress_to_free) {
                    slow5_press_free(compress);
                }
                return NULL;
            }
            free(read->raw_signal);
            read->len_raw_signal = bytes_raw_sig;
            read->raw_signal = (int16_t *) raw_sig_svb;
        }

        memcpy(mem + curr_len, &read->len_raw_signal, sizeof read->len_raw_signal);
        curr_len += sizeof read->len_raw_signal;
        memcpy(mem + curr_len, read->raw_signal, bytes_raw_sig);
        curr_len += bytes_raw_sig;

        // Auxiliary fields
        if (read->aux_map != NULL && aux_meta != NULL) { // TODO error if one is NULL but not another
            for (uint64_t i = 0; i < aux_meta->num; ++ i) {
                struct slow5_rec_aux_data aux_data = { 0 };

                bool malloc_flag = false;

                khint_t pos = kh_get(slow5_s2a, read->aux_map, aux_meta->attrs[i]);
                if (pos != kh_end(read->aux_map)) {
                    aux_data = kh_value(read->aux_map, pos);
                } else if (!SLOW5_IS_PTR(aux_meta->types[i])) {
                    aux_data.len = 1;
                    aux_data.bytes = SLOW5_AUX_TYPE_META[aux_meta->types[i]].size;
                    aux_data.data = (uint8_t *) malloc(aux_data.bytes);
                    malloc_flag = true;
                    slow5_memcpy_null_type(aux_data.data, aux_meta->types[i]);
                }

                if (SLOW5_IS_PTR(aux_meta->types[i])) {
                    // Realloc if necessary
                    if (curr_len + sizeof aux_data.len + aux_data.bytes >= cap) {
                        cap *= 2;
                        mem = (char *) realloc(mem, cap);
                        SLOW5_MALLOC_CHK(mem);
                    }

                    memcpy(mem + curr_len, &aux_data.len, sizeof aux_data.len);
                    curr_len += sizeof aux_data.len;

                } else if (curr_len + aux_data.bytes >= cap) { // Realloc if necessary
                    cap *= 2;
                    mem = (char *) realloc(mem, cap);
                    SLOW5_MALLOC_CHK(mem);
                }

                if (aux_data.len != 0) {
                    memcpy(mem + curr_len, aux_data.data, aux_data.bytes);
                }
                curr_len += aux_data.bytes;

                if (malloc_flag) {
                    free(aux_data.data);
                }
            }
        }

        slow5_compress_footer_next(compress->record_press);
        slow5_rec_size_t record_size;

        size_t record_sizet;
        void *comp_mem = slow5_ptr_compress(compress->record_press, mem, curr_len, &record_sizet);
        record_size = record_sizet;
        free(mem);

        if (comp_mem != NULL) {
            uint8_t *comp_mem_full = (uint8_t *) malloc(sizeof record_size + record_size);
            SLOW5_MALLOC_CHK(comp_mem_full);
            // Copy size of compressed record
            memcpy(comp_mem_full, &record_size, sizeof record_size);
            // Copy compressed record
            memcpy(comp_mem_full + sizeof record_size, comp_mem, record_size);
            free(comp_mem);

            curr_len = sizeof record_size + record_size;
            mem = (char *) comp_mem_full;
        } else {
            free(mem);
            curr_len = 0;
            mem = NULL;
        }

        if (compress_to_free) {
            slow5_press_free(compress);
        }
    }

    if (n != NULL) {
        *n = curr_len;
    }

    return (void *) mem;
}

int slow5_encode(void **mem, size_t *bytes, slow5_rec_t *read, slow5_file_t *s5p){

    slow5_press_t *press_ptr = NULL;

    if(s5p->compress){

        //todo - check error
        //assert(sf->compress->record_press!=NULL);
        //assert(sf->compress->signal_press!=NULL);

        slow5_press_method_t press_out = {s5p->compress->record_press->method, s5p->compress->signal_press->method};
        press_ptr = slow5_press_init(press_out);
        if(!press_ptr){
            SLOW5_ERROR("Could not initialize the slow5 compression method%s","");
            return -1;
        }
    }

    *mem = slow5_rec_to_mem(read, s5p->header->aux_meta, s5p->format, press_ptr, bytes);
    slow5_press_free(press_ptr);

    if(*mem == NULL){
        SLOW5_ERROR("Could not encode the slow5 record%s","");
        return -1;
    }

    return 0;

}

/*
 * free a slow5 record
 * if slow5_rec_free(read) has already been called undefined behaviour occurs
 * DONE
 */
void slow5_rec_free(struct slow5_rec *read) {
    if (read != NULL) {
        free(read->read_id);
        free(read->raw_signal);
        slow5_rec_aux_free(read->aux_map);
        free(read);
    }
}


/**
 * Create the index file for slow5 file.
 * Overwrites if already exists.
 *
 * <0 on error,
 * >=0 on success.
 *
 * @param   s5p slow5 file structure
 * @return  error codes described above
 */
int slow5_idx_create(struct slow5_file *s5p) {
    char *index_pathname;
    if (s5p == NULL || s5p->meta.pathname == NULL ||
        (index_pathname = slow5_get_idx_path(s5p->meta.pathname)) == NULL) {
        return -1;
    } else if (slow5_idx_to(s5p, index_pathname) == -1) {
        free(index_pathname);
        return -1;
    }

    free(index_pathname);
    return 0;
}

/**
 * Loads the index file for slow5 file.
 * Creates the index if not found.
 *
 * <0 on error,
 * >=0 on success.
 *
 * @param   s5p slow5 file structure
 * @return  error codes described above
 */
int slow5_idx_load(struct slow5_file *s5p) {
    s5p->index = slow5_idx_init(s5p);
    if (s5p->index) {
        return 0;
    } else {
        return -1;
    }
}

void slow5_idx_unload(struct slow5_file *s5p) {
    slow5_idx_free(s5p->index);
    s5p->index = NULL;
    return;
}

/**
 * Print the binary end of file to a file pointer.
 *
 * On success, the number of bytes written is returned.
 * On error, slow5_errno is set to SLOW5_ERR_IO and returned. Errno should be checked for details.
 *
 * @param   fp      output file pointer
 * @return  number of bytes written, SLOW5_ERR_IO on error
 */
ssize_t slow5_eof_fwrite(FILE *fp) {
    const char eof[] = SLOW5_BINARY_EOF;

    size_t n;
    if ((n = fwrite(eof, sizeof *eof, sizeof eof, fp)) != sizeof eof) {
        SLOW5_ERROR("%s", "Could not write blow5 end of file.")
        return slow5_errno = SLOW5_ERR_IO;
    } else {
        return n;
    }
}


/*
 * get the slow5 format from the string name
 * returns SLOW5_FORMAT_UNKNOWN if name is NULL or has no match
 */
enum slow5_fmt slow5_name_get_fmt(const char *name) {
    enum slow5_fmt format = SLOW5_FORMAT_UNKNOWN;

    if (name != NULL) {
        for (size_t i = 0; i < SLOW5_LENGTH(SLOW5_FORMAT_META); ++ i) {
            const struct slow5_fmt_meta meta = SLOW5_FORMAT_META[i];
            if (strcmp(meta.name, name) == 0) {
                format = meta.format;
                break;
            }
        }
    }

    return format;
}

/*
 * get the slow5 format by parsing the file path extension
 * returns SLOW5_FORMAT_UNKNOWN if path is NULL or not found
 */
enum slow5_fmt slow5_path_get_fmt(const char *path) {
    enum slow5_fmt format = SLOW5_FORMAT_UNKNOWN;

    if (path) {
        const char *dot = strrchr(path, '.');
        if (dot && path[strlen(path) - 1] != '.') {
            format = slow5_name_get_fmt(dot + 1);
        }
    }

    return format;
}

/*
 * get the unmalloced string representation of the slow5 format
 * do not try to modify the returned string
 */
const char *slow5_fmt_get_name(enum slow5_fmt format) {
    const char *str = NULL;

    for (size_t i = 0; i < SLOW5_LENGTH(SLOW5_FORMAT_META); ++ i) {
        const struct slow5_fmt_meta meta = SLOW5_FORMAT_META[i];
        if (meta.format == format) {
            str = meta.name;
            break;
        }
    }

    return str;
}

/*
 * get the index file path from the slow5 file path
 * path cannot be NULL
 */
char *slow5_get_idx_path(const char *path) {
    size_t new_len = strlen(path) + strlen(SLOW5_INDEX_EXTENSION);
    char *str = (char *) malloc((new_len + 1) * sizeof *str); // +1 for '\0'
    SLOW5_MALLOC_CHK(str);
    memcpy(str, path, strlen(path));
    strcpy(str + strlen(path), SLOW5_INDEX_EXTENSION);

    return str;
}


// Return
// 0    success
// -1   input invalid
// -2   failure
int slow5_convert(struct slow5_file *from, FILE *to_fp, enum slow5_fmt to_format, slow5_press_method_t to_compress) {
    if (from == NULL || to_fp == NULL || to_format == SLOW5_FORMAT_UNKNOWN) {
        return -1;
    }

    if (slow5_hdr_fwrite(to_fp, from->header, to_format, to_compress) == -1) {
        return -2;
    }

    struct slow5_rec *read = NULL;
    int ret;
    struct slow5_press *press_ptr = slow5_press_init(to_compress);
    if (press_ptr == NULL) {
        return -2;
    }
    while ((ret = slow5_get_next(&read, from)) == 0) {
        if (slow5_rec_fwrite(to_fp, read, from->header->aux_meta, to_format, press_ptr) == -1) {
            slow5_press_free(press_ptr);
            slow5_rec_free(read);
            return -2;
        }
    }
    slow5_press_free(press_ptr);
    slow5_rec_free(read);
    if (ret != SLOW5_ERR_EOF) {
        return -2;
    }

    if (to_format == SLOW5_FORMAT_BINARY) {
        if (slow5_eof_fwrite(to_fp) == -1) {
            return -2;
        }
    }

    return 0;
}

#define STR_TO_AUX_TYPE(str, err, type, len, raw_type, prim_type, ptr_type) \
    (SLOW5_IS_TYPE_TRUNC(str, raw_type)) { \
        if (len == sizeof (# raw_type) - 1) { \
            type = prim_type; \
            *err = 0; \
        } else if (len == sizeof (# raw_type) && str[len - 1] == '*') { \
            type = ptr_type; \
            *err = 0; \
        } else if (prim_type == SLOW5_ENUM) { /* allow extra characters for enum type */ \
            if (str[sizeof (# raw_type) - 1] != '*') { \
                type = prim_type; \
            } else { \
                type = ptr_type; \
            } \
            *err = 0; \
        } else { \
            *err = -1; \
        } \
    }
/*
 * get slow5 type from english string representation
 * str, err cannot be NULL
 * *err set to 0 on success, -1 on failure
 * TODO could return an error type -1
 */
enum slow5_aux_type slow5_str_to_aux_type(const char *str, int *err) {
    enum slow5_aux_type type = SLOW5_INT8_T;

    size_t len = strlen(str);

    if STR_TO_AUX_TYPE(str, err, type, len,         int8_t,     SLOW5_INT8_T,   SLOW5_INT8_T_ARRAY)
    else if STR_TO_AUX_TYPE(str, err, type, len,    int16_t,    SLOW5_INT16_T,  SLOW5_INT16_T_ARRAY)
    else if STR_TO_AUX_TYPE(str, err, type, len,    int32_t,    SLOW5_INT32_T,  SLOW5_INT32_T_ARRAY)
    else if STR_TO_AUX_TYPE(str, err, type, len,    int64_t,    SLOW5_INT64_T,  SLOW5_INT64_T_ARRAY)
    else if STR_TO_AUX_TYPE(str, err, type, len,    uint8_t,    SLOW5_UINT8_T,  SLOW5_UINT8_T_ARRAY)
    else if STR_TO_AUX_TYPE(str, err, type, len,    uint16_t,   SLOW5_UINT16_T, SLOW5_UINT16_T_ARRAY)
    else if STR_TO_AUX_TYPE(str, err, type, len,    uint32_t,   SLOW5_UINT32_T, SLOW5_UINT32_T_ARRAY)
    else if STR_TO_AUX_TYPE(str, err, type, len,    uint64_t,   SLOW5_UINT64_T, SLOW5_UINT64_T_ARRAY)
    else if STR_TO_AUX_TYPE(str, err, type, len,    float,      SLOW5_FLOAT,    SLOW5_FLOAT_ARRAY)
    else if STR_TO_AUX_TYPE(str, err, type, len,    double,     SLOW5_DOUBLE,   SLOW5_DOUBLE_ARRAY)
    else if STR_TO_AUX_TYPE(str, err, type, len,    char,       SLOW5_CHAR,     SLOW5_STRING)
    else if STR_TO_AUX_TYPE(str, err, type, len,    enum,       SLOW5_ENUM,     SLOW5_ENUM_ARRAY)
    else {
        *err = -1;
    }

    return type;
}

#define MEMCPY_TYPE_FROM_STR(data, value, type, err, raw_type, slow5_type, conv) \
    (type == slow5_type) { \
        raw_type value_conv = conv(value, &err); \
        if (err != -1) { \
            memcpy(data, &value_conv, sizeof value_conv); \
        } \
    }
/*
 * memcpy to data the type equivalent of the string value
 * data, value cannot be NULL
 * type shouldn't be a pointer
 * undefined behaviour if data doesn't have enough space for its type
 * returns 0 on success, -1 on error
 */
int slow5_memcpy_type_from_str(uint8_t *data, const char *value, enum slow5_aux_type type) {
    int err = -1;

    if (strcmp(value, SLOW5_ASCII_MISSING) == 0) {
        slow5_memcpy_null_type(data, type);
        err = 0;
    } else if MEMCPY_TYPE_FROM_STR(data, value, type, err,  int8_t,     SLOW5_INT8_T,   slow5_ato_int8)
    else if MEMCPY_TYPE_FROM_STR(data, value, type, err,    uint8_t,    SLOW5_UINT8_T,  slow5_ato_uint8)
    else if MEMCPY_TYPE_FROM_STR(data, value, type, err,    int16_t,    SLOW5_INT16_T,  slow5_ato_int16)
    else if MEMCPY_TYPE_FROM_STR(data, value, type, err,    uint16_t,   SLOW5_UINT16_T, slow5_ato_uint16)
    else if MEMCPY_TYPE_FROM_STR(data, value, type, err,    int32_t,    SLOW5_INT32_T,  slow5_ato_int32)
    else if MEMCPY_TYPE_FROM_STR(data, value, type, err,    uint32_t,   SLOW5_UINT32_T, slow5_ato_uint32)
    else if MEMCPY_TYPE_FROM_STR(data, value, type, err,    int64_t,    SLOW5_INT64_T,  slow5_ato_int64)
    else if MEMCPY_TYPE_FROM_STR(data, value, type, err,    uint64_t,   SLOW5_UINT64_T, slow5_ato_uint64)
    else if MEMCPY_TYPE_FROM_STR(data, value, type, err,    float,      SLOW5_FLOAT,    slow5_strtof_check)
    else if MEMCPY_TYPE_FROM_STR(data, value, type, err,    double,     SLOW5_DOUBLE,   slow5_strtod_check)
    else if (type == SLOW5_CHAR &&
            strnlen(value, 2) == 1) { /* don't allow empty string "" or "\0" or more than one character */
        err = 0;
        char value_conv = value[0];
        memcpy(data, &value_conv, sizeof value_conv);
    } else if MEMCPY_TYPE_FROM_STR(data, value, type, err,  uint8_t,    SLOW5_ENUM,     slow5_ato_uint8)

    return err;
}

#define MEMCPY_NULL_TYPE(data, type, raw_type, slow5_type, null_val) \
    (type == slow5_type) { \
        raw_type null = null_val; \
        memcpy(data, &null, sizeof null); \
    }
/*
 * write to data the null representative of the type
 * data cannot be NULL
 * type should only be primitive
 * undefined behaviour if data doesn't have enough space for its type
 */
void slow5_memcpy_null_type(uint8_t *data, enum slow5_aux_type type) {

    if MEMCPY_NULL_TYPE(data, type,         int8_t,     SLOW5_INT8_T,   SLOW5_INT8_T_NULL)
    else if MEMCPY_NULL_TYPE(data, type,    uint8_t,    SLOW5_UINT8_T,  SLOW5_UINT8_T_NULL)
    else if MEMCPY_NULL_TYPE(data, type,    int16_t,    SLOW5_INT16_T,  SLOW5_INT16_T_NULL)
    else if MEMCPY_NULL_TYPE(data, type,    uint16_t,   SLOW5_UINT16_T, SLOW5_UINT16_T_NULL)
    else if MEMCPY_NULL_TYPE(data, type,    int32_t,    SLOW5_INT32_T,  SLOW5_INT32_T_NULL)
    else if MEMCPY_NULL_TYPE(data, type,    uint32_t,   SLOW5_UINT32_T, SLOW5_UINT32_T_NULL)
    else if MEMCPY_NULL_TYPE(data, type,    int64_t,    SLOW5_INT64_T,  SLOW5_INT64_T_NULL)
    else if MEMCPY_NULL_TYPE(data, type,    uint64_t,   SLOW5_UINT64_T, SLOW5_UINT64_T_NULL)
    else if MEMCPY_NULL_TYPE(data, type,    float,      SLOW5_FLOAT,    SLOW5_FLOAT_NULL)
    else if MEMCPY_NULL_TYPE(data, type,    double,     SLOW5_DOUBLE,   SLOW5_DOUBLE_NULL)
    else if MEMCPY_NULL_TYPE(data, type,    char,       SLOW5_CHAR,     SLOW5_CHAR_NULL)
    else if MEMCPY_NULL_TYPE(data, type,    uint8_t,    SLOW5_ENUM,     SLOW5_ENUM_NULL)
}


/*
 * get malloced string representing missing data
 * len cannot be NULL
 * on malloc error sets *len to -1 and returns NULL
 */
static char *get_missing_str(size_t *len) {
    char *str = malloc(2 * sizeof *str);
    SLOW5_MALLOC_CHK(str);
    if (str == NULL) {
        *len = -1;
        return str;
    }

    str[0] = SLOW5_ASCII_MISSING_CHAR;
    str[1] = '\0';
    *len = 1;

    return str;
}

/*
 * convert data to string representation
 * ptr_len is the number of elements of data
 * data, str_len cannot be NULL
 * undefined behaviour if type doesn't match or ptr_len is too long for data given the type is a pointer
 * sets str_len to -1 on error and returns NULL
 * TODO use macro to clean up
 */
char *slow5_data_to_str(uint8_t *data, enum slow5_aux_type type, uint64_t ptr_len, size_t *str_len) {
    char *str = NULL;

    if (type == SLOW5_INT8_T) {
        if (*(int8_t *) data == SLOW5_INT8_T_NULL) {
            str = get_missing_str(str_len);
            // sets *str_len to -1 and str to NULL on malloc error
        } else {
            *str_len = slow5_asprintf(&str, "%" PRId8, *(int8_t *) data);
            // sets *str_len to negative and doesn't change str on error
        }
    } else if (type == SLOW5_UINT8_T) {
        if (*(uint8_t *) data == SLOW5_UINT8_T_NULL) {
            str = get_missing_str(str_len);
        } else {
            *str_len = slow5_asprintf(&str, "%" PRIu8, *(uint8_t *) data);
        }
    } else if (type == SLOW5_INT16_T) {
        if (*(int16_t *) data == SLOW5_INT16_T_NULL) {
            str = get_missing_str(str_len);
        } else {
            *str_len = slow5_asprintf(&str, "%" PRId16, *(int16_t *) data);
        }
    } else if (type == SLOW5_UINT16_T) {
        if (*(uint16_t *) data == SLOW5_UINT16_T_NULL) {
            str = get_missing_str(str_len);
        } else {
            *str_len = slow5_asprintf(&str, "%" PRIu16, *(uint16_t *) data);
        }
    } else if (type == SLOW5_INT32_T) {
        if (*(int32_t *) data == SLOW5_INT32_T_NULL) {
            str = get_missing_str(str_len);
        } else {
            *str_len = slow5_asprintf(&str, "%" PRId32, *(int32_t *) data);
        }
    } else if (type == SLOW5_UINT32_T) {
        if (*(uint32_t *) data == SLOW5_UINT32_T_NULL) {
            str = get_missing_str(str_len);
        } else {
            *str_len = slow5_asprintf(&str, "%" PRIu32, *(uint32_t *) data);
        }
    } else if (type == SLOW5_INT64_T) {
        if (*(int64_t *) data == SLOW5_INT64_T_NULL) {
            str = get_missing_str(str_len);
        } else {
            *str_len = slow5_asprintf(&str, "%" PRId64, *(int64_t *) data);
        }
    } else if (type == SLOW5_UINT64_T) {
        if (*(uint64_t *) data == SLOW5_UINT64_T_NULL) {
            str = get_missing_str(str_len);
        } else {
            *str_len = slow5_asprintf(&str, "%" PRIu64, *(uint64_t *) data);
        }
    } else if (type == SLOW5_FLOAT) {
        if (isnan(*(float *) data)) {
            str = get_missing_str(str_len);
        } else {
            str = slow5_float_to_str(*(float *) data, str_len);
        }
    } else if (type == SLOW5_DOUBLE) {
        if (isnan(*(double *) data)) {
            str = get_missing_str(str_len);
        } else {
            str = slow5_double_to_str(*(double *) data, str_len);
        }
    } else if (type == SLOW5_CHAR) {
        if (*(char *) data == SLOW5_CHAR_NULL) {
            str = get_missing_str(str_len);
        } else {
            str = (char *) malloc(sizeof (char));
            SLOW5_MALLOC_CHK(str);

            if (str == NULL) {
                *str_len = -1;
            } else {
                *str_len = sizeof (char);
                str[0] = *(char *) data;
            }
        }
    } else if (type == SLOW5_ENUM) {
        if (*(uint8_t *) data == SLOW5_ENUM_NULL) {
            str = get_missing_str(str_len);
        } else {
            *str_len = slow5_asprintf(&str, "%" PRIu8, *(uint8_t *) data);
        }
    } else if (type == SLOW5_STRING) {
        if (ptr_len == 0) {
            str = get_missing_str(str_len);
        } else {
            str = strdup((char *) data);
            *str_len = strlen(str);
        }
    } else if (SLOW5_IS_PTR(type)) {
        if (ptr_len == 0) {
            str = get_missing_str(str_len);
        } else {
            size_t str_cap = SLOW5_AUX_ARRAY_STR_CAP_INIT;
            str = (char *) malloc(str_cap * sizeof *str);
            SLOW5_MALLOC_CHK(str);

            size_t str_cur = 0;
            uint64_t i;
            for (i = 0; i < ptr_len - 1; ++ i) {
                size_t str_len_sep;
                char *str_sep = slow5_data_to_str(data + i * SLOW5_AUX_TYPE_META[type].size, SLOW5_TO_PRIM_TYPE(type), 1, &str_len_sep);

                if (str_cur + str_len_sep + 1 > str_cap) { // +1 for SLOW5_SEP_ARRAY_CHAR
                    // Realloc
                    str_cap = str_cap << 1;
                    str = (char *) realloc(str, str_cap * sizeof *str);
                    SLOW5_MALLOC_CHK(str);
                }

                memcpy(str + str_cur, str_sep, str_len_sep);
                str[str_cur + str_len_sep] = SLOW5_SEP_ARRAY_CHAR;
                str_cur += str_len_sep + 1;

                free(str_sep);
            }
            size_t str_len_sep;
            char *str_sep = slow5_data_to_str(data + i * SLOW5_AUX_TYPE_META[type].size, SLOW5_TO_PRIM_TYPE(type), 1, &str_len_sep);

            if (str_cur + str_len_sep + 1 > str_cap) { // +1 for '\0'
                // Realloc
                str_cap = str_cap << 1;
                str = (char *) realloc(str, str_cap * sizeof *str);
                SLOW5_MALLOC_CHK(str);
            }

            memcpy(str + str_cur, str_sep, str_len_sep);
            str[str_cur + str_len_sep] = '\0';

            *str_len = str_cur + str_len_sep;
            free(str_sep);
        }
    }

    //if (str == NULL)

    return str;
}

/* set the logging level */
void slow5_set_log_level(enum slow5_log_level_opt log_level) {
    slow5_log_level = log_level;
}

/* set the condition of exit (on error, warning, neither) */
void slow5_set_exit_condition(enum slow5_exit_condition_opt exit_condition) {
    slow5_exit_condition = exit_condition;
}

/*
 * is slow5 file at end? seek back, read and find out
 * return 0 if not at end and file pointer left unchanged
 * return 1 if at end
 * return -1 and set slow5_errno on error
 * return -2 and set slow5_errno=SLOW5_ERR_TRUNC if eof found but not at end
 **/
int slow5_is_eof(FILE *fp, const char *eof, size_t n) {
    if (!fp) {
        SLOW5_ERROR("Argument '%s' cannot be NULL.", SLOW5_TO_STR(fp));
        slow5_errno = SLOW5_ERR_ARG;
        return -1;
    }

    char *buf_eof = (char *) malloc(sizeof *eof * n);
    if (!buf_eof) {
        SLOW5_MALLOC_ERROR();
        slow5_errno = SLOW5_ERR_MEM;
        return -1;
    }

    /* go back */
    if (fseek(fp, - sizeof *eof * n, SEEK_CUR) != 0) {
        SLOW5_ERROR("Seeking back '%zu' bytes failed: %s.",
                sizeof *eof * n, strerror(errno));
        free(buf_eof);
        slow5_errno = SLOW5_ERR_IO;
        return -1;
    }

    size_t itms_read = fread(buf_eof, sizeof *eof, n, fp);
    if (itms_read == n) {
        if (memcmp(eof, buf_eof, sizeof *eof * n) == 0) { /* eof marker found */
            if (getc(fp) != EOF || !feof(fp)) { /* not actually at end of file */
                free(buf_eof);
                slow5_errno = SLOW5_ERR_TRUNC;
                return -2;
            }
            free(buf_eof);
            return 1;
        }
    }

    free(buf_eof);
    return 0;
}

/* following can be moved to version.c file later when it gows out of control */
/**
 * compare file versions
 * return <0 if x < y
 * return 0 if x == y
 * return >0 if x > y
 */
int slow5_version_cmp(struct slow5_version x, struct slow5_version y) {
    if (x.major > y.major ||
            (x.major == y.major && x.minor > y.minor) ||
            (x.major == y.major && x.minor == y.minor && x.patch > y.patch)) {
        return 1;
    } else if (x.major == y.major && x.minor == y.minor && x.patch == y.patch) {
        return 0;
    } else {
        return -1;
    }
}

// if the slow5 version >= 0.2.0 where signal_press was introduced
int slow5_signal_press_version_cmp(struct slow5_version current){
    struct slow5_version signal_press_version = { .major = 0, .minor = 2, .patch = 0 }; //signal_press related flag in blow5 header was introduced in version 0.2.0
    return slow5_version_cmp(current,signal_press_version);
}

int slow5_version_sanity(struct slow5_hdr *hdr){
    struct slow5_version current = hdr->version;
    struct slow5_version enum_version = { .major = 0, .minor = 2, .patch = 0 }; //enums were introduced in version 0.2.0
    if(slow5_version_cmp(current, enum_version) < 0 && hdr->aux_meta &&
        (hdr->aux_meta->enum_labels != NULL || hdr->aux_meta->enum_num_labels != NULL))
    {
        SLOW5_WARNING("You file version '" SLOW5_VERSION_STRING_FORMAT "' has an enum datatype which was only introduced in version '" SLOW5_VERSION_STRING_FORMAT "'",
                      current.major, current.minor, current.patch, enum_version.major, enum_version.minor, enum_version.patch);
        return 1;

    }
    return 0;
}

//bump the slow5 file version to a newer version if a compression method unavailable in the current file version was requested
struct slow5_version slow5_press_version_bump(struct slow5_version current, slow5_press_method_t method){

    struct slow5_version signal_press_version = { .major = 0, .minor = 2, .patch = 0 };
    if(slow5_version_cmp(current,signal_press_version) < 0 &&
        (method.record_method == SLOW5_COMPRESS_SVB_ZD || method.record_method == SLOW5_COMPRESS_ZSTD ||
         method.signal_method == SLOW5_COMPRESS_SVB_ZD || method.signal_method == SLOW5_COMPRESS_ZSTD )
    ){
        SLOW5_INFO("SLOW5 version updated to '" SLOW5_VERSION_STRING_FORMAT "' as the requested compression option is unavailable in the current fileversion '" SLOW5_VERSION_STRING_FORMAT "'",
            signal_press_version.major, signal_press_version.minor, signal_press_version.patch,
            current.major, current.minor, current.patch);
        return signal_press_version;
    }
    else{
        return current;
    }

}

//int main(void) {

/*
slow5_f = slow5_open("../test/data/out/a.out/test.slow5", "w");
slow5_write_hdr(slow5_f);
slow5_close(slow5_f);
*/

/*
slow5 = slow5_open("../test/data/out/a.out/test.slow5", "r");
slow5_read_hdr(slow5);
slow5_close(slow5);
*/

/*
 * slow5_file *f_in = slow5_fopen("hi.slow5", "r");
 * slow5_file *f_out = slow5_fopen("hi.blow5", "w");
 */
//FILE *f_in = fopen("../test/data/err/version_too_large.slow5", "r");
//FILE *f_out = fopen("hi.blow5", "w");

//   struct slow5_file *s5p = slow5_open("../test/data/exp/one_fast5/exp_1.slow5", "r");
//struct SLOW5 *slow5 = slow5_init_empty(void);
//slow5_read(slow5, SLOW5_FORMAT_ASCII, f_in);

/*
 * slow5_write(f_in, f_out);
 */
/*
struct SLOW5WriteConf *conf = slow5_wconf_init(SLOW5_FORMAT_BINARY, SLOW5_COMPRESS_ZLIB);
slow5_write(slow5, conf, f_out);

slow5_wconf_destroy(&conf);
*/
/*
    slow5_hdr_print(s5p->header);
    struct slow5_rec *rec = NULL;
    slow5_get("a649a4ae-c43d-492a-b6a1-a5b8b8076be4", &rec, s5p);
    slow5_rec_free(rec);
    slow5_close(s5p);
    */

//fclose(f_in);
//fclose(f_out);
//}

/*
slow5_convert(Format from_format, Format to_format, FILE *from_file, FILE *to_file) {
}

f2s_single(FILE *fast5, FILE *slow5) {
}

f2s() {
    f2s_single();
    merge();
}
*/
