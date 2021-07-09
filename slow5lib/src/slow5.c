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
//#include "slow5_error.h"
//#include "slow5_press.h"
#include "slow5_misc.h"
#include "klib/ksort.h"
//#include "klib/kvec.h"
//#include "klib/khash.h"

KSORT_INIT_STR

// TODO fail with getline if end of file occurs on a non-empty line
// TODO (void) cast if ignoring return value
// TODO sizeof of macros at compile time rather than strlen
// TODO Make sure all mallocs are checked for success
// TODO replace SLOW5_ASSERT with proper error messages and handling

// Initial string buffer capacity for parsing the data header
#define SLOW5_HEADER_DATA_BUF_INIT_CAP (1024) // 2^10 TODO is this too much? Or put to a page length
// Max length is 6 (−32768) for a int16_t
#define INT16_MAX_LENGTH (6)
// Max length is 10 (4294967295) for a uint32_t
#define UINT32_MAX_LENGTH (10)
// Fixed string buffer capacity for storing signal
#define SLOW5_SIGNAL_BUF_FIXED_CAP (8) // 2^3 since INT16_MAX_LENGTH=6
// Initial capacity for converting the header to a string
#define SLOW5_HEADER_STR_INIT_CAP (1024) // 2^10 TODO is this good? Or put to a page length
// Initial capacity for the number of auxiliary fields
#define SLOW5_AUX_META_CAP_INIT (32) // 2^5 TODO is this good? Too small?
// Initial capacity for parsing auxiliary array
#define SLOW5_AUX_ARRAY_CAP_INIT (256) // 2^8
// Initial capacity for storing auxiliary array string
#define SLOW5_AUX_ARRAY_STR_CAP_INIT (1024) // 2^10

static inline void slow5_free(struct slow5_file *s5p);
static inline khash_t(slow5_s2a) *slow5_rec_aux_init(void);
static inline void slow5_rec_set_aux_map(khash_t(slow5_s2a) *aux_map, const char *field, const uint8_t *data, size_t len, uint64_t bytes, enum slow5_aux_type type);


enum slow5_log_level_opt  slow5_log_level = SLOW5_LOG_INFO;
enum slow5_exit_condition_opt  slow5_exit_condition = SLOW5_EXIT_OFF;

__thread int slow5_errno = 0;

/* Definitions */

// slow5 file


/*
* core function for opening a SLOW5 file
* determines if SLOW5 ASCII or binary based on the file extension
* parses and populates the header
* initialises decompression buffers
*/
struct slow5_file *slow5_init(FILE *fp, const char *pathname, enum slow5_fmt format) {
    // Pathname cannot be NULL at this point
    if (fp == NULL) {
        SLOW5_WARNING("%s","Cannot initialise a SLOW5 file with a NULL file pointer.");
        return NULL;
    }

    if (format == SLOW5_FORMAT_UNKNOWN) {

        // Attempt to determine format
        // from pathname
        if ((format = slow5_path_get_fmt(pathname)) == SLOW5_FORMAT_UNKNOWN) {
            fclose(fp);
            SLOW5_ERROR("%s","Could not determine SLOW5 file format. Extension must be .slow5 or .blow5");
            return NULL;
        }
    }

    struct slow5_file *s5p;
    slow5_press_method_t method;
    struct slow5_hdr *header = slow5_hdr_init(fp, format, &method);
    if (header == NULL) {
        fclose(fp);
        SLOW5_ERROR("%s","Error parsing SLOW5 header.");
        s5p = NULL;
    } else {
        s5p = (struct slow5_file *) calloc(1, sizeof *s5p);

        s5p->fp = fp;
        s5p->format = format;
        s5p->header = header;
        s5p->compress = slow5_press_init(method);

        if ((s5p->meta.fd = fileno(fp)) == -1) {
            SLOW5_ERROR("%s","Obtaining the fileno() for the file stream failed. Not a valid file stream.");
            slow5_close(s5p);
            s5p = NULL;
        }
        s5p->meta.pathname = pathname;
        s5p->meta.start_rec_offset = ftello(fp);
        if(s5p->meta.start_rec_offset == -1 ){
            SLOW5_ERROR("ftello() failed. %s",strerror(errno));
            slow5_close(s5p);
            s5p = NULL;
        }
    }

    return s5p;
}

// TODO this needs to be refined: talk to Sasha (he wrote this here)
struct slow5_file *slow5_init_empty(FILE *fp, const char *pathname, enum slow5_fmt format) {

    if(slow5_is_big_endian()){
        SLOW5_ERROR("%s","Big endian machine detected. SLOW5lib only support little endian at this time. Please open a github issue stating your machine spec.");
        return NULL;
    }
    // Pathname cannot be NULL at this point
    if (fp == NULL) {
        return NULL;
    }

    if (format == SLOW5_FORMAT_UNKNOWN) {

        // Attempt to determine format
        // from pathname
        if ((format = slow5_path_get_fmt(pathname)) == SLOW5_FORMAT_UNKNOWN) {
            fclose(fp);
            return NULL;
        }
    }

    struct slow5_file *s5p;
    struct slow5_hdr *header = slow5_hdr_init_empty();
    header->version = SLOW5_ASCII_VERSION_STRUCT;
    s5p = (struct slow5_file *) calloc(1, sizeof *s5p);

    s5p->fp = fp;
    s5p->format = format;
    s5p->header = header;

    if ((s5p->meta.fd = fileno(fp)) == -1) {
        slow5_close(s5p);
        s5p = NULL;
    }
    s5p->meta.pathname = pathname;
    s5p->meta.start_rec_offset = ftello(fp);

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
 * The user at the moment is expected to give "r"
 * TODO : Make "r" into "rb" if BLOW5 - this is not an issue for POSIX systems as mode b is not used [https://man7.org/linux/man-pages/man3/fopen.3.html]
 *
 * @param   pathname    relative or absolute path to slow5 file
 * @param   mode        only "r" for the moment for reading
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
 * @param   pathname    relative or absolute path to slow5 file
 * @param   mode        only "r" for the moment
 * @param   format      format of the slow5 file
 * @return              slow5 file structure
 */
struct slow5_file *slow5_open_with(const char *pathname, const char *mode, enum slow5_fmt format) {
    if(slow5_is_big_endian()){
        SLOW5_ERROR("%s","Big endian machine detected. SLOW5lib only support little endian at this time. Please open a github issue stating your machine spec.");
        return NULL;
    }
    if (pathname == NULL) {
        SLOW5_ERROR("%s","Argument 'pathname' cannot be NULL.");
        slow5_errno = SLOW5_ERR_ARG;
        return NULL;
    }
    if (mode == NULL) {
        SLOW5_ERROR("%s","Argument 'mode' cannot be NULL.");
        slow5_errno = SLOW5_ERR_ARG;
        return NULL;
    }

    if(strcmp(mode,"r")){
        SLOW5_WARNING("Currently, the only supported mode is 'r'. You entered %s",mode);
    }

    FILE *fp = fopen(pathname, mode);
    if(fp==NULL){
        SLOW5_ERROR("Error opening file %s: %s",pathname, strerror(errno));
    }
    return slow5_init(fp, pathname, format);

}

// Close a slow5 file
int slow5_close(struct slow5_file *s5p) {
    int ret;

    if (s5p == NULL) {
        ret = EOF;
    } else {
        ret = fclose(s5p->fp);
        slow5_free(s5p);
    }

    return ret;
}

static inline void slow5_free(struct slow5_file *s5p) {
    if (s5p != NULL) {
        slow5_press_free(s5p->compress);
        slow5_hdr_free(s5p->header);

        if (s5p->index != NULL) {

            if(s5p->index->dirty){ //if the index has been changed, write it back
                if (s5p->index->fp != NULL) {
                    SLOW5_ASSERT(fclose(s5p->index->fp) == 0);
                }

                s5p->index->fp = fopen(s5p->index->pathname, "wb");
                slow5_idx_write(s5p->index);
            }
            slow5_idx_free(s5p->index);
        }

        free(s5p);
    }
}


// slow5 header

struct slow5_hdr *slow5_hdr_init_empty(void) {

    struct slow5_hdr *header = (struct slow5_hdr *) calloc(1, sizeof *(header));

    return header;
}



static inline int slow5_is_version_compatible(struct slow5_version file_version){
    if(file_version.major > SLOW5_ASCII_VERSION_STRUCT.major){
        return 0;
    }
    else if (file_version.minor > SLOW5_ASCII_VERSION_STRUCT.minor){
        return 0;
    }
    else if (file_version.patch > SLOW5_ASCII_VERSION_STRUCT.patch){
        return 0;
    }
    else{
        return 1;
    }
}

// parses a slow5 header
struct slow5_hdr *slow5_hdr_init(FILE *fp, enum slow5_fmt format, slow5_press_method_t *method) {

    struct slow5_hdr *header = (struct slow5_hdr *) calloc(1, sizeof *(header));
    char *buf = NULL;

    // Parse slow5 header

    if (format == SLOW5_FORMAT_ASCII) {

        *method = SLOW5_COMPRESS_NONE;

        // Buffer for file parsing
        size_t cap = SLOW5_HEADER_DATA_BUF_INIT_CAP;
        buf = (char *) malloc(cap * sizeof *buf);
        char *bufp;
        ssize_t buf_len;
        int err;

        // 1st line - slow5_version
        if ((buf_len = getline(&buf, &cap, fp)) == -1) {
            SLOW5_WARNING("%s", "'getline()' failed. Is this an empty file?");
            free(buf);
            free(header);
            return NULL;
        }
        if(buf_len==0){
            SLOW5_WARNING("%s", "'Malformed SLOW5 header. Why is the very first line empty?");
            free(buf);
            free(header);
            return NULL;
        }
        buf[buf_len - 1] = '\0'; // Remove newline for later parsing
        // "#slow5_version"
        bufp = buf;
        char *tok = slow5_strsep(&bufp, SLOW5_SEP_COL);
        if (strcmp(tok, SLOW5_HEADER_FILE_VERSION_ID) != 0) {
            SLOW5_WARNING("Malformed SLOW5 header. Expected '%s', instead found '%s'", SLOW5_HEADER_FILE_VERSION_ID, tok);
            free(buf);
            free(header);
            return NULL;
        }
        // Parse file version
        tok = slow5_strsep(&bufp, SLOW5_SEP_COL);
        char *toksub;
        if ((toksub = slow5_strsep(&tok, ".")) == NULL) { // Major version
            SLOW5_WARNING("%s", "Malformed SLOW5 header. Version string is expected to be in the format x.y.z");
            free(buf);
            free(header);
            return NULL;
        }
        header->version.major = slow5_ato_uint8(toksub, &err);
        if (err == -1 || (toksub = slow5_strsep(&tok, ".")) == NULL) { // Minor version
            SLOW5_WARNING("%s", "Malformed SLOW5 header. Version string is expected to be in the format x.y.z");
            free(buf);
            free(header);
            return NULL;
        }
        header->version.minor = slow5_ato_uint8(toksub, &err);
        if (err == -1 || (toksub = slow5_strsep(&tok, ".")) == NULL) { // Patch version
            SLOW5_WARNING("%s", "Malformed SLOW5 header. Version string is expected to be in the format x.y.z");
            free(buf);
            free(header);
            return NULL;
        }
        header->version.patch = slow5_ato_uint8(toksub, &err);
        if (err == -1 || slow5_strsep(&tok, ".") != NULL) { // No more tokenators
            SLOW5_WARNING("%s", "Malformed SLOW5 header. Version string is expected to be in the format x.y.z");
            free(buf);
            free(header);
            return NULL;
        }
        if(slow5_is_version_compatible(header->version)==0){
            SLOW5_ERROR("file version (%d.%d.%d) in your slow5 file is higher than the maximally compatible version (%d.%d.%d) by this slow5lib. Please use a newer version of slow5lib",
                header->version.major, header->version.minor, header->version.patch,
                     SLOW5_ASCII_VERSION_STRUCT.major,  SLOW5_ASCII_VERSION_STRUCT.minor,  SLOW5_ASCII_VERSION_STRUCT.patch);
            free(buf);
            free(header);
            return NULL;
        }

        // 2nd line - num_read_groups
        if ((buf_len = getline(&buf, &cap, fp)) <= 0) {
            SLOW5_WARNING("%s", "'getline()' failed. Where is the num_read_groups?");
            free(buf);
            free(header);
            return NULL;
        }
        buf[buf_len - 1] = '\0'; // Remove newline for later parsing
        // "#num_read_groups"
        bufp = buf;
        tok = slow5_strsep(&bufp, SLOW5_SEP_COL);
        if (strcmp(tok, SLOW5_HEADER_NUM_GROUPS_ID) != 0) {
            SLOW5_WARNING("Malformed SLOW5 header. Expected '%s', instead found '%s'", SLOW5_HEADER_NUM_GROUPS_ID, tok);
            free(buf);
            free(header);
            return NULL;
        }
        // Parse num read groups
        tok = slow5_strsep(&bufp, SLOW5_SEP_COL);
        header->num_read_groups = slow5_ato_uint32(tok, &err);
        if (err == -1) {
            SLOW5_WARNING("Malformed SLOW5 header. Invalid number of read groups - '%s'", tok);
            free(buf);
            free(header);
            return NULL;
        }

        //parse the data header (TODO: assertions inside this function must be replaced properly)
        if(slow5_hdr_data_init(fp, buf, &cap, header, NULL) != 0){
            SLOW5_WARNING("%s","Parsing the data header failed.");
            free(buf);
            free(header);
            return NULL;
        }
        //parse the datatype and column name fields (TODO: error checking inside this function is really bad - fix it)
        header->aux_meta = slow5_aux_meta_init(fp, buf, &cap, NULL);

    } else if (format == SLOW5_FORMAT_BINARY) {
        const char magic[] = SLOW5_BINARY_MAGIC_NUMBER;

        char buf_magic[sizeof magic];
        uint32_t header_size;

        // TODO pack and do one read

        if (fread(buf_magic, sizeof *magic, sizeof magic, fp) != sizeof magic ||
                memcmp(magic, buf_magic, sizeof *magic * sizeof magic) != 0 ||
                fread(&header->version.major, sizeof header->version.major, 1, fp) != 1 ||
                fread(&header->version.minor, sizeof header->version.minor, 1, fp) != 1 ||
                fread(&header->version.patch, sizeof header->version.patch, 1, fp) != 1 ||
                fread(method, sizeof *method, 1, fp) != 1 ||
                fread(&header->num_read_groups, sizeof header->num_read_groups, 1, fp) != 1 ||
                fseek(fp, SLOW5_BINARY_HEADER_SIZE_OFFSET, SEEK_SET) == -1 ||
                fread(&header_size, sizeof header_size, 1, fp) != 1) {

            SLOW5_WARNING("%s","Malformed BLOW5 header");
            free(header);
            return NULL;
        }

        if(slow5_is_version_compatible(header->version)==0){
            SLOW5_ERROR("file version (%d.%d.%d) in your slow5 file is higher than the maximally compatible version (%d.%d.%d) by this slow5lib. Please use a newer version of slow5lib",
                header->version.major, header->version.minor, header->version.patch,
                     SLOW5_ASCII_VERSION_STRUCT.major,  SLOW5_ASCII_VERSION_STRUCT.minor,  SLOW5_ASCII_VERSION_STRUCT.patch);
            free(header);
            return NULL;
        }

        size_t cap = SLOW5_HEADER_DATA_BUF_INIT_CAP;
        buf = (char *) malloc(cap * sizeof *buf);

        // Header data
        uint32_t header_act_size;
        SLOW5_ASSERT(slow5_hdr_data_init(fp, buf, &cap, header, &header_act_size) == 0);
        header->aux_meta = slow5_aux_meta_init(fp, buf, &cap, &header_act_size);
        if (header_act_size != header_size) {
            slow5_hdr_free(header);
            return NULL;
        }
    }

    free(buf);

    return header;
}


//get the list of hdr data keys in sorted order (only the returned pointer must be freed, not the ones inside - subject to change)
//currently used in pyslow5
//len is the number of elements in the return char** list
//returns null if there are no data header attributes
const char **slow5_get_hdr_keys(const slow5_hdr_t *header,uint64_t *len){
    if(len!=NULL){
        *len = header->data.num_attrs;
    }
    if(header->data.num_attrs == 0){
        return NULL;
    }
    const char **data_attrs = (const char **) malloc(header->data.num_attrs * sizeof *data_attrs);
    uint32_t i = 0;
    for (khint_t j = kh_begin(header->data.attrs); j != kh_end(header->data.attrs); ++ j) {
        if (kh_exist(header->data.attrs, j)) {
            data_attrs[i] = kh_key(header->data.attrs, j);
            ++ i;
        }
    }

    // Sort header data attributes alphabetically
    ks_mergesort(str, header->data.num_attrs, data_attrs, 0);

    return data_attrs;
}

/**
 * Get the header as a string in the specified format.
 *
 * Returns NULL if s5p is NULL
 * or format is SLOW5_FORMAT_UNKNOWN
 * or an internal error occurs.
 *
 * @param   s5p     slow5 file
 * @param   format  slow5 format to write the entry in
 * @param   comp    compression method
 * @param   n       number of bytes written to the returned buffer
 * @return  malloced memory storing the slow5 header representation,
 *          to use free() on afterwards
 */
// TODO don't allow comp of SLOW5_COMPRESS_GZIP for SLOW5_FORMAT_ASCII

//  flattened header returned as a void * (incase of BLOW5 magic number is also included)

void *slow5_hdr_to_mem(struct slow5_hdr *header, enum slow5_fmt format, slow5_press_method_t comp, size_t *n) {
    char *mem = NULL;

    if (header == NULL || format == SLOW5_FORMAT_UNKNOWN) {
        return mem;
    }

    size_t len = 0;
    size_t cap = SLOW5_HEADER_STR_INIT_CAP;
    mem = (char *) malloc(cap * sizeof *mem);
    uint32_t header_size;

    if (format == SLOW5_FORMAT_ASCII) {

        struct slow5_version *version = &header->version;

        // Relies on SLOW5_HEADER_DATA_BUF_INIT_CAP being bigger than
        // strlen(ASCII_SLOW5_HEADER) + UINT32_MAX_LENGTH + strlen("\0")
        int len_ret = sprintf(mem, SLOW5_ASCII_SLOW5_HEADER_FORMAT,
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

        struct slow5_version *version = &header->version;

        // Relies on SLOW5_HEADER_DATA_BUF_INIT_CAP
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
        memcpy(mem + len, &comp, sizeof comp);
        len += sizeof comp;
        memcpy(mem + len, &header->num_read_groups, sizeof header->num_read_groups);
        len += sizeof header->num_read_groups;

        memset(mem + len, '\0', SLOW5_BINARY_HEADER_SIZE_OFFSET - len);
        len = SLOW5_BINARY_HEADER_SIZE_OFFSET;

        // Skip header size for later
        len += sizeof header_size;
    }

    size_t len_to_cp;
    // Get unsorted list of header data attributes.
    if (header->data.num_attrs != 0) {
        const char **data_attrs = slow5_get_hdr_keys(header,NULL);

        // Write header data attributes to string
        for (size_t i = 0; i < header->data.num_attrs; ++ i) {
            const char *attr = data_attrs[i];

            // Realloc if necessary
            if (len + 1 + strlen(attr) >= cap) { // + 1 for SLOW5_HEADER_DATA_PREFIX_CHAR
                cap *= 2;
                mem = (char *) realloc(mem, cap * sizeof *mem);
            }

            mem[len] = SLOW5_HEADER_DATA_PREFIX_CHAR;
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
                    }

                    mem[len] = SLOW5_SEP_COL_CHAR;
                    ++ len;

                    if (value != NULL) {
                        len_to_cp = strlen(value);

                        //special case for "."
                        if(strlen(value)==0){
                            len_to_cp++;
                        }

                        // Realloc if necessary
                        if (len + len_to_cp >= cap) {
                            cap *= 2;
                            mem = (char *) realloc(mem, cap * sizeof *mem);
                        }

                        if(strlen(value)==0){ //special case for "."
                            memcpy(mem + len, ".", len_to_cp);
                        }
                        else {
                            memcpy(mem + len, value, len_to_cp);
                        }
                        len += len_to_cp;
                    }
                } else {
                    // Realloc if necessary
                    if (len + 1 >= cap) { // +1 for SLOW5_SEP_COL_CHAR
                        cap *= 2;
                        mem = (char *) realloc(mem, cap * sizeof *mem);
                    }

                    mem[len] = SLOW5_SEP_COL_CHAR;
                    ++ len;
                }
            }

            // Realloc if necessary
            if (len + 1 >= cap) { // +1 for '\n'
                cap *= 2;
                mem = (char *) realloc(mem, cap * sizeof *mem);
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
    }
    memcpy(mem + len, str_to_cp, len_to_cp);
    len += len_to_cp;
    free(str_to_cp);

    if (format == SLOW5_FORMAT_ASCII) {
        // Realloc if necessary
        if (len + 1 >= cap) { // +1 for '\0'
            cap *= 2;
            mem = (char *) realloc(mem, cap * sizeof *mem);
        }

        mem[len] = '\0';
    } else if (format == SLOW5_FORMAT_BINARY) { //write the header size in bytes (which was skipped previously)
        header_size = len - (SLOW5_BINARY_HEADER_SIZE_OFFSET + sizeof header_size);
        memcpy(mem + SLOW5_BINARY_HEADER_SIZE_OFFSET, &header_size, sizeof header_size);
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
        size_t types_cap = SLOW5_HEADER_STR_INIT_CAP;
        types = (char *) malloc(types_cap);

        // Assumption that SLOW5_HEADER_STR_INIT_CAP > strlen(SLOW5_ASCII_TYPE_HEADER_MIN)
        const char *str_to_cp = SLOW5_ASCII_TYPE_HEADER_MIN;
        size_t len_to_cp = strlen(str_to_cp);
        memcpy(types, str_to_cp, len_to_cp);
        types_len += len_to_cp;

        for (uint16_t i = 0; i < aux_meta->num; ++ i) {
            const char *str_to_cp = SLOW5_AUX_TYPE_META[aux_meta->types[i]].type_str;
            len_to_cp = strlen(str_to_cp);

            if (types_len + len_to_cp + 1 >= types_cap) { // +1 for SLOW5_SEP_COL_CHAR
                types_cap *= 2;
                types = (char *) realloc(types, types_cap);
            }

            types[types_len ++] = SLOW5_SEP_COL_CHAR;
            memcpy(types + types_len, str_to_cp, len_to_cp);
            types_len += len_to_cp;
        }

        if (types_len + 2 >= types_cap) { // +2 for '\n' and '\0'
            types_cap *= 2;
            types = (char *) realloc(types, types_cap);
        }
        types[types_len ++] = '\n';
        types[types_len] = '\0';

    } else {
        types = strdup(SLOW5_ASCII_TYPE_HEADER_MIN "\n");
        types_len = strlen(types);
    }

    *len = types_len;

    return types;
}

char *slow5_hdr_attrs_to_str(struct slow5_aux_meta *aux_meta, size_t *len) {
    char *attrs = NULL;
    size_t attrs_len = 0;

    if (aux_meta != NULL) {
        size_t attrs_cap = SLOW5_HEADER_STR_INIT_CAP;
        attrs = (char *) malloc(attrs_cap);

        // Assumption that SLOW5_HEADER_STR_INIT_CAP > strlen(SLOW5_ASCII_TYPE_HEADER_MIN)
        const char *str_to_cp = SLOW5_ASCII_COLUMN_HEADER_MIN;
        size_t len_to_cp = strlen(str_to_cp);
        memcpy(attrs, str_to_cp, len_to_cp);
        attrs_len += len_to_cp;

        for (uint16_t i = 0; i < aux_meta->num; ++ i) {
            str_to_cp = aux_meta->attrs[i];
            len_to_cp = strlen(str_to_cp);

            if (attrs_len + len_to_cp + 1 >= attrs_cap) { // +1 for SLOW5_SEP_COL_CHAR
                attrs_cap *= 2;
                attrs = (char *) realloc(attrs, attrs_cap);
            }

            attrs[attrs_len ++] = SLOW5_SEP_COL_CHAR;
            memcpy(attrs + attrs_len, str_to_cp, len_to_cp);
            attrs_len += len_to_cp;
        }

        if (attrs_len + 2 >= attrs_cap) { // +2 for '\n' and '\0'
            attrs_cap *= 2;
            attrs = (char *) realloc(attrs, attrs_cap);
        }
        attrs[attrs_len ++] = '\n';
        attrs[attrs_len] = '\0';

    } else {
        attrs = strdup(SLOW5_ASCII_COLUMN_HEADER_MIN "\n");
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
 * @param   fp      output file pointer
 * @param   s5p     slow5_rec pointer
 * @param   format  slow5 format to write the entry in
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

char **slow5_get_aux_names(const slow5_hdr_t *header,uint64_t *len){

    u_int64_t my_len = (header->aux_meta == NULL) ? 0 : header->aux_meta->num;

    if(len){
        *len = my_len;
    }

    if(my_len==0){
        return NULL;
    }
    else {
        return (header->aux_meta->attrs);
    }
}


enum slow5_aux_type *slow5_get_aux_types(const slow5_hdr_t *header,uint64_t *len){
    u_int64_t my_len = (header->aux_meta == NULL) ? 0 : header->aux_meta->num;

    if(len){
        *len = my_len;
    }

    if(my_len==0){
        return NULL;
    }
    else {
        return (header->aux_meta->types);
    }
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

            SLOW5_ASSERT(slow5_hdr_add_attr(attr, header) != -3);
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
    if (header != NULL) {
        slow5_hdr_data_free(header);
        slow5_aux_meta_free(header->aux_meta);

        free(header);
    }
}


/************************************* slow5 header data *************************************/

//read slow5 data header from file
int slow5_hdr_data_init(FILE *fp, char *buf, size_t *cap, struct slow5_hdr *header, uint32_t *hdr_len) {

    int ret = 0;

    char *buf_orig = buf;
    uint32_t hdr_len_tmp = 0;

    kv_init(header->data.maps);
    kv_resize(khash_t(slow5_s2s) *, header->data.maps, header->num_read_groups);

    for (uint64_t i = 0; i < (uint64_t) header->num_read_groups; ++ i) {
        kv_A(header->data.maps, i) = kh_init(slow5_s2s);
        ++ header->data.maps.n;
    }

    khash_t(slow5_s) *data_attrs = kh_init(slow5_s);

    // Parse slow5 header data

    ssize_t buf_len;

    // Get first line of header data
    SLOW5_ASSERT((buf_len = getline(&buf, cap, fp)) != -1);
    buf[buf_len - 1] = '\0'; // Remove newline for later parsing
    if (hdr_len != NULL) {
        hdr_len_tmp += buf_len;
    }

    uint32_t num_data_attrs = 0;
    // While the column header hasn't been reached
    while (strncmp(buf, SLOW5_ASCII_TYPE_HEADER_MIN, strlen(SLOW5_ASCII_TYPE_HEADER_MIN)) != 0) {

        // Ensure prefix is there
        SLOW5_ASSERT(buf[0] == SLOW5_HEADER_DATA_PREFIX_CHAR);
        char *shift = buf + strlen(SLOW5_HEADER_DATA_PREFIX); // Remove prefix

        // Get the attribute name
        char *attr = strdup(slow5_strsep(&shift, SLOW5_SEP_COL));
        char *val;

        int ret;
        kh_put(slow5_s, data_attrs, attr, &ret);
        ++ num_data_attrs;
        SLOW5_ASSERT(!(ret == -1 || ret == 0));

        // Iterate through the values
        uint32_t i = 0;
        while ((val = slow5_strsep(&shift, SLOW5_SEP_COL)) != NULL && i <= header->num_read_groups - 1) {

            // Set key
            int absent;
            khint_t pos = kh_put(slow5_s2s, header->data.maps.a[i], attr, &absent);
            SLOW5_ASSERT(absent != -1);

            //if the value is ".", we store an empty string
            char *val_dup = strdup(val);
            if(strcmp(val_dup,".") == 0 ){
                val_dup[0]='\0';
            }

            // Set value
            kh_val(header->data.maps.a[i], pos) = val_dup;

            ++ i;
        }
        // Ensure that read group number of entries are read
        SLOW5_ASSERT(i == header->num_read_groups);

        // Get next line
        SLOW5_ASSERT((buf_len = getline(&buf, cap, fp)) != -1);
        buf[buf_len - 1] = '\0'; // Remove newline for later parsing
        if (hdr_len != NULL) {
            hdr_len_tmp += buf_len;
        }
    }

    SLOW5_ASSERT(*cap <= SLOW5_HEADER_DATA_BUF_INIT_CAP); // TESTING to see if getline has to realloc (if this fails often maybe put a larger buffer size)

    if (hdr_len != NULL) {
        *hdr_len = hdr_len_tmp;
    }

    if (buf != buf_orig) {
        free(buf);
    }

    if (ret == 0) {
        header->data.num_attrs = num_data_attrs;
        header->data.attrs = data_attrs;
    }

    return ret;
}

struct slow5_aux_meta *slow5_aux_meta_init_empty(void) {
    struct slow5_aux_meta *aux_meta = (struct slow5_aux_meta *) calloc(1, sizeof *aux_meta);

    aux_meta->cap = SLOW5_AUX_META_CAP_INIT;
    aux_meta->attrs = (char **) malloc(aux_meta->cap * sizeof *(aux_meta->attrs));
    aux_meta->types = (enum slow5_aux_type *) malloc(aux_meta->cap * sizeof *(aux_meta->types));
    aux_meta->sizes = (uint8_t *) malloc(aux_meta->cap * sizeof *(aux_meta->sizes));

    return aux_meta;
}

//reads the data type row in the header (from file) and populates slow5_aux_meta. note: the buffer should be the one used for slow5_hdr_init
struct slow5_aux_meta *slow5_aux_meta_init(FILE *fp, char *buf, size_t *cap, uint32_t *hdr_len) {

    struct slow5_aux_meta *aux_meta = NULL;

    ssize_t buf_len = strlen(buf);

    // Parse data types and deduce the sizes
    if (buf_len != strlen(SLOW5_ASCII_TYPE_HEADER_MIN)) {
        aux_meta = (struct slow5_aux_meta *) calloc(1, sizeof *aux_meta);
        char *shift = buf += strlen(SLOW5_ASCII_TYPE_HEADER_MIN);

        char *tok = slow5_strsep(&shift, SLOW5_SEP_COL);
        SLOW5_ASSERT(strcmp(tok, "") == 0);

        aux_meta->cap = SLOW5_AUX_META_CAP_INIT;
        aux_meta->types = (enum slow5_aux_type *) malloc(aux_meta->cap * sizeof *(aux_meta->types));
        aux_meta->sizes = (uint8_t *) malloc(aux_meta->cap * sizeof *(aux_meta->sizes));

        aux_meta->num = 0;
        while ((tok = slow5_strsep(&shift, SLOW5_SEP_COL)) != NULL) {
            int err;
            enum slow5_aux_type type = slow5_str_to_aux_type(tok, &err);

            if (err == -1) {
                free(aux_meta->types);
                free(aux_meta->sizes);
                free(aux_meta);
                return NULL;
            }

            aux_meta->types[aux_meta->num] = type;
            aux_meta->sizes[aux_meta->num] = SLOW5_AUX_TYPE_META[type].size;

            ++ aux_meta->num;
            if (aux_meta->num > aux_meta->cap) {
                aux_meta->cap = aux_meta->cap << 1; // TODO is this ok?
                aux_meta->types = (enum slow5_aux_type *) realloc(aux_meta->types, aux_meta->cap * sizeof *(aux_meta->types));
                aux_meta->sizes = (uint8_t *) realloc(aux_meta->sizes, aux_meta->cap * sizeof *(aux_meta->sizes));
            }
        }
    }

    // Get column names
    SLOW5_ASSERT((buf_len = getline(&buf, cap, fp)) != -1);
    if (hdr_len != NULL) {
        *hdr_len += buf_len;
    }
    buf[-- buf_len] = '\0'; // Remove newline for later parsing

    if (strncmp(buf, SLOW5_ASCII_COLUMN_HEADER_MIN, strlen(SLOW5_ASCII_COLUMN_HEADER_MIN)) != 0) {
        if (aux_meta != NULL) {
            free(aux_meta->types);
            free(aux_meta->sizes);
            free(aux_meta);
            return NULL;
        }
        return NULL;
    }

    // Parse auxiliary attributes
    if (buf_len != strlen(SLOW5_ASCII_COLUMN_HEADER_MIN)) {
        if (aux_meta == NULL) {
            return NULL;
        }
        char *shift = buf += strlen(SLOW5_ASCII_COLUMN_HEADER_MIN);

        char *tok = slow5_strsep(&shift, SLOW5_SEP_COL);
        SLOW5_ASSERT(strcmp(tok, "") == 0);

        aux_meta->attr_to_pos = kh_init(slow5_s2ui32);
        aux_meta->attrs = (char **) malloc(aux_meta->cap * sizeof *(aux_meta->attrs));

        for (uint64_t i = 0; i < aux_meta->num; ++ i) {
            if ((tok = slow5_strsep(&shift, SLOW5_SEP_COL)) == NULL) {
                for (uint64_t j = 0; j < i; ++ j) {
                    free(aux_meta->attrs[j]);
                }
                kh_destroy(slow5_s2ui32, aux_meta->attr_to_pos);
                free(aux_meta->attrs);
                free(aux_meta->types);
                free(aux_meta->sizes);
                free(aux_meta);
                return NULL;
            }

            aux_meta->attrs[i] = strdup(tok);

            int absent;
            khint_t pos = kh_put(slow5_s2ui32, aux_meta->attr_to_pos, aux_meta->attrs[i], &absent);
            if (absent == -1 || absent == -2) {
                for (uint64_t j = 0; j <= i; ++ j) {
                    free(aux_meta->attrs[j]);
                }
                kh_destroy(slow5_s2ui32, aux_meta->attr_to_pos);
                free(aux_meta->attrs);
                free(aux_meta->types);
                free(aux_meta->sizes);
                free(aux_meta);
                return NULL;
            }
            kh_value(aux_meta->attr_to_pos, pos) = i;
        }
        if ((tok = slow5_strsep(&shift, SLOW5_SEP_COL)) != NULL) {
            slow5_aux_meta_free(aux_meta);
            return NULL;
        }
    }

    return aux_meta;
}

// Return
// 0    success
// -1   null input
// -2   other failure
int slow5_aux_meta_add(struct slow5_aux_meta *aux_meta, const char *attr, enum slow5_aux_type type) {
    if (aux_meta == NULL || attr == NULL) {
        return -1;
    }

    if (aux_meta->attr_to_pos == NULL) {
        aux_meta->attr_to_pos = kh_init(slow5_s2ui32);
    }

    if (aux_meta->num == aux_meta->cap) {
        aux_meta->cap = aux_meta->cap << 1; // TODO is this ok?
        aux_meta->attrs = (char **) realloc(aux_meta->attrs, aux_meta->cap * sizeof *(aux_meta->attrs));
        aux_meta->types = (enum slow5_aux_type *) realloc(aux_meta->types, aux_meta->cap * sizeof *(aux_meta->types));
        aux_meta->sizes = (uint8_t *) realloc(aux_meta->sizes, aux_meta->cap * sizeof *(aux_meta->sizes));
    }

    aux_meta->attrs[aux_meta->num] = strdup(attr);

    int absent;
    khint_t pos = kh_put(slow5_s2ui32, aux_meta->attr_to_pos, aux_meta->attrs[aux_meta->num], &absent);
    if (absent == -1 || absent == -2) {
        free(aux_meta->attrs[aux_meta->num]);
        return -2;
    }
    kh_value(aux_meta->attr_to_pos, pos) = aux_meta->num;

    aux_meta->types[aux_meta->num] = type;
    aux_meta->sizes[aux_meta->num] = SLOW5_AUX_TYPE_META[type].size;

    ++ aux_meta->num;

    return 0;
}

void slow5_aux_meta_free(struct slow5_aux_meta *aux_meta) {
    if (aux_meta != NULL) {
        for (uint64_t i = 0; i < aux_meta->num; ++ i) {
            free(aux_meta->attrs[i]);
        }
        free(aux_meta->attrs);
        kh_destroy(slow5_s2ui32, aux_meta->attr_to_pos);
        free(aux_meta->types);
        free(aux_meta->sizes);
        free(aux_meta);
    }
}

void slow5_hdr_data_free(struct slow5_hdr *header) {

    if (header->data.attrs != NULL && header->data.maps.a != NULL) {

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

        // Free header data map
        for (size_t j = 0; j < kv_size(header->data.maps); ++ j) {
            kh_destroy(slow5_s2s, header->data.maps.a[j]);
        }

        kh_destroy(slow5_s, header->data.attrs);
        kv_destroy(header->data.maps);
    }
}


// slow5 record
/**
 * Get a read entry from a slow5 file corresponding to a read_id.
 *
 * Allocates memory for *read if it is NULL.
 * Otherwise, the data in *read is freed and overwritten.
 * slow5_rec_free() should always be called when finished with the structure.
 *
 * Require the slow5 index to be loaded using slow5_idx_load
 *
 * Return:
 *  >=0   the read was successfully found and stored
 *  <0   error code
 *
 * Errors:
 * SLOW5_ERR_NOTFOUND   read_id was not found in the index
 * SLOW5_ERR_ARG        read_id, read or s5p is NULL
 * SLOW5_ERR_IO         other error when reading the slow5 file
 * SLOW5_ERR_RECPARSE   parsing error
 * SLOW5_ERR_NOIDX      the index has not been loaded
 *
 * @param   read_id the read identifier
 * @param   read    address of a slow5_rec pointer
 * @param   s5p     slow5 file
 * @return  error code described above
 */
int slow5_get(const char *read_id, struct slow5_rec **read, struct slow5_file *s5p) {
    if (read_id == NULL || read == NULL || s5p == NULL) {
        SLOW5_WARNING("%s","read_id, read and s5p cannot be NULL.");
        return SLOW5_ERR_ARG;
    }

    int ret = 0;
    char *read_mem = NULL;
    ssize_t bytes_to_read = -1;

    // index must be loaded
    if (s5p->index == NULL) {
        // index not loaded
        SLOW5_WARNING("%s","SLOW5 index should have been loaded using slow5_idx_load() before calling slow5_get().");
        return SLOW5_ERR_NOIDX;
    }

    // Get index record
    struct slow5_rec_idx read_index;
    if (slow5_idx_get(s5p->index, read_id, &read_index) == -1) {
        // read_id not found in index
        return SLOW5_ERR_NOTFOUND;
    }

    if (s5p->format == SLOW5_FORMAT_ASCII) {

        // Malloc string to hold the read
        read_mem = (char *) malloc(read_index.size * sizeof *read_mem);

        // Read into the string
        // Don't read in newline for parsing
        ssize_t bytes_to_read = (read_index.size - 1) * sizeof *read_mem;
        if (pread(s5p->meta.fd, read_mem, bytes_to_read, read_index.offset) != bytes_to_read) {
            free(read_mem);
            // reading error
            SLOW5_WARNING("pread could not read %ld bytes as expected.",(long)bytes_to_read);
            return SLOW5_ERR_IO;
        }

        // Null terminate
        read_mem[read_index.size - 1] = '\0';

    } else if (s5p->format == SLOW5_FORMAT_BINARY) {

        // Read into the string and miss the preceding size
        size_t bytes_to_read_sizet;
        read_mem = (char *) slow5_pread_depress_multi(s5p->compress->method, s5p->meta.fd,
                                                read_index.size - sizeof (slow5_rec_size_t),
                                                read_index.offset + sizeof (slow5_rec_size_t),
                                                &bytes_to_read_sizet);
        bytes_to_read = bytes_to_read_sizet;
        if (read_mem == NULL) {
            SLOW5_WARNING("%s","slow5_pread_depress_multi failed.");
            // reading error
            return SLOW5_ERR_IO;
        }
        /*
        read_mem = (char *) malloc(read_size);
        pread(s5p->meta.fd, read_mem, read_size, read_index.offset + sizeof read_size);
        printf("printing read_mem comp:\n"); // TESTING
        fwrite(read_mem, read_size, 1, stdout); // TESTING
        size_t decomp_size = 0;
        void *read_decomp = slow5_ptr_depress(s5p->compress, read_mem, read_size, &decomp_size);
        printf("\nprinting read_mem decomp:\n"); // TESTING
        fwrite(read_decomp, read_size, 1, stdout); // TESTING
        */
    }

    if (*read == NULL) {
        // Allocate memory for read
        *read = (struct slow5_rec *) calloc(1, sizeof **read);
    } else {
        // Free previously allocated strings
        free((*read)->read_id);
        (*read)->read_id = NULL;
        if ((*read)->aux_map != NULL) {
            // Free previously allocated auxiliary data
            slow5_rec_aux_free((*read)->aux_map);
            (*read)->aux_map = NULL;
        }
    }

    if (slow5_rec_parse(read_mem, bytes_to_read, read_id, *read, s5p->format, s5p->header->aux_meta) == -1) {
        SLOW5_WARNING("%s","SLOW5 record parsing failed.");
        ret = SLOW5_ERR_RECPARSE;
    }
    free(read_mem);

    return ret;
}

// Return -1 on failure to parse
int slow5_rec_parse(char *read_mem, size_t read_size, const char *read_id, struct slow5_rec *read, enum slow5_fmt format, struct slow5_aux_meta *aux_meta) {
    int ret = 0;
    uint64_t prev_len_raw_signal = 0;

    if (format == SLOW5_FORMAT_ASCII) {

        char *tok;
        if ((tok = slow5_strsep(&read_mem, SLOW5_SEP_COL)) == NULL) {
            SLOW5_WARNING("%s","Error when parsing readID.");
            return -1;
        }

        uint8_t i = 0;
        bool more_to_parse = false;
        int err;
        do {
            switch (i) {
                case COL_read_id:
                    // Ensure line matches requested id
                    if (read_id != NULL) {
                        if (strcmp(tok, read_id) != 0) {
                            ret = -1;
                            SLOW5_WARNING("Requested read ID [%s] does not match the read ID in fetched record [%s].",read_id,tok);
                            break;
                        }
                    }
                    read->read_id_len = strlen(tok);
                    read->read_id = strdup(tok);
                    break;

                case COL_read_group:
                    read->read_group = slow5_ato_uint32(tok, &err);
                    if (err == -1) {
                        SLOW5_WARNING("%s","Error parsing read group.");
                        ret = -1;
                    }
                    break;

                case COL_digitisation:
                    read->digitisation = slow5_strtod_check(tok, &err);
                    if (err == -1) {
                        SLOW5_WARNING("%s","Error parsing digitisation.");
                        ret = -1;
                    }
                    break;

                case COL_offset:
                    read->offset = slow5_strtod_check(tok, &err);
                    if (err == -1) {
                        SLOW5_WARNING("%s","Error parsing offset.");
                        ret = -1;
                    }
                    break;

                case COL_range:
                    read->range = slow5_strtod_check(tok, &err);
                    if (err == -1) {
                        SLOW5_WARNING("%s","Error parsing range.");
                        ret = -1;
                    }
                    break;

                case COL_sampling_rate:
                    read->sampling_rate = slow5_strtod_check(tok, &err);
                    if (err == -1) {
                        SLOW5_WARNING("%s","Error parsing sampling rate.");
                        ret = -1;
                    }
                    break;

                case COL_len_raw_signal:
                    if (read->len_raw_signal != 0) {
                        prev_len_raw_signal = read->len_raw_signal;
                    }
                    read->len_raw_signal = slow5_ato_uint64(tok, &err);
                    if (err == -1) {
                        SLOW5_WARNING("%s","Error parsing raw signal length.");
                        ret = -1;
                    }
                    break;

                case COL_raw_signal: {
                    if (read->raw_signal == NULL) {
                        read->raw_signal = (int16_t *) malloc(read->len_raw_signal * sizeof *(read->raw_signal));
                    } else if (prev_len_raw_signal < read->len_raw_signal) {
                        read->raw_signal = (int16_t *) realloc(read->raw_signal, read->len_raw_signal * sizeof *(read->raw_signal));
                    }

                    char *signal_tok;
                    if ((signal_tok = slow5_strsep(&tok, SLOW5_SEP_ARRAY)) == NULL) {
                        // 0 signals
                        SLOW5_WARNING("%s","Error locating the raw signal.");
                        ret = -1;
                        break;
                    }

                    uint64_t j = 0;

                    // Parse raw signal
                    do {
                        (read->raw_signal)[j] = slow5_ato_int16(signal_tok, &err);
                        if (err == -1) {
                            SLOW5_WARNING("%s","Error parsing the raw signal.");
                            ret = -1;
                            break;
                        }
                        ++ j;
                    } while ((signal_tok = slow5_strsep(&tok, SLOW5_SEP_ARRAY)) != NULL);
                    if (ret != -1 && j != read->len_raw_signal) {
                        SLOW5_WARNING("%s","Raw signal is potentially truncated.");
                        ret = -1;
                    }

                } break;

                    // All columns parsed
                default:
                    more_to_parse = true;
                    break;

            }
            ++ i;

        } while (ret != -1 &&
                 !more_to_parse &&
                 (tok = slow5_strsep(&read_mem, SLOW5_SEP_COL)) != NULL);

        if (i < SLOW5_COLS_NUM) {
            // Not all main columns parsed
            SLOW5_WARNING("%s","All SLOW5 main columns were not parsed.");
            ret = -1;
        } else if (!more_to_parse && aux_meta != NULL){
            SLOW5_WARNING("%s","Header contained auxiliary feilds but the record does not.");
            ret = -1;
        } else if (more_to_parse) {
            // Main columns parsed and more to parse
            if (aux_meta == NULL) {
                SLOW5_WARNING("%s","Auxiliary fields are missing in header, but present in record.");
                ret = -1;
            } else {

                // TODO abstract to function (slow5_rec_aux_parse)
                khash_t(slow5_s2a) *aux_map = slow5_rec_aux_init();

                for (i = 0; i < aux_meta->num; ++ i) {
                    if (tok == NULL) {
                        SLOW5_WARNING("Auxiliary field [%s] is missing.",aux_meta->attrs[i]);
                        slow5_rec_aux_free(aux_map);
                        return -1;
                    }

                    uint64_t bytes = 0;
                    uint64_t len = 1;
                    uint8_t *data = NULL;
                    if (SLOW5_IS_PTR(aux_meta->types[i])) {
                        // Type is an array
                        if (aux_meta->types[i] == SLOW5_STRING) {
                            if(strcmp(tok,".")==0){
                                len = 0;
                                data = 0;
                            }
                            else{
                                len = strlen(tok);
                                data = (uint8_t *) malloc((len + 1) * aux_meta->sizes[i]);
                                memcpy(data, tok, (len + 1) * aux_meta->sizes[i]);
                            }
                        } else {
                            if(strcmp(tok,".")==0){
                                len = 0;
                                data = 0;
                            }
                            else{
                                // Split tok by SLOW5_SEP_ARRAY
                                char *tok_sep;
                                uint64_t array_cap = SLOW5_AUX_ARRAY_CAP_INIT;
                                uint64_t array_i = 0;
                                data = (uint8_t *) malloc(array_cap * aux_meta->sizes[i]);

                                while ((tok_sep = slow5_strsep(&tok, SLOW5_SEP_ARRAY)) != NULL) {
                                    // Storing comma-separated array
                                    // Dynamic array creation

                                    // Memcpy giving the primitive type not the array type (that's why minus SLOW5_CHAR)
                                    if (slow5_memcpy_type_from_str(data + (aux_meta->sizes[i] * array_i), tok_sep, SLOW5_TO_PRIM_TYPE(aux_meta->types[i])) == -1) {
                                        free(data);
                                        slow5_rec_aux_free(aux_map);
                                        SLOW5_WARNING("Auxiliary fields [%s] parsing failed.",aux_meta->attrs[i]);
                                        return -1;
                                    }

                                    ++ array_i;

                                    if (array_i >= array_cap) {
                                        array_cap = array_cap << 1;
                                        data = (uint8_t *) realloc(data, array_cap * aux_meta->sizes[i]);
                                    }
                                }

                                len = array_i;

                            }
                        }
                        bytes = len * aux_meta->sizes[i];

                    } else {
                        data = (uint8_t *) malloc(aux_meta->sizes[i]);
                        if(strcmp(tok,".")==0){
                            slow5_memcpy_type_from_null_str(data, aux_meta->types[i]);
                        }
                        else{
                            if (slow5_memcpy_type_from_str(data, tok, aux_meta->types[i]) == -1) {
                                free(data);
                                slow5_rec_aux_free(aux_map);
                                SLOW5_WARNING("Auxiliary fields [%s] parsing failed.",aux_meta->attrs[i]);
                                return -1;
                            }
                        }
                        bytes = aux_meta->sizes[i];
                    }

                    int absent;
                    khint_t pos = kh_put(slow5_s2a, aux_map, aux_meta->attrs[i], &absent);
                    if (absent == -1 || absent == -2) {
                        slow5_rec_aux_free(aux_map);
                        SLOW5_WARNING("Auxiliary fields [%s] is duplicated.",aux_meta->attrs[i]);
                        return -1;
                    }
                    struct slow5_rec_aux_data *aux_data = &kh_value(aux_map, pos);
                    aux_data->len = len;
                    aux_data->bytes = bytes;
                    aux_data->data = data;
                    aux_data->type = aux_meta->types[i];

                    tok = slow5_strsep(&read_mem, SLOW5_SEP_COL);
                }
                // Ensure line ends
                if (tok != NULL) {
                    kh_destroy(slow5_s2a, aux_map);
                    SLOW5_WARNING("%s","The parsing prematurely ended while some more data remaining");
                    return -1;
                } else {
                    read->aux_map = aux_map;
                }
            }
        }

    } else if (format == SLOW5_FORMAT_BINARY) {

        int64_t i = 0;
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
                    // Ensure line matches requested id
                    if (read_id != NULL) {
                        if (strcmp(read->read_id, read_id) != 0) {
                            SLOW5_WARNING("Requested read ID [%s] does not match the read ID in fetched record [%s].",read_id,read->read_id);
                            ret = -1;
                            break;
                        }
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

                case COL_raw_signal:
                    size = read->len_raw_signal * sizeof *read->raw_signal;
                    read->raw_signal = (int16_t *) realloc(read->raw_signal, size);
                    if(read->raw_signal == NULL){
                        SLOW5_ERROR("%s", "Raw signal returned was NULL");
                    }
                    memcpy(read->raw_signal, read_mem + offset, size);
                    offset += size;
                    break;

                    // All columns parsed
                default:
                    main_cols_parsed = true;
                    break;
            }

            ++ i;
        }

        if (offset > read_size){// Read too much
            SLOW5_WARNING("Malformed record. offset %ld, read_size %ld",(long)offset, (long)read_size);
            ret = -1;
        }
        else if(aux_meta == NULL && offset < read_size){ // More to read but no auxiliary meta
            SLOW5_WARNING("More to read but no auxiliary meta. offset %ld, read_size %ld",(long)offset, (long)read_size);
            ret = -1;
        }
        else if (aux_meta != NULL && offset == read_size) { // No more to read but auxiliary meta expected
            SLOW5_WARNING("No more to read but auxiliary meta expected. offset %ld, read_size %ld",(long)offset, (long)read_size);
            ret = -1;
        }
        else if (aux_meta != NULL) {
            // Parse auxiliary data

            khash_t(slow5_s2a) *aux_map = slow5_rec_aux_init();

            for (i = 0; i < aux_meta->num; ++ i) {
                if (offset >= read_size) {
                    slow5_rec_aux_free(aux_map);
                    SLOW5_WARNING("Parsing auxiliary field [%s] failed",aux_meta->attrs[i]);
                    return -1;
                }

                uint64_t len = 1;
                if (SLOW5_IS_PTR(aux_meta->types[i])) {
                    // Type is an array
                    size = sizeof len;
                    memcpy(&len, read_mem + offset, size);
                    offset += size;
                }

                uint8_t *data;
                uint64_t bytes = len * aux_meta->sizes[i];

                if (len==0){
                    data = NULL;
                }
                else{
                    if (aux_meta->types[i] == SLOW5_STRING) {
                        data = (uint8_t *) malloc(bytes + 1);
                        memcpy(data, read_mem + offset, bytes);
                        offset += bytes;
                        data[bytes] = '\0';
                    } else {
                        data = (uint8_t *) malloc(bytes);
                        memcpy(data, read_mem + offset, bytes);
                        offset += bytes;
                    }
                }

                int absent;
                khint_t pos = kh_put(slow5_s2a, aux_map, aux_meta->attrs[i], &absent);
                if (absent == -1 || absent == -2) {
                    slow5_rec_aux_free(aux_map);
                    SLOW5_WARNING("Auxiliary fields [%s] is duplicated.",aux_meta->attrs[i]);
                    return -1;
                }
                struct slow5_rec_aux_data *aux_data = &kh_value(aux_map, pos);
                aux_data->len = len;
                aux_data->bytes = bytes;
                aux_data->data = data;
                aux_data->type = aux_meta->types[i];
            }

            if (offset != read_size) {
                slow5_rec_aux_free(aux_map);
                SLOW5_WARNING("Malformed record. offset %ld, read_size %ld",(long)offset, (long)read_size);
                return -1;
            } else {
                read->aux_map = aux_map;
            }
        }
    }

    return ret;
}

static inline khash_t(slow5_s2a) *slow5_rec_aux_init(void) {
    khash_t(slow5_s2a) *aux_map = kh_init(slow5_s2a);
    return aux_map;
}

void slow5_rec_aux_free(khash_t(slow5_s2a) *aux_map) {
    if (aux_map != NULL) {
        for (khint_t i = kh_begin(aux_map); i != kh_end(aux_map); ++ i) {
            if (kh_exist(aux_map, i)) {
                kh_del(slow5_s2a, aux_map, i);
                struct slow5_rec_aux_data *aux_data = &kh_value(aux_map, i);
                free(aux_data->data);
            }
        }

        // Using aux_meta?
        // Then no independence between different slow5_file pointers
        // Hence, always close slow5 files after slow5 records?
        // Is this a restriction on the library or ok?
        /*
        for (uint16_t i = 0; i < aux_meta->num; ++ i) {
            char *attr = aux_meta->attrs[i];

            khint_t pos = kh_get(slow5_s2a, aux_map, attr);

            if (kh_exist(aux_map, pos)) {
                free((void *) kh_key(aux_map, pos)); // TODO avoid void *
                kh_del(slow5_s2a, aux_map, pos);
                struct slow5_rec_aux_data *aux_data = &kh_value(aux_map, pos);
                free(aux_data->data);
            }
        }
        */
        kh_destroy(slow5_s2a, aux_map);
    }
}

/**
 * Get the read entry under the current file pointer of a slow5 file.
 *
 * Allocates memory for *read if it is NULL.
 * Otherwise, the data in *read is freed and overwritten.
 * slow5_rec_free() should be called when finished with the structure.
 *
 * Return
 * >=0  the read was successfully found and stored
 * <0   error code
 * Error
 * SLOW5_ERR_EOF       EOF reached
 * SLOW5_ERR_ARG       read_id, read or s5p is NULL
 * SLOW5_ERR_IO      ither reading error when reading the slow5 file
 * SLOW5_ERR_RECPARSE    record parsing error
 *
 * @param   read    address of a slow5_rec pointer
 * @param   s5p     slow5 file
 * @return  error code described above
 */
int slow5_get_next(struct slow5_rec **read, struct slow5_file *s5p) {
    if (read == NULL || s5p == NULL) {
        return SLOW5_ERR_ARG;
    }

    int ret = 0;
    char *read_mem = NULL;

    if (s5p->format == SLOW5_FORMAT_ASCII) {
        size_t cap = 0;
        ssize_t read_len;
        if ((read_len = getline(&read_mem, &cap, s5p->fp)) == -1) {
            free(read_mem);
            if(feof(s5p->fp)){
                return SLOW5_ERR_EOF;
            }
            return SLOW5_ERR_IO;
        }
        read_mem[-- read_len] = '\0'; // Remove newline for parsing

        if (*read == NULL) {
            // Allocate memory for read
            *read = (struct slow5_rec *) calloc(1, sizeof **read);
        } else {
            // Free previously allocated read id
            free((*read)->read_id);
            (*read)->read_id = NULL;
            if ((*read)->aux_map != NULL) {
                // Free previously allocated auxiliary data
                slow5_rec_aux_free((*read)->aux_map);
                (*read)->aux_map = NULL;
            }
        }

        if (slow5_rec_parse(read_mem, read_len, NULL, *read, s5p->format, s5p->header->aux_meta) == -1) {
            ret = SLOW5_ERR_RECPARSE;
        }

        free(read_mem);

    } else if (s5p->format == SLOW5_FORMAT_BINARY) {
        if (*read == NULL) {
            // Allocate memory for read
            *read = (struct slow5_rec *) calloc(1, sizeof **read);
        } else {
            // Free previously allocated read id
            free((*read)->read_id);
            (*read)->read_id = NULL;
            if ((*read)->aux_map != NULL) {
                // Free previously allocated auxiliary data
                slow5_rec_aux_free((*read)->aux_map);
                (*read)->aux_map = NULL;
            }
        }

        slow5_rec_size_t record_size;
        if (fread(&record_size, sizeof record_size, 1, s5p->fp) != 1) {
           if(feof(s5p->fp)){
                return SLOW5_ERR_EOF;
            }
            return SLOW5_ERR_IO;
        }

        size_t size_decomp;
        char *rec_decomp = (char *) slow5_fread_depress(s5p->compress, record_size, s5p->fp, &size_decomp);
        if (rec_decomp == NULL) {
            return SLOW5_ERR_IO;
        }

        if (slow5_rec_parse(rec_decomp, size_decomp, NULL, *read, s5p->format, s5p->header->aux_meta) == -1) {
            ret = SLOW5_ERR_RECPARSE;
        }

        free(rec_decomp);
    }

    return ret;
}

// For non-array types
// Return
// -1   input invalid
// -2   attr not found
// -3   type is an array type
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

    if (read->aux_map == NULL) {
        read->aux_map = kh_init(slow5_s2a);
    }
    slow5_rec_set_aux_map(read->aux_map, attr, (uint8_t *) data, 1, aux_meta->sizes[i], aux_meta->types[i]);

    return 0;
}

// For array types
// Return
// -1   input invalid
// -2   attr not found
// -3   type is not an array type
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

    if (read->aux_map == NULL) {
        read->aux_map = kh_init(slow5_s2a);
    }
    slow5_rec_set_aux_map(read->aux_map, attr, (uint8_t *) data, len, aux_meta->sizes[i] * len, aux_meta->types[i]);

    return 0;
}

static inline void slow5_rec_set_aux_map(khash_t(slow5_s2a) *aux_map, const char *attr, const uint8_t *data, size_t len, uint64_t bytes, enum slow5_aux_type type) {
    khint_t pos = kh_get(slow5_s2a, aux_map, attr);
    struct slow5_rec_aux_data *aux_data;
    if (pos != kh_end(aux_map)) {
        aux_data = &kh_value(aux_map, pos);
    } else {
        int ret;
        pos = kh_put(slow5_s2a, aux_map, attr, &ret);
        SLOW5_ASSERT(ret != -1);
        aux_data = &kh_value(aux_map, pos);
    }
    aux_data->len = len;
    aux_data->bytes = bytes;
    aux_data->type = type;
    aux_data->data = (uint8_t *) malloc(bytes);
    memcpy(aux_data->data, data, bytes);
}

int8_t slow5_aux_get_int8(const struct slow5_rec *read, const char *field, int *err) {
    int8_t val = INT8_MAX;
    int tmp_err = -1;

    if (read != NULL && field != NULL && read->aux_map != NULL) {
        khint_t pos = kh_get(slow5_s2a, read->aux_map, field);
        if (pos != kh_end(read->aux_map)) {
            struct slow5_rec_aux_data aux_data = kh_value(read->aux_map, pos);
            if (aux_data.type == SLOW5_INT8_T) {
                val = *((int8_t*) aux_data.data);
                tmp_err = 0;
            }
        }
    }

    if (err != NULL) {
        *err = tmp_err;
    }
    return val;
}
int16_t slow5_aux_get_int16(const struct slow5_rec *read, const char *field, int *err) {
    int16_t val = INT16_MAX;
    int tmp_err = -1;

    if (read != NULL && field != NULL && read->aux_map != NULL) {
        khint_t pos = kh_get(slow5_s2a, read->aux_map, field);
        if (pos != kh_end(read->aux_map)) {
            struct slow5_rec_aux_data aux_data = kh_value(read->aux_map, pos);
            if (aux_data.type == SLOW5_INT16_T) {
                val = *((int16_t*) aux_data.data);
                tmp_err = 0;
            }
        }
    }

    if (err != NULL) {
        *err = tmp_err;
    }
    return val;
}
int32_t slow5_aux_get_int32(const struct slow5_rec *read, const char *field, int *err) {
    int32_t val = INT32_MAX;
    int tmp_err = -1;

    if (read != NULL && field != NULL && read->aux_map != NULL) {
        khint_t pos = kh_get(slow5_s2a, read->aux_map, field);
        if (pos != kh_end(read->aux_map)) {
            struct slow5_rec_aux_data aux_data = kh_value(read->aux_map, pos);
            if (aux_data.type == SLOW5_INT32_T) {
                val = *((int32_t*) aux_data.data);
                tmp_err = 0;
            }
        }
    }

    if (err != NULL) {
        *err = tmp_err;
    }
    return val;
}
int64_t slow5_aux_get_int64(const struct slow5_rec *read, const char *field, int *err) {
    int64_t val = INT64_MAX;
    int tmp_err = -1;

    if (read != NULL && field != NULL && read->aux_map != NULL) {
        khint_t pos = kh_get(slow5_s2a, read->aux_map, field);
        if (pos != kh_end(read->aux_map)) {
            struct slow5_rec_aux_data aux_data = kh_value(read->aux_map, pos);
            if (aux_data.type == SLOW5_INT64_T) {
                val = *((int64_t*) aux_data.data);
                tmp_err = 0;
            }
        }
    }

    if (err != NULL) {
        *err = tmp_err;
    }
    return val;
}
uint8_t slow5_aux_get_uint8(const struct slow5_rec *read, const char *field, int *err) {
    uint8_t val = UINT8_MAX;
    int tmp_err = -1;

    if (read != NULL && field != NULL && read->aux_map != NULL) {
        khint_t pos = kh_get(slow5_s2a, read->aux_map, field);
        if (pos != kh_end(read->aux_map)) {
            struct slow5_rec_aux_data aux_data = kh_value(read->aux_map, pos);
            if (aux_data.type == SLOW5_UINT8_T) {
                val = *((uint8_t*) aux_data.data);
                tmp_err = 0;
            }
        }
    }

    if (err != NULL) {
        *err = tmp_err;
    }
    return val;
}
uint16_t slow5_aux_get_uint16(const struct slow5_rec *read, const char *field, int *err) {
    uint16_t val = UINT16_MAX;
    int tmp_err = -1;

    if (read != NULL && field != NULL && read->aux_map != NULL) {
        khint_t pos = kh_get(slow5_s2a, read->aux_map, field);
        if (pos != kh_end(read->aux_map)) {
            struct slow5_rec_aux_data aux_data = kh_value(read->aux_map, pos);
            if (aux_data.type == SLOW5_UINT16_T) {
                val = *((uint16_t*) aux_data.data);
                tmp_err = 0;
            }
        }
    }

    if (err != NULL) {
        *err = tmp_err;
    }
    return val;
}
uint32_t slow5_aux_get_uint32(const struct slow5_rec *read, const char *field, int *err) {
    uint32_t val = UINT32_MAX;
    int tmp_err = -1;

    if (read != NULL && field != NULL && read->aux_map != NULL) {
        khint_t pos = kh_get(slow5_s2a, read->aux_map, field);
        if (pos != kh_end(read->aux_map)) {
            struct slow5_rec_aux_data aux_data = kh_value(read->aux_map, pos);
            if (aux_data.type == SLOW5_UINT32_T) {
                val = *((uint32_t*) aux_data.data);
                tmp_err = 0;
            }
        }
    }

    if (err != NULL) {
        *err = tmp_err;
    }
    return val;
}
uint64_t slow5_aux_get_uint64(const struct slow5_rec *read, const char *field, int *err) {
    uint64_t val = UINT64_MAX;
    int tmp_err = -1;

    if (read != NULL && field != NULL && read->aux_map != NULL) {
        khint_t pos = kh_get(slow5_s2a, read->aux_map, field);
        if (pos != kh_end(read->aux_map)) {
            struct slow5_rec_aux_data aux_data = kh_value(read->aux_map, pos);
            if (aux_data.type == SLOW5_UINT64_T) {
                val = *((uint64_t*) aux_data.data);
                tmp_err = 0;
            }
        }
    }

    if (err != NULL) {
        *err = tmp_err;
    }
    return val;
}
float slow5_aux_get_float(const struct slow5_rec *read, const char *field, int *err) {
    float val = FLT_MAX;
    int tmp_err = -1;

    if (read != NULL && field != NULL && read->aux_map != NULL) {
        khint_t pos = kh_get(slow5_s2a, read->aux_map, field);
        if (pos != kh_end(read->aux_map)) {
            struct slow5_rec_aux_data aux_data = kh_value(read->aux_map, pos);
            if (aux_data.type == SLOW5_FLOAT) {
                val = *((float*) aux_data.data);
                tmp_err = 0;
            }
        }
    }

    if (err != NULL) {
        *err = tmp_err;
    }
    return val;
}
double slow5_aux_get_double(const struct slow5_rec *read, const char *field, int *err) {
    double val = DBL_MAX;
    int tmp_err = -1;

    if (read != NULL && field != NULL && read->aux_map != NULL) {
        khint_t pos = kh_get(slow5_s2a, read->aux_map, field);
        if (pos != kh_end(read->aux_map)) {
            struct slow5_rec_aux_data aux_data = kh_value(read->aux_map, pos);
            if (aux_data.type == SLOW5_DOUBLE) {
                val = *((double*) aux_data.data);
                tmp_err = 0;
            }
        }
    }

    if (err != NULL) {
        *err = tmp_err;
    }
    return val;
}
char slow5_aux_get_char(const struct slow5_rec *read, const char *field, int *err) {
    char val = CHAR_MAX;
    int tmp_err = -1;

    if (read != NULL && field != NULL && read->aux_map != NULL) {
        khint_t pos = kh_get(slow5_s2a, read->aux_map, field);
        if (pos != kh_end(read->aux_map)) {
            struct slow5_rec_aux_data aux_data = kh_value(read->aux_map, pos);
            if (aux_data.type == SLOW5_CHAR) {
                val = *((char*) aux_data.data);
                tmp_err = 0;
            }
        }
    }

    if (err != NULL) {
        *err = tmp_err;
    }
    return val;
}
int8_t *slow5_aux_get_int8_array(const struct slow5_rec *read, const char *field, uint64_t *len, int *err) {
    int8_t *val = NULL;
    int tmp_err = -1;

    if (read != NULL && field != NULL && read->aux_map != NULL) {
        khint_t pos = kh_get(slow5_s2a, read->aux_map, field);
        if (pos != kh_end(read->aux_map)) {
            struct slow5_rec_aux_data aux_data = kh_value(read->aux_map, pos);
            if (aux_data.type == SLOW5_INT8_T_ARRAY) {
                val = (int8_t*) aux_data.data;
                *len = aux_data.len;
                tmp_err = 0;
            }
        }
    }

    if (err != NULL) {
        *err = tmp_err;
    }
    return val;
}
int16_t *slow5_aux_get_int16_array(const struct slow5_rec *read, const char *field, uint64_t *len, int *err) {
    int16_t *val = NULL;
    int tmp_err = -1;

    if (read != NULL && field != NULL && read->aux_map != NULL) {
        khint_t pos = kh_get(slow5_s2a, read->aux_map, field);
        if (pos != kh_end(read->aux_map)) {
            struct slow5_rec_aux_data aux_data = kh_value(read->aux_map, pos);
            if (aux_data.type == SLOW5_INT16_T_ARRAY) {
                val = (int16_t *) aux_data.data;
                *len = aux_data.len;
                tmp_err = 0;
            }
        }
    }

    if (err != NULL) {
        *err = tmp_err;
    }
    return val;
}
int32_t *slow5_aux_get_int32_array(const struct slow5_rec *read, const char *field, uint64_t *len, int *err) {
    int32_t *val = NULL;
    int tmp_err = -1;

    if (read != NULL && field != NULL && read->aux_map != NULL) {
        khint_t pos = kh_get(slow5_s2a, read->aux_map, field);
        if (pos != kh_end(read->aux_map)) {
            struct slow5_rec_aux_data aux_data = kh_value(read->aux_map, pos);
            if (aux_data.type == SLOW5_INT32_T_ARRAY) {
                val = (int32_t *) aux_data.data;
                *len = aux_data.len;
                tmp_err = 0;
            }
        }
    }

    if (err != NULL) {
        *err = tmp_err;
    }
    return val;
}
int64_t *slow5_aux_get_int64_array(const struct slow5_rec *read, const char *field, uint64_t *len, int *err) {
    int64_t *val = NULL;
    int tmp_err = -1;

    if (read != NULL && field != NULL && read->aux_map != NULL) {
        khint_t pos = kh_get(slow5_s2a, read->aux_map, field);
        if (pos != kh_end(read->aux_map)) {
            struct slow5_rec_aux_data aux_data = kh_value(read->aux_map, pos);
            if (aux_data.type == SLOW5_INT64_T_ARRAY) {
                val = (int64_t *) aux_data.data;
                *len = aux_data.len;
                tmp_err = 0;
            }
        }
    }

    if (err != NULL) {
        *err = tmp_err;
    }
    return val;
}
uint8_t *slow5_aux_get_uint8_array(const struct slow5_rec *read, const char *field, uint64_t *len, int *err) {
    uint8_t *val = NULL;
    int tmp_err = -1;

    if (read != NULL && field != NULL && read->aux_map != NULL) {
        khint_t pos = kh_get(slow5_s2a, read->aux_map, field);
        if (pos != kh_end(read->aux_map)) {
            struct slow5_rec_aux_data aux_data = kh_value(read->aux_map, pos);
            if (aux_data.type == SLOW5_UINT8_T_ARRAY) {
                val = (uint8_t *) aux_data.data;
                *len = aux_data.len;
                tmp_err = 0;
            }
        }
    }

    if (err != NULL) {
        *err = tmp_err;
    }
    return val;
}
uint16_t *slow5_aux_get_uint16_array(const struct slow5_rec *read, const char *field, uint64_t *len, int *err) {
    uint16_t *val = NULL;
    int tmp_err = -1;

    if (read != NULL && field != NULL && read->aux_map != NULL) {
        khint_t pos = kh_get(slow5_s2a, read->aux_map, field);
        if (pos != kh_end(read->aux_map)) {
            struct slow5_rec_aux_data aux_data = kh_value(read->aux_map, pos);
            if (aux_data.type == SLOW5_UINT16_T_ARRAY) {
                val = (uint16_t *) aux_data.data;
                *len = aux_data.len;
                tmp_err = 0;
            }
        }
    }

    if (err != NULL) {
        *err = tmp_err;
    }
    return val;
}
uint32_t *slow5_aux_get_uint32_array(const struct slow5_rec *read, const char *field, uint64_t *len, int *err) {
    uint32_t *val = NULL;
    int tmp_err = -1;

    if (read != NULL && field != NULL && read->aux_map != NULL) {
        khint_t pos = kh_get(slow5_s2a, read->aux_map, field);
        if (pos != kh_end(read->aux_map)) {
            struct slow5_rec_aux_data aux_data = kh_value(read->aux_map, pos);
            if (aux_data.type == SLOW5_UINT32_T_ARRAY) {
                val = (uint32_t *) aux_data.data;
                *len = aux_data.len;
                tmp_err = 0;
            }
        }
    }

    if (err != NULL) {
        *err = tmp_err;
    }
    return val;
}
uint64_t *slow5_aux_get_uint64_array(const struct slow5_rec *read, const char *field, uint64_t *len, int *err) {
    uint64_t *val = NULL;
    int tmp_err = -1;

    if (read != NULL && field != NULL && read->aux_map != NULL) {
        khint_t pos = kh_get(slow5_s2a, read->aux_map, field);
        if (pos != kh_end(read->aux_map)) {
            struct slow5_rec_aux_data aux_data = kh_value(read->aux_map, pos);
            if (aux_data.type == SLOW5_UINT64_T_ARRAY) {
                val = (uint64_t *) aux_data.data;
                *len = aux_data.len;
                tmp_err = 0;
            }
        }
    }

    if (err != NULL) {
        *err = tmp_err;
    }
    return val;
}
float *slow5_aux_get_float_array(const struct slow5_rec *read, const char *field, uint64_t *len, int *err) {
    float *val = NULL;
    int tmp_err = -1;

    if (read != NULL && field != NULL && read->aux_map != NULL) {
        khint_t pos = kh_get(slow5_s2a, read->aux_map, field);
        if (pos != kh_end(read->aux_map)) {
            struct slow5_rec_aux_data aux_data = kh_value(read->aux_map, pos);
            if (aux_data.type == SLOW5_FLOAT_ARRAY) {
                val = (float *) aux_data.data;
                *len = aux_data.len;
                tmp_err = 0;
            }
        }
    }

    if (err != NULL) {
        *err = tmp_err;
    }
    return val;
}
double *slow5_aux_get_double_array(const struct slow5_rec *read, const char *field, uint64_t *len, int *err) {
    double *val = NULL;
    int tmp_err = -1;

    if (read != NULL && field != NULL && read->aux_map != NULL) {
        khint_t pos = kh_get(slow5_s2a, read->aux_map, field);
        if (pos != kh_end(read->aux_map)) {
            struct slow5_rec_aux_data aux_data = kh_value(read->aux_map, pos);
            if (aux_data.type == SLOW5_DOUBLE_ARRAY) {
                val = (double *) aux_data.data;
                *len = aux_data.len;
                tmp_err = 0;
            }
        }
    }

    if (err != NULL) {
        *err = tmp_err;
    }
    return val;
}
char *slow5_aux_get_string(const struct slow5_rec *read, const char *field, uint64_t *len, int *err) {
    char *val = NULL;
    int tmp_err = -1;

    if (read != NULL && field != NULL && read->aux_map != NULL) {
        khint_t pos = kh_get(slow5_s2a, read->aux_map, field);
        if (pos != kh_end(read->aux_map)) {
            struct slow5_rec_aux_data aux_data = kh_value(read->aux_map, pos);
            if (aux_data.type == SLOW5_STRING) {
                val = (char *) aux_data.data;
                if (len != NULL) {
                    *len = aux_data.len;
                }
                tmp_err = 0;
            }
        }
    }

    if (err != NULL) {
        *err = tmp_err;
    }
    return val;
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

        char sig_buf[SLOW5_SIGNAL_BUF_FIXED_CAP];

        uint64_t i;
        for (i = 1; i < read->len_raw_signal; ++ i) {
            int sig_len = sprintf(sig_buf, SLOW5_FORMAT_STRING_RAW_SIGNAL SLOW5_SEP_ARRAY, read->raw_signal[i - 1]);

            memcpy(mem + curr_len, sig_buf, sig_len);
            curr_len += sig_len;
        }
        if (read->len_raw_signal > 0) {
            // Trailing signal
            int len_to_cp = sprintf(sig_buf, SLOW5_FORMAT_STRING_RAW_SIGNAL, read->raw_signal[i - 1]);
            memcpy(mem + curr_len, sig_buf, len_to_cp);
            curr_len += len_to_cp;
        }

        // Auxiliary fields
        size_t cap = max_len;
        if (read->aux_map != NULL && aux_meta != NULL) {
            for (uint16_t i = 0; i < aux_meta->num; ++ i) {

                // Realloc if necessary
                if (curr_len + 1 >= cap) { // +1 for '\t'
                    cap *= 2;
                    mem = (char *) realloc(mem, cap);
                }
                mem[curr_len ++] = '\t';

                khint_t pos = kh_get(slow5_s2a, read->aux_map, aux_meta->attrs[i]);
                if (pos != kh_end(read->aux_map)) {
                    size_t type_len;
                    struct slow5_rec_aux_data aux_data = kh_value(read->aux_map, pos);
                    char *type_str = slow5_data_to_str(aux_data.data, aux_data.type, aux_data.len, &type_len);

                    // Realloc if necessary
                    if (curr_len + type_len >= cap) {
                        cap *= 2;
                        mem = (char *) realloc(mem, cap);
                    }

                    memcpy(mem + curr_len, type_str, type_len);
                    curr_len += type_len;

                    free(type_str);
                }

            }
        }

        // Trailing newline
        // Realloc if necessary
        if (curr_len + 2 >= cap) { // +2 for '\n' and '\0'
            cap *= 2;
            mem = (char *) realloc(mem, cap);
        }
        strcpy(mem + curr_len, "\n"); // Copies null byte as well
        curr_len += 1;

    } else if (format == SLOW5_FORMAT_BINARY) {

        bool compress_to_free = false;
        if (compress == NULL) {
            compress = slow5_press_init(SLOW5_COMPRESS_NONE);
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
        memcpy(mem + curr_len, &read->len_raw_signal, sizeof read->len_raw_signal);
        curr_len += sizeof read->len_raw_signal;
        memcpy(mem + curr_len, read->raw_signal, read->len_raw_signal * sizeof *read->raw_signal);
        curr_len += read->len_raw_signal * sizeof *read->raw_signal;

        // Auxiliary fields
        if (read->aux_map != NULL && aux_meta != NULL) { // TODO error if one is NULL but not another
            for (uint16_t i = 0; i < aux_meta->num; ++ i) {
                khint_t pos = kh_get(slow5_s2a, read->aux_map, aux_meta->attrs[i]);
                if (pos != kh_end(read->aux_map)) {
                    struct slow5_rec_aux_data aux_data = kh_value(read->aux_map, pos);

                    if (SLOW5_IS_PTR(aux_meta->types[i])) {
                        // Realloc if necessary
                        if (curr_len + sizeof aux_data.len + aux_data.bytes >= cap) {
                            cap *= 2;
                            mem = (char *) realloc(mem, cap);
                        }

                        memcpy(mem + curr_len, &aux_data.len, sizeof aux_data.len);
                        curr_len += sizeof aux_data.len;

                    } else {
                        // Realloc if necessary
                        if (curr_len + aux_data.bytes >= cap) {
                            cap *= 2;
                            mem = (char *) realloc(mem, cap);
                        }
                    }
                    if(aux_data.len != 0 ){
                        memcpy(mem + curr_len, aux_data.data, aux_data.bytes);
                    }
                    curr_len += aux_data.bytes;
                }
            }
        }

        slow5_compress_footer_next(compress);
        slow5_rec_size_t record_size;

        size_t record_sizet;
        void *comp_mem = slow5_ptr_compress(compress, mem, curr_len, &record_sizet);
        record_size = record_sizet;
        free(mem);

        if (comp_mem != NULL) {
            uint8_t *comp_mem_full = (uint8_t *) malloc(sizeof record_size + record_size);
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
    if (s5p->index){
        return 0;
    }
    else{
        return -1;
    }
}

void slow5_idx_unload(struct slow5_file *s5p){
    slow5_idx_free(s5p->index);
    s5p->index = NULL;
    return;
}

/**
 * Print the binary end of file to a file pointer.
 *
 * On success, the number of bytes written is returned.
 * On error, -1 is returned.
 *
 * @param   fp      output file pointer
 * @return  number of bytes written, -1 on error
 */
ssize_t slow5_eof_fwrite(FILE *fp) {
    const char eof[] = SLOW5_BINARY_EOF;

    size_t n;
    if ((n = fwrite(eof, sizeof *eof, sizeof eof, fp)) != sizeof eof) {
        return -1;
    } else {
        return n;
    }
}


// slow5 extension parsing

enum slow5_fmt slow5_name_get_fmt(const char *name) {
    enum slow5_fmt format = SLOW5_FORMAT_UNKNOWN;

    if (name != NULL) {
        for (size_t i = 0; i < sizeof SLOW5_FORMAT_META / sizeof SLOW5_FORMAT_META[0]; ++ i) {
            const struct slow5_fmt_meta meta = SLOW5_FORMAT_META[i];
            if (strcmp(meta.name, name) == 0) {
                format = meta.format;
                break;
            }
        }
    }

    return format;
}

enum slow5_fmt slow5_path_get_fmt(const char *path) {
    enum slow5_fmt format = SLOW5_FORMAT_UNKNOWN;

    if (path != NULL) {
        int64_t i;
        for (i = strlen(path) - 1; i >= 0; -- i) {
            if (path[i] == '.') {
                const char *ext = path + i + 1;
                format = slow5_name_get_fmt(ext);
                break;
            }
        }
    }

    return format;
}

// Get the slow5 format name from the format
const char *slow5_fmt_get_name(enum slow5_fmt format) {
    const char *str = NULL;

    for (size_t i = 0; i < sizeof SLOW5_FORMAT_META / sizeof SLOW5_FORMAT_META[0]; ++ i) {
        const struct slow5_fmt_meta meta = SLOW5_FORMAT_META[i];
        if (meta.format == format) {
            str = meta.name;
            break;
        }
    }

    return str;
}

char *slow5_get_idx_path(const char *path) {
    size_t new_len = strlen(path) + strlen(SLOW5_INDEX_EXTENSION);
    char *str = (char *) malloc((new_len + 1) * sizeof *str); // +1 for '\0'
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

#define STR_TO_AUX_TYPE(str, len, raw_type, prim_do, array_do) \
    (SLOW5_IS_TYPE_TRUNC(str, raw_type)) { \
        if (len == sizeof (# raw_type) - 1) { \
            prim_do; \
        } else if (len == sizeof (# raw_type) && str[len - 1] == '*') { \
            array_do; \
        } else { \
            *err = -1; \
            return type; \
        } \
    }

enum slow5_aux_type slow5_str_to_aux_type(const char *str, int *err) {
    enum slow5_aux_type type = SLOW5_INT8_T;

    size_t len = strlen(str);

    if STR_TO_AUX_TYPE(str, len, int8_t, type = SLOW5_INT8_T, type = SLOW5_INT8_T_ARRAY)
    else if STR_TO_AUX_TYPE(str, len, int16_t, type = SLOW5_INT16_T, type = SLOW5_INT16_T_ARRAY)
    else if STR_TO_AUX_TYPE(str, len, int32_t, type = SLOW5_INT32_T, type = SLOW5_INT32_T_ARRAY)
    else if STR_TO_AUX_TYPE(str, len, int64_t, type = SLOW5_INT64_T, type = SLOW5_INT64_T_ARRAY)
    else if STR_TO_AUX_TYPE(str, len, uint8_t, type = SLOW5_UINT8_T, type = SLOW5_UINT8_T_ARRAY)
    else if STR_TO_AUX_TYPE(str, len, uint16_t, type = SLOW5_UINT16_T, type = SLOW5_UINT16_T_ARRAY)
    else if STR_TO_AUX_TYPE(str, len, uint32_t, type = SLOW5_UINT32_T, type = SLOW5_UINT32_T_ARRAY)
    else if STR_TO_AUX_TYPE(str, len, uint64_t, type = SLOW5_UINT64_T, type = SLOW5_UINT64_T_ARRAY)
    else if STR_TO_AUX_TYPE(str, len, float, type = SLOW5_FLOAT, type = SLOW5_FLOAT_ARRAY)
    else if STR_TO_AUX_TYPE(str, len, double, type = SLOW5_DOUBLE, type = SLOW5_DOUBLE_ARRAY)
    else if STR_TO_AUX_TYPE(str, len, char, type = SLOW5_CHAR, type = SLOW5_STRING)
    else {
        *err = -1;
        return type;
    }

    *err = 0;
    return type;
}

int slow5_memcpy_type_from_str(uint8_t *data, const char *value, enum slow5_aux_type type) {
    int err = -1;

    // TODO fix this is disgusting :(
    if (type == SLOW5_INT8_T) {
        int8_t value_conv = slow5_ato_int8(value, &err);
        if (err != -1) {
            memcpy(data, &value_conv, sizeof value_conv);
        }
    } else if (type == SLOW5_UINT8_T) {
        uint8_t value_conv = slow5_ato_uint8(value, &err);
        if (err != -1) {
            memcpy(data, &value_conv, sizeof value_conv);
        }
    } else if (type == SLOW5_INT16_T) {
        int16_t value_conv = slow5_ato_int16(value, &err);
        if (err != -1) {
            memcpy(data, &value_conv, sizeof value_conv);
        }
    } else if (type == SLOW5_UINT16_T) {
        uint16_t value_conv = slow5_ato_uint16(value, &err);
        if (err != -1) {
            memcpy(data, &value_conv, sizeof value_conv);
        }
    } else if (type == SLOW5_INT32_T) {
        int32_t value_conv = slow5_ato_int32(value, &err);
        if (err != -1) {
            memcpy(data, &value_conv, sizeof value_conv);
        }
    } else if (type == SLOW5_UINT32_T) {
        uint32_t value_conv = slow5_ato_uint32(value, &err);
        if (err != -1) {
            memcpy(data, &value_conv, sizeof value_conv);
        }
    } else if (type == SLOW5_INT64_T) {
        int64_t value_conv = slow5_ato_int64(value, &err);
        if (err != -1) {
            memcpy(data, &value_conv, sizeof value_conv);
        }
    } else if (type == SLOW5_UINT64_T) {
        uint64_t value_conv = slow5_ato_uint64(value, &err);
        if (err != -1) {
            memcpy(data, &value_conv, sizeof value_conv);
        }
    } else if (type == SLOW5_FLOAT) {
        float value_conv = slow5_strtof_check(value, &err);
        if (err != -1) {
            memcpy(data, &value_conv, sizeof value_conv);
        }
    } else if (type == SLOW5_DOUBLE) {
        double value_conv = slow5_strtod_check(value, &err);
        if (err != -1) {
            memcpy(data, &value_conv, sizeof value_conv);
        }
    } else if (type == SLOW5_CHAR) {
        char value_conv = value[0];
        if (err != -1) {
            memcpy(data, &value_conv, sizeof value_conv);
        }
    }

    return err;
}

//if a .
void slow5_memcpy_type_from_null_str(uint8_t *data, enum slow5_aux_type type) {

    // TODO fix this is disgusting :(
    if (type == SLOW5_INT8_T) {
        int8_t value_conv = SLOW5_INT8_T_NULL;
        memcpy(data, &value_conv, sizeof value_conv);
    } else if (type == SLOW5_UINT8_T) {
        uint8_t value_conv = SLOW5_UINT8_T_NULL;
        memcpy(data, &value_conv, sizeof value_conv);
    } else if (type == SLOW5_INT16_T) {
        int16_t value_conv = SLOW5_INT16_T_NULL;
        memcpy(data, &value_conv, sizeof value_conv);
    } else if (type == SLOW5_UINT16_T) {
        uint16_t value_conv = SLOW5_UINT16_T_NULL;
        memcpy(data, &value_conv, sizeof value_conv);
    } else if (type == SLOW5_INT32_T) {
        int32_t value_conv = SLOW5_INT32_T_NULL;
        memcpy(data, &value_conv, sizeof value_conv);
    } else if (type == SLOW5_UINT32_T) {
        uint32_t value_conv = SLOW5_UINT32_T_NULL;
        memcpy(data, &value_conv, sizeof value_conv);
    } else if (type == SLOW5_INT64_T) {
        int64_t value_conv = SLOW5_INT64_T_NULL;
        memcpy(data, &value_conv, sizeof value_conv);
    } else if (type == SLOW5_UINT64_T) {
        uint64_t value_conv = SLOW5_UINT64_T_NULL;
        memcpy(data, &value_conv, sizeof value_conv);
    } else if (type == SLOW5_FLOAT) {
        float value_conv = SLOW5_FLOAT_NULL;
        memcpy(data, &value_conv, sizeof value_conv);
    } else if (type == SLOW5_DOUBLE) {
        double value_conv = SLOW5_DOUBLE_NULL;
        memcpy(data, &value_conv, sizeof value_conv);
    } else if (type == SLOW5_CHAR) {
        char value_conv = SLOW5_CHAR_NULL;
        memcpy(data, &value_conv, sizeof value_conv);
    }

    return;
}



char *slow5_data_to_str(uint8_t *data, enum slow5_aux_type type, uint64_t len, size_t *str_len) {
    char *str = NULL;

    if (type == SLOW5_INT8_T) {
        if(*(int8_t *) data == SLOW5_INT8_T_NULL){
            *str_len = slow5_asprintf(&str, ".");
            fprintf(stderr,"|%s|\n",str);
        }
        else {
            *str_len = slow5_asprintf(&str, "%" PRId8, *(int8_t *) data);
        }
    } else if (type == SLOW5_UINT8_T) {
        if(*(uint8_t *) data == SLOW5_UINT8_T_NULL){
            *str_len = slow5_asprintf(&str, ".");
        }
        else {
            *str_len = slow5_asprintf(&str, "%" PRIu8, *(uint8_t *) data);
        }
    } else if (type == SLOW5_INT16_T) {
        if(*(int16_t *) data == SLOW5_INT16_T_NULL){
            *str_len = slow5_asprintf(&str, ".");
        }
        else {
            *str_len = slow5_asprintf(&str, "%" PRId16, *(int16_t *) data);
        }
    } else if (type == SLOW5_UINT16_T) {
        if(*(uint16_t *) data == SLOW5_UINT16_T_NULL){
            *str_len = slow5_asprintf(&str, ".");
        }
        else {
            *str_len = slow5_asprintf(&str, "%" PRIu16, *(uint16_t *) data);
        }
    } else if (type == SLOW5_INT32_T) {
        if(*(int32_t *) data == SLOW5_INT32_T_NULL){
            *str_len = slow5_asprintf(&str, ".");
        }
        else {
            *str_len = slow5_asprintf(&str, "%" PRId32, *(int32_t *) data);
        }
    } else if (type == SLOW5_UINT32_T) {
        if(*(uint32_t *) data == SLOW5_UINT32_T_NULL){
            *str_len = slow5_asprintf(&str, ".");
        }
        else {
            *str_len = slow5_asprintf(&str, "%" PRIu32, *(uint32_t *) data);
        }
    } else if (type == SLOW5_INT64_T) {
        if(*(int64_t *) data == SLOW5_INT64_T_NULL){
            *str_len = slow5_asprintf(&str, ".");
        }
        else {
            *str_len = slow5_asprintf(&str, "%" PRId64, *(int64_t *) data);
        }
    } else if (type == SLOW5_UINT64_T) {
        if(*(uint64_t *) data == SLOW5_UINT64_T_NULL){
            *str_len = slow5_asprintf(&str, ".");
        }
        else {
            *str_len = slow5_asprintf(&str, "%" PRIu64, *(uint64_t *) data);
        }
    } else if (type == SLOW5_FLOAT) {
        if(isnan(*(float *) data)){
            *str_len = slow5_asprintf(&str, ".");
        }
        else {
            str = slow5_float_to_str(*(float *) data, str_len);
        }
    } else if (type == SLOW5_DOUBLE) {
        if(isnan(*(double *) data)){
            *str_len = slow5_asprintf(&str, ".");
        }
        else {
            str = slow5_double_to_str(*(double *) data, str_len);
        }
    } else if (type == SLOW5_CHAR) {
        if(*(char *) data == SLOW5_CHAR_NULL){
            *str_len = slow5_asprintf(&str, ".");
        }
        else {
            *str_len = sizeof (char);
            str = (char *) malloc(sizeof (char));
            str[0] = *(char *) data;
        }
    } else if (type == SLOW5_STRING) {
        if(len==0){
            *str_len = slow5_asprintf(&str, ".");
        }
        else{
            str = strdup((char *) data);
            *str_len = strlen(str);
        }
    } else if (SLOW5_IS_PTR(type)) {

        if(len==0){
            *str_len = slow5_asprintf(&str, ".");
        }
        else{
            size_t str_cap = SLOW5_AUX_ARRAY_STR_CAP_INIT;
            str = (char *) malloc(str_cap * sizeof *str);

            size_t str_cur = 0;
            uint64_t i;
            for (i = 0; i < len - 1; ++ i) {
                size_t str_len_sep;
                char *str_sep = slow5_data_to_str(data + i * SLOW5_AUX_TYPE_META[type].size, SLOW5_TO_PRIM_TYPE(type), 1, &str_len_sep);

                if (str_cur + str_len_sep + 1 > str_cap) { // +1 for SLOW5_SEP_ARRAY_CHAR
                    // Realloc
                    str_cap = str_cap << 1;
                    str = (char *) realloc(str, str_cap * sizeof *str);
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
            }

            memcpy(str + str_cur, str_sep, str_len_sep);
            str[str_cur + str_len_sep] = '\0';

            *str_len = str_cur + str_len_sep;
            free(str_sep);
        }
    }

    return str;
}

void slow5_set_log_level(enum slow5_log_level_opt log_level){
    slow5_log_level = log_level;
}
void slow5_set_exit_condition(enum slow5_exit_condition_opt exit_condition){
    slow5_exit_condition = exit_condition;
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
struct SLOW5WriteConf *conf = slow5_wconf_init(SLOW5_FORMAT_BINARY, SLOW5_COMPRESS_GZIP);
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
