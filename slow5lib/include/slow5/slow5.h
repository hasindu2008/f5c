/**
 * @file slow5.h
 * @brief SLOW5 API
 * @author Sasha Jenner (jenner.sasha@gmail.com), Hasindu Gamaarachchi (hasindu@garvan.org.au), Hiruna Samarakoon
 * @date 27/02/2021
 */

/* IMPORTANT: The comments in this header file are NOT the API documentation
The API documentation is available at https://hasindu2008.github.io/slow5tools/
The comments here are possibly not upto date
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

// Header with slow5 file definitions
// TODO structure pack to min size
// TODO fix and add function descriptions
// TODO remove unnessary header inclusions

#ifndef SLOW5_H
#define SLOW5_H

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include "klib/khash.h"
#include "klib/kvec.h"
#include "slow5_press.h"
#include "slow5_defs.h"
#include "slow5_error.h"

#ifdef __cplusplus
extern "C" {
#endif

/**************************************************************************************************
 ***  struct definitions **************************************************************************
 **************************************************************************************************/

/**
* @enum slow5_fmt
* File formats we are dealing with
*/
enum slow5_fmt {
    SLOW5_FORMAT_UNKNOWN, ///< the format is unknown, usually the case before parsing the file extension
    SLOW5_FORMAT_ASCII,   ///< the format is ASCII SLOW5
    SLOW5_FORMAT_BINARY   ///< the format is binary SLOW5 (that is, BLOW5)
};

/**
* @struct slow5_fmt_meta
* SLOW5 file meta data
*/
struct slow5_fmt_meta {
    const char *name;       ///< path of the SLOW5 file
    enum slow5_fmt format;  ///< the format of the SLOW5 file
};

static const struct slow5_fmt_meta SLOW5_FORMAT_META[] = {
    { SLOW5_ASCII_NAME,   SLOW5_FORMAT_ASCII    },
    { SLOW5_BINARY_NAME,  SLOW5_FORMAT_BINARY   }
};

/*** SLOW5 Header *********************************************************************************/

/**
* @struct slow5_version
* @brief SLOW5 file version
*/
struct slow5_version {
    uint8_t major;  ///< major version
    uint8_t minor;  ///< minor version
    uint8_t patch;  ///< patch version
};

static const struct slow5_version SLOW5_VERSION_STRUCT = SLOW5_VERSION_ARRAY;

// SLOW5 auxiliary types
// DO NOT rearrange! See subtracting SLOW5_INT8_T_ARRAY in SLOW5_TO_PRIM_TYPE
// if adding more in future, primitive types must be added after SLOW5_CHAR and arrays after SLOW5_STRING
// both the primitive type and the array type must be simultaneously added
enum slow5_aux_type {
    SLOW5_INT8_T = 0,
    SLOW5_INT16_T,
    SLOW5_INT32_T,
    SLOW5_INT64_T,
    SLOW5_UINT8_T,
    SLOW5_UINT16_T,
    SLOW5_UINT32_T,
    SLOW5_UINT64_T,
    SLOW5_FLOAT,
    SLOW5_DOUBLE,
    SLOW5_CHAR,
    SLOW5_ENUM,

    SLOW5_INT8_T_ARRAY,
    SLOW5_INT16_T_ARRAY,
    SLOW5_INT32_T_ARRAY,
    SLOW5_INT64_T_ARRAY,
    SLOW5_UINT8_T_ARRAY,
    SLOW5_UINT16_T_ARRAY,
    SLOW5_UINT32_T_ARRAY,
    SLOW5_UINT64_T_ARRAY,
    SLOW5_FLOAT_ARRAY,
    SLOW5_DOUBLE_ARRAY,
    SLOW5_STRING,
    SLOW5_ENUM_ARRAY,
};

#define SLOW5_IS_PTR(type)              ((type) >= SLOW5_INT8_T_ARRAY)
#define SLOW5_TO_PRIM_TYPE(ptr_type)    ((enum slow5_aux_type) ((ptr_type) - SLOW5_INT8_T_ARRAY))

/* NULL (missing value) representation */
#define SLOW5_INT8_T_NULL     INT8_MAX
#define SLOW5_INT16_T_NULL    INT16_MAX
#define SLOW5_INT32_T_NULL    INT32_MAX
#define SLOW5_INT64_T_NULL    INT64_MAX
#define SLOW5_UINT8_T_NULL    UINT8_MAX
#define SLOW5_UINT16_T_NULL   UINT16_MAX
#define SLOW5_UINT32_T_NULL   UINT32_MAX
#define SLOW5_UINT64_T_NULL   UINT64_MAX
#define SLOW5_FLOAT_NULL      nanf("")
#define SLOW5_DOUBLE_NULL     nan("")
#define SLOW5_CHAR_NULL       0
#define SLOW5_ENUM_NULL       SLOW5_UINT8_T_NULL


// Type with corresponding size
struct slow5_aux_type_meta {
    enum slow5_aux_type type;
    uint8_t size;
    const char *type_str;
};

#define SLOW5_TO_STR(x) #x
#define SLOW5_AUX_TYPE_META_PRIM(aux_type, raw_type) { aux_type, sizeof (raw_type), SLOW5_TO_STR(raw_type) }
#define SLOW5_AUX_TYPE_META_ARRAY(aux_type, raw_type) { aux_type, sizeof (raw_type), SLOW5_TO_STR(raw_type) "*" }
//any modifications to slow5_aux_type should follow by appropriate modifications to this.
//the order should be identical to that in slow5_aux_type
static const struct slow5_aux_type_meta SLOW5_AUX_TYPE_META[] = {
    // Needs to be the same order as the enum definition
    SLOW5_AUX_TYPE_META_PRIM(   SLOW5_INT8_T,           int8_t      ),
    SLOW5_AUX_TYPE_META_PRIM(   SLOW5_INT16_T,          int16_t     ),
    SLOW5_AUX_TYPE_META_PRIM(   SLOW5_INT32_T,          int32_t     ),
    SLOW5_AUX_TYPE_META_PRIM(   SLOW5_INT64_T,          int64_t     ),
    SLOW5_AUX_TYPE_META_PRIM(   SLOW5_UINT8_T,          uint8_t     ),
    SLOW5_AUX_TYPE_META_PRIM(   SLOW5_UINT16_T,         uint16_t    ),
    SLOW5_AUX_TYPE_META_PRIM(   SLOW5_UINT32_T,         uint32_t    ),
    SLOW5_AUX_TYPE_META_PRIM(   SLOW5_UINT64_T,         uint64_t    ),
    SLOW5_AUX_TYPE_META_PRIM(   SLOW5_FLOAT,            float       ),
    SLOW5_AUX_TYPE_META_PRIM(   SLOW5_DOUBLE,           double      ),
    SLOW5_AUX_TYPE_META_PRIM(   SLOW5_CHAR,             char        ),
    { SLOW5_ENUM,       sizeof (uint8_t),   "enum"  },

    SLOW5_AUX_TYPE_META_ARRAY(  SLOW5_INT8_T_ARRAY,     int8_t      ),
    SLOW5_AUX_TYPE_META_ARRAY(  SLOW5_INT16_T_ARRAY,    int16_t     ),
    SLOW5_AUX_TYPE_META_ARRAY(  SLOW5_INT32_T_ARRAY,    int32_t     ),
    SLOW5_AUX_TYPE_META_ARRAY(  SLOW5_INT64_T_ARRAY,    int64_t     ),
    SLOW5_AUX_TYPE_META_ARRAY(  SLOW5_UINT8_T_ARRAY,    uint8_t     ),
    SLOW5_AUX_TYPE_META_ARRAY(  SLOW5_UINT16_T_ARRAY,   uint16_t    ),
    SLOW5_AUX_TYPE_META_ARRAY(  SLOW5_UINT32_T_ARRAY,   uint32_t    ),
    SLOW5_AUX_TYPE_META_ARRAY(  SLOW5_UINT64_T_ARRAY,   uint64_t    ),
    SLOW5_AUX_TYPE_META_ARRAY(  SLOW5_FLOAT_ARRAY,      float       ),
    SLOW5_AUX_TYPE_META_ARRAY(  SLOW5_DOUBLE_ARRAY,     double      ),
    SLOW5_AUX_TYPE_META_ARRAY(  SLOW5_STRING,           char        ),
    { SLOW5_ENUM_ARRAY, sizeof (uint8_t),   "enum*" },
};

// Auxiliary attribute to position map: attribute string -> index position
KHASH_MAP_INIT_STR(slow5_s2ui32, uint32_t)

/**
* @struct slow5_aux_meta
* SLOW5 auxiliary field metadata (information available in the SLOW5 header)
*/
typedef struct slow5_aux_meta {
    uint32_t num;                       ///< number of auxiliary fields
    size_t cap;                         ///< capacity of the arrays: attrs, types and sizes

    khash_t(slow5_s2ui32) *attr_to_pos; ///< hash table that maps field name string -> index position in the following arrays.
    char **attrs;                       ///< field names
    enum slow5_aux_type *types;         ///< field datatypes
    uint8_t *sizes;                     ///< field datatype sizes, for arrays this stores the size (in bytes) of the corresponding primitive type (TODO: this is probably redundant)

    /* the following are NULL if no auxiliary datatype is an enum type */
    char ***enum_labels;                ///< array of enum labels stored as strings
    uint8_t *enum_num_labels;           ///< array of number of enum labels
} slow5_aux_meta_t;

// Header data map: attribute string -> data string
KHASH_MAP_INIT_STR(slow5_s2s, char *)
// Header data attributes set
KHASH_SET_INIT_STR(slow5_s)

/**
* @struct slow5_hdr_data
* SLOW5 header data (constant attributes in FAST5 files)
*/
struct slow5_hdr_data {
    uint32_t num_attrs;                   ///< Number of data attributes
    khash_t(slow5_s) *attrs;              ///< Set of the data attribute keys for internal access(incase of multiple read groups, the union of keys from all read groups)
    kvec_t(khash_t(slow5_s2s) *) maps;    ///< Dynamic vector of hash maps (attribute key string -> attribute value string). Length of the vector is requal to  num_read_groups. Index in the vector corresponds to the read group number. The keys that are not relevant to a particular read group are not stored in this hash map.
};
typedef struct slow5_hdr_data slow5_hdr_data_t;

/**
* @struct slow5_hdr
* SLOW5 header
*/
struct slow5_hdr {
	struct slow5_version version;       ///< SLOW5 file version
    uint32_t num_read_groups;           ///< Number of read groups
    slow5_hdr_data_t data;              ///< Header data (constant fields in FAST5 files). Not to be directly accessed, use provided functions instead.
    slow5_aux_meta_t *aux_meta;         ///< Auxiliary field metadata. Not to be directly accessed, use provided functions instead.
};
typedef struct slow5_hdr slow5_hdr_t;

/*** SLOW5 record *********************************************************************************/

// SLOW5 primary record columns stored as an enum to keep  the order of the columns.
// TODO: make the first one is set to zero
enum slow5_cols {
    SLOW5_COLS_FOREACH(SLOW5_GENERATE_ENUM)
    SLOW5_COLS_NUM
};

/**
* @struct slow5_rec_aux_data
* SLOW5 auxiliary field data (represents a single SLOW5 auxiliary field of a particular read record)
*/
struct slow5_rec_aux_data {
    uint64_t len;       ///< number of elements in a array (if a primitive type this is always 1)
    uint64_t bytes;     ///< total number of bytes in data (currently, the allocated size which is equal to the amount of data in it)
    enum slow5_aux_type type; ///< data type of the auxiliary attribute
    uint8_t *data;      ///< raw data
};

// Header data map: auxiliary attribute string -> auxiliary data
KHASH_MAP_INIT_STR(slow5_s2a, struct slow5_rec_aux_data)

typedef uint64_t slow5_rec_size_t; //size of the whole record (in bytes)
typedef uint16_t slow5_rid_len_t;  //length of the read ID string (does not include null character)
typedef khash_t(slow5_s2a) slow5_aux_data_t;  //Auxiliary field name string -> auxiliary field data value

/**
* @struct slow5_rec
* SLOW5 record data struct (represents a single SLOW5 record)
*/
struct slow5_rec {
    slow5_rid_len_t read_id_len;        ///< length of the read ID string (does not include null character)
    SLOW5_COLS_FOREACH(SLOW5_GENERATE_STRUCT) ///< macro magic that generates the primary fields (see below)
                                        ///< char* read_id;
                                        ///< uint32_t read_group;
                                        ///< double digitisation;
                                        ///< double offset;
                                        ///< double range;
                                        ///< double sampling_rate;
                                        ///< uint64_t len_raw_signal;
                                        ///< int16_t* raw_signal;
    slow5_aux_data_t *aux_map;               ///< Auxiliary field name string -> auxiliary field data value. Not to be directly accessed, use provided functions instead.
};
typedef struct slow5_rec slow5_rec_t;

/*** SLOW5 file handler ***************************************************************************/

/**
* @struct slow5_file
* SLOW5 file meta data
*/
struct slow5_file_meta {
    const char *pathname;       ///< file path
    int fd;                     ///< file descriptor
    uint64_t start_rec_offset;  ///< offset (in bytes) of the first SLOW5 record (skipping the SLOW5 header; used for indexing)
    char *fread_buffer;         ///< buffer for fread
    const char *mode;           ///< file mode
};
typedef struct slow5_file_meta slow5_file_meta_t;

typedef struct slow5_idx slow5_idx_t;

/**
* @struct slow5_file
* SLOW5 file structure
*/
struct slow5_file {
    FILE *fp;                   ///< file pointer
    enum slow5_fmt format;      ///< whether SLOW5, BLOW5 etc...
    slow5_press_t *compress;    ///< compression related metadata
    slow5_hdr_t *header;        ///< SLOW5 header
    slow5_idx_t *index;         ///< SLOW5 index (NULL if not applicable)
    slow5_file_meta_t meta;     ///< file metadata
};
typedef struct slow5_file slow5_file_t;


/**************************************************************************************************
 ***  High-level API ******************************************************************************
 **************************************************************************************************/

/*
IMPORTANT: The high-level is stable
existing function prototypes must NOT be changed as such changes affects the backward compatibility
newer functions can be added while keeping the existing ones intact
*/

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
slow5_file_t *slow5_open(const char *pathname, const char *mode);

/**
 * Close a slow5 file and free its memory.
 * If the file had been opened for writing or appending, EOF marker will be written
 *
 * @param   s5p slow5 file structure
 * @return      same as fclose()
 */
int slow5_close(slow5_file_t *s5p);

/**
 * Create the index file for slow5 file.
 * Overwrites if already exists.
 *
 * @param   s5p slow5 file structure
 * @return  0 if successful,  <-1> on error
 */
int slow5_idx_create(slow5_file_t *s5p);

/**
 * Loads the index file for slow5 file.
 * Creates the index if not found.
 *
 * Return -1 on error,
 * 0 on success.
 *
 * @param   s5p slow5 file structure
 * @return  error codes described above
 */
int slow5_idx_load(slow5_file_t *s5p);

/**
 * Unloads an index associted to a slow5_file_t using slow5_idx_load and free the memory.
 *
 * @param   s5p slow5 file structure
 */
void slow5_idx_unload(slow5_file_t *s5p);

/**
 * Get a header data attribute for a particular read_group.
 *
 * Returns NULL if the attribute name doesn't exist
 * or the read group is out of range
 * or an input parameter is NULL.
 *
 * @param   attr        attribute name
 * @param   read_group  the read group
 * @param   header      slow5 header
 * @return  the attribute's value, or NULL on error
 */
char *slow5_hdr_get(const char *attr, uint32_t read_group, const slow5_hdr_t *header);

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
int slow5_get(const char *read_id, slow5_rec_t **read, slow5_file_t *s5p);

/**
 * Get the read entry under the current file pointer of a slow5 file.
 *
 * Allocates memory for *read if it is NULL.
 * Otherwise, the data in *read is freed and overwritten.
 * slow5_rec_free() should be called when finished with the structure.
 *
 * Return value:
 * >=0  the read was successfully found and stored
 * <0   error code
 *
 * Errors:
 * SLOW5_ERR_EOF        EOF reached
 * SLOW5_ERR_ARG        read_id, read or s5p is NULL
 * SLOW5_ERR_IO         other error when reading the slow5 file
 * SLOW5_ERR_RECPARSE   record parsing error
 *
 * @param   read    address of a slow5_rec_t pointer
 * @param   s5p     slow5 file
 * @return  error code described above
 */
int slow5_get_next(slow5_rec_t **read, slow5_file_t *s5p);


/**
 * Free a slow5 record.
 *
 * @param   read    address of a slow5_rec_t pointer
 */
void slow5_rec_free(slow5_rec_t *read);


/**
 * Get an auxiliary field in a SLOW5 record as an 8-bit signed integer.
 *
 * @param   read    address of a slow5_rec_t pointer
 * @param   field   auxiliary field name
 * @param   err     error code, 0 on success, <0 on failure and slow5_errno is set
 *                  SLOW5_ERR_ARG   if read or field is NULL
 *                  SLOW5_ERR_NOAUX if no auxiliary hash map for the record
 *                  SLOW5_ERR_NOFLD if the field was not found
 *                  SLOW5_ERR_TYPE if the desired return type does not match the field's type
 * @return  field data value or SLOW5_INT8_T_NULL on failure
 */
int8_t slow5_aux_get_int8(const slow5_rec_t *read, const char *field, int *err);
int16_t slow5_aux_get_int16(const slow5_rec_t *read, const char *field, int *err);
int32_t slow5_aux_get_int32(const slow5_rec_t *read, const char *field, int *err);
int64_t slow5_aux_get_int64(const slow5_rec_t *read, const char *field, int *err);
uint8_t slow5_aux_get_uint8(const slow5_rec_t *read, const char *field, int *err);
uint16_t slow5_aux_get_uint16(const slow5_rec_t *read, const char *field, int *err);
uint32_t slow5_aux_get_uint32(const slow5_rec_t *read, const char *field, int *err);
uint64_t slow5_aux_get_uint64(const slow5_rec_t *read, const char *field, int *err);
float slow5_aux_get_float(const slow5_rec_t *read, const char *field, int *err);
double slow5_aux_get_double(const slow5_rec_t *read, const char *field, int *err);
char slow5_aux_get_char(const slow5_rec_t *read, const char *field, int *err);
uint8_t slow5_aux_get_enum(const slow5_rec_t *read, const char *field, int *err);

/**
 * Get an auxiliary field in a SLOW5 record as an 8-bit signed integer array.
 *
 * @param   read    address of a slow5_rec_t pointer
 * @param   field   auxiliary field name
 * @param   len     number of data values in the returned array
 * @param   err     error code, 0 on success, <0 on failure and slow5_errno is set
 *                  SLOW5_ERR_ARG   if read or field is NULL
 *                  SLOW5_ERR_NOAUX if no auxiliary hash map for the record
 *                  SLOW5_ERR_NOFLD if the field was not found
 *                  SLOW5_ERR_TYPE if the desired return type does not match the field's type
 * @return  pointer to the array of data values or NULL on error
 */
int8_t *slow5_aux_get_int8_array(const slow5_rec_t *read, const char *field, uint64_t *len, int *err);
int16_t *slow5_aux_get_int16_array(const slow5_rec_t *read, const char *field, uint64_t *len, int *err);
int32_t *slow5_aux_get_int32_array(const slow5_rec_t *read, const char *field, uint64_t *len, int *err);
int64_t *slow5_aux_get_int64_array(const slow5_rec_t *read, const char *field, uint64_t *len, int *err);
uint8_t *slow5_aux_get_uint8_array(const slow5_rec_t *read, const char *field, uint64_t *len, int *err);
uint16_t *slow5_aux_get_uint16_array(const slow5_rec_t *read, const char *field, uint64_t *len, int *err);
uint32_t *slow5_aux_get_uint32_array(const slow5_rec_t *read, const char *field, uint64_t *len, int *err);
uint64_t *slow5_aux_get_uint64_array(const slow5_rec_t *read, const char *field, uint64_t *len, int *err);
float *slow5_aux_get_float_array(const slow5_rec_t *read, const char *field, uint64_t *len, int *err);
double *slow5_aux_get_double_array(const slow5_rec_t *read, const char *field, uint64_t *len, int *err);
char *slow5_aux_get_string(const slow5_rec_t *read, const char *field, uint64_t *len, int *err);
uint8_t *slow5_aux_get_enum_array(const slow5_rec_t *read, const char *field, uint64_t *len, int *err);

/****** Writing SLOW5 files ******.
 * This is just around the corner.
 * However, this is being procrastinated until someone requests. If anyone is interested please open a GitHub issue.
***/

/**
 * Adds a new header data attribute.
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
int slow5_hdr_add(const char *attr, slow5_hdr_t *header);

/**
 * Sets a header data attribute for a particular read_group.
 *
 * Doesn't take memory ownership of the value given.
 *
 * Returns -1 if the attribute name doesn't exist
 * or the read group is out of range
 * or an input parameter is NULL.
 * Returns 0 other.
 *
 * @param   attr        attribute name
 * @param   value       new attribute value
 * @param   read_group  the read group
 * @param   header      slow5 header
 * @return  0 on success, -1 on error
 */
int slow5_hdr_set(const char *attr, const char *value, uint32_t read_group, slow5_hdr_t *header);

/**
 * Adds an auxiliary field to a SLOW5 header.
 * Return
 *
 * 0    success
 * -1   null input
 * -2   other failure
 * -3   use slow5_aux_meta_add_enum instead if type is SLOW5_ENUM or SLOW5_ENUM_ARRAY
 * TODO this error checking is bad, reorder parameters
 * @param   field       field name
 * @param   type        slow5 data type
 * @param   header      pointer to the header
 * @return  0 on success, <0 on error
 */
int slow5_aux_add(const char *field, enum slow5_aux_type type, slow5_hdr_t *header);

/**
 * Writes the associated SLOW5 header to a SLOW5 file.
 *
 * On success, the number of bytes written is returned.
 * On error, -1 is returned.
 *
 * @param   s5p              slow5 file structure
 * @return  number of bytes written, -1 on error
 */
int slow5_hdr_write(slow5_file_t *s5p);

/**
 * Initialises an empty SLOW5 record.
 * To be freed with slow5_rec_free().
 *
 * @return  ptr to the record
 */
static inline slow5_rec_t *slow5_rec_init(void) {
    slow5_rec_t *read = (slow5_rec_t *) calloc(1, sizeof *read);

    return read;
}

// sets an auxiliary field (a primitive datatype) of a SLOW5 record
// For non-array types
// Return
// -1   input invalid
// -2   field not found
// -3   type is an array type
// -4   data is invalid (eg enum out of range)
// these two functions require the "field" to be lifetime scope. Can be fixed by pointing to header aux attribute list.
// https://github.com/hasindu2008/slow5lib/commit/63c81a25689b608275277b66300dde824d371bf7
int slow5_aux_set(slow5_rec_t *read, const char *field, const void *data, slow5_hdr_t *header);
//sets an auxiliary field (string datatype) of a SLOW5 record
int slow5_aux_set_string(slow5_rec_t *read, const char *field, const char *data, slow5_hdr_t *header);

/**
 * Writes a SLOW5 record to a SLOW5 file.
 *
 *
 * On success, the number of bytes written is returned.
 * On error, -1 is returned.
 *
 * @param   read             slow5_rec pointer
 * @param   s5p              slow5 file structure
 * @return  number of bytes written, -1 on error
 */
int slow5_write(slow5_rec_t *read, slow5_file_t *s5p);

//set compression for writing
//should be immediately done after opening a blow5 for writing (mode 'w')
int slow5_set_press(slow5_file_t *s5p, enum slow5_press_method rec_press, enum slow5_press_method sig_press);


/**************************************************************************************************
 ***  Low-level API *******************************************************************************
 **************************************************************************************************/

//set the log verbosity level. the log is printed to the standard error.
//sets a global variable, so not thread safe
void slow5_set_log_level(enum slow5_log_level_opt log_level);

//set the exit condition for slow5lib
//sets a global variable, so not thread safe
void slow5_set_exit_condition(enum slow5_exit_condition_opt exit_condition);


//get the list of hdr data keys in sorted order (only the returned pointer must be freed, not the ones inside - subjet to change)
//len is the numberof elements
//returns null if no attributes
const char **slow5_get_hdr_keys(const slow5_hdr_t *header,uint64_t *len);

//gets the list of read ids from the SLOW5 index
//the list of read is is a pointer and must not be freed by user
//*len will have the number of read ids
//NULL will be returned in case of error
char **slow5_get_rids(const slow5_file_t *s5p, uint64_t *len);

//get the pointer to auxilliary field names
char **slow5_get_aux_names(const slow5_hdr_t *header,uint64_t *len);
//get the pointer to auxilliary field types
enum slow5_aux_type *slow5_get_aux_types(const slow5_hdr_t *header,uint64_t *len);
/**
 * get the enum labels for a specific auxiliary field and set the number of labels in *n
 * return NULL on error and slow5_errno set to
 * SLOW5_ERR_ARG    if header, field NULL, n can be NULL
 * SLOW5_ERR_NOAUX  if auxiliary header is NULL
 * SLOW5_ERR_TYPE   if the enum labels or num_labels array is NULL, or the field type is not an enum type
 * SLOW5_ERR_NOFLD  if the auxiliary field was not found
 * SLOW5_ERR_MEM    memory allocation error
 */
char **slow5_get_aux_enum_labels(const slow5_hdr_t *header, const char *field, uint8_t *n);

int slow5_get_next_bytes(void **mem, size_t *bytes, slow5_file_t *s5p);

int slow_decode(void **mem, size_t *bytes, slow5_rec_t **read, slow5_file_t *s5p);

int slow5_encode(void **mem, size_t *bytes, slow5_rec_t *read, slow5_file_t *s5p);

int slow5_write_bytes(void *mem, size_t bytes, slow5_file_t *s5p);

/*
IMPORTANT: The following low-level API functions are not yet finalised or documented, until someone requests.
If anyone is interested, please open a GitHub issue, rather than trying to figure out from the code.
Function prototypes can be changed without notice or completely removed. So do NOT use these functions in your code.
these functions are used by slow5tools and pyslow5 - so any change to a function here means slow5tools and pyslow5 must be fixed.
*/

/**
 * Open a slow5 file of a specific format with a mode given it's pathname.
 *
 * Return NULL if pathname or mode is NULL, or if the format specified doesn't match the file.
 * slow5_open_with(pathname, mode, SLOW5_FORMAT_UNKNOWN) is equivalent to slow5_open(pathname, mode).
 *
 * Otherwise, return a slow5 file structure with the header parsed.
 * slow5_close() should be called when finished with the structure.
 *
 * TODO: same issues as in slow5_open are applicable to this
 *
 * @param   pathname    relative or absolute path to slow5 file
 * @param   mode        same mode as in fopen()
 * @param   format      format of the slow5 file
 * @return              slow5 file structure
 */
slow5_file_t *slow5_open_with(const char *pathname, const char *mode, enum slow5_fmt format);

/**
 * Add a new header data attribute.
 *
 * Returns -1 if an input parameter is NULL.
 * Returns -2 if the attribute already exists.
 * Returns 0 other.
 *
 * @param   attr    attribute name
 * @param   header  slow5 header
 * @return  0 on success, <0 on error as described above
 */
int slow5_hdr_add_attr(const char *attr, slow5_hdr_t *header);

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
// TODO check return type but should be large enough to return -1 and the largest read group
int64_t slow5_hdr_add_rg(slow5_hdr_t *header);


/**
 * Get the header in the specified format.
 *
 * Returns NULL if s5p is NULL
 * or format is SLOW5_FORMAT_UNKNOWN
 * or an internal error occurs.
 *
 * @param   header          slow5 header
 * @param   format          slow5 format to write the entry in
 * @param   comp            compression method
 * @param   written number of bytes written to the returned buffer
 * @return  malloced memory storing the slow5 header representation,
 *          to use free() on afterwards
 */
void *slow5_hdr_to_mem(slow5_hdr_t *header, enum slow5_fmt format, slow5_press_method_t comp, size_t *written);

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
int slow5_hdr_fwrite(FILE *fp, slow5_hdr_t *header, enum slow5_fmt format, slow5_press_method_t comp);
static inline int slow5_hdr_print(slow5_hdr_t *header, enum slow5_fmt format, slow5_press_method_t comp) {
    return slow5_hdr_fwrite(stdout, header, format, comp);
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
 * @param   written     number of bytes written to the returned buffer
 * @param   compress    compress structure
 * @return  malloced string to use free() on, NULL on error
 */
void *slow5_rec_to_mem(slow5_rec_t *read, slow5_aux_meta_t *aux_meta, enum slow5_fmt format, slow5_press_t *compress, size_t *n);

/**
 * Print a read entry in the specified format to a file pointer.
 *
 * On success, the number of bytes written is returned.
 * On error, -1 is returned.
 *
 * @param   fp      output file pointer
 * @param   read    slow5_rec pointer
 * @param   format  slow5 format to write entry in
 * @param   compress
 * @return  number of bytes written, -1 on error
 */
int slow5_rec_fwrite(FILE *fp, slow5_rec_t *read, slow5_aux_meta_t *aux_meta, enum slow5_fmt format, slow5_press_t *compress);
static inline int slow5_rec_print(slow5_rec_t *read, slow5_aux_meta_t *aux_meta, enum slow5_fmt format, slow5_press_t *compress) {
    return slow5_rec_fwrite(stdout, read, aux_meta, format, compress);
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
ssize_t slow5_eof_fwrite(FILE *fp);
static inline ssize_t slow5_eof_print(void) {
    return slow5_eof_fwrite(stdout);
}

// Return
// 0    success
// -1   input invalid
// -2   failure
int slow5_convert(slow5_file_t *from, FILE *to_fp, enum slow5_fmt to_format, slow5_press_method_t to_compress);



#ifdef __cplusplus
}
#endif


#endif
