/**
 * @file slow5_defs.h
 * @brief SLOW5 macro definitions
 * @author Sasha Jenner (jenner.sasha@gmail.com)
 * @date 27/02/2021
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

#ifndef SLOW5_DEFS_H
#define SLOW5_DEFS_H

#include <stdint.h>
#include <inttypes.h>

#ifdef __cplusplus
extern "C" {
#endif

/* This is for internal use only - do not use any of the following directly unless they are in the API documentation
The API documentation is available at https://hasindu2008.github.io/slow5tools/
*/

// library version
#define SLOW5_LIB_VERSION "0.6.0"

// maximum file version supported by this library - independent of slow5 library version above
// if updating change all 4 below
#define SLOW5_VERSION_MAJOR (0)
#define SLOW5_VERSION_MINOR (2)
#define SLOW5_VERSION_PATCH (0)
#define SLOW5_VERSION_STRING "0.2.0"

// file version helpers
#define SLOW5_VERSION_STRING_FORMAT \
    "%" PRIu8 SLOW5_HDR_FILE_VERSION_SEP \
    "%" PRIu8 SLOW5_HDR_FILE_VERSION_SEP \
    "%" PRIu8
#define SLOW5_VERSION_ARRAY { .major = SLOW5_VERSION_MAJOR, .minor = SLOW5_VERSION_MINOR, .patch = SLOW5_VERSION_PATCH }

// SLOW5 format specs
#define SLOW5_HDR_PREFIX            "#"
#define SLOW5_HDR_DATA_PREFIX       "@"
#define SLOW5_HDR_DATA_PREFIX_CHAR  '@'
#define SLOW5_COLUMN_HDR_PREFIX     "#"
#define SLOW5_SEP_COL               "\t"
#define SLOW5_SEP_COL_CHAR          '\t'
#define SLOW5_SEP_COL_NAME          "tab"
#define SLOW5_SEP_ARRAY             ","
#define SLOW5_SEP_ARRAY_CHAR        ','
#define SLOW5_HDR_FILE_VERSION      "slow5_version"
#define SLOW5_HDR_FILE_VERSION_SEP  "."
#define SLOW5_HDR_NUM_GROUPS        "num_read_groups"
#define SLOW5_HDR_NUM_GROUPS_INIT   (1)
#define SLOW5_HDR_ENUM_LABELS_BEGIN '{'
#define SLOW5_HDR_ENUM_LABELS_END   '}'

// Order, format string and type of main SLOW5 columns
// NOTE if this is changed, also edit:
//      slow5.c:slow5_rec_to_str SLOW5_FORMAT_BINARY
//      slow5idx_clean.c
#define SLOW5_COLS(col, end)                    \
    col(char*,      "%s",       read_id)        /* A malloced string */ \
    col(uint32_t,   "%" PRIu32, read_group)     \
    col(double,     "%s",       digitisation)   \
    col(double,     "%s",       offset)         \
    col(double,     "%s",       range)          \
    col(double,     "%s",       sampling_rate)  \
    col(uint64_t,   "%" PRIu64, len_raw_signal) \
    end(int16_t*,   ,           raw_signal) // Use end() for last column
// Format string of raw signal
#define SLOW5_FORMAT_STRING_RAW_SIGNAL "%" PRId16

// Apply the same function to each column including the last one
#define SLOW5_COLS_FOREACH(foo) SLOW5_COLS(foo, foo)

#define SLOW5_GENERATE_STRUCT(type, fmt, name)              type name;
#define SLOW5_GENERATE_ENUM(type, fmt, name)                COL_ ## name,
#define SLOW5_GENERATE_NAME_STRING(type, fmt, name)         #name
#define SLOW5_GENERATE_NAME_STRING_SEP(type, fmt, name)     SLOW5_GENERATE_NAME_STRING(type, fmt, name) SLOW5_SEP_COL
#define SLOW5_GENERATE_TYPE_STRING(type, fmt, name)         #type
#define SLOW5_GENERATE_TYPE_STRING_SEP(type, fmt, name)     SLOW5_GENERATE_TYPE_STRING(type, fmt, name) SLOW5_SEP_COL
#define SLOW5_GENERATE_FORMAT_STRING(type, fmt, name)       fmt
#define SLOW5_GENERATE_FORMAT_STRING_SEP(type, fmt, name)   SLOW5_GENERATE_FORMAT_STRING(type, fmt, name) SLOW5_SEP_COL
#define SLOW5_GENERATE_NULL(type, fmt, name)

// More SLOW5 specs
#define SLOW5_HDR_ID(header_name)           SLOW5_HDR_PREFIX header_name
#define SLOW5_HDR_FILE_VERSION_ID           SLOW5_HDR_ID(SLOW5_HDR_FILE_VERSION)
#define SLOW5_HDR_NUM_GROUPS_ID             SLOW5_HDR_ID(SLOW5_HDR_NUM_GROUPS)

#define SLOW5_HDR_ENTRY(header_name, data)  SLOW5_HDR_ID(header_name) SLOW5_SEP_COL data "\n"

// ASCII SLOW5 specs
#define SLOW5_ASCII_NAME                    "slow5"
#define SLOW5_ASCII_EXTENSION               "." SLOW5_ASCII_NAME
#define SLOW5_ASCII_NUM_GROUPS_FORMAT       "%" PRIu32
#define SLOW5_ASCII_ENTRY_VERSION           SLOW5_HDR_ENTRY(SLOW5_HDR_FILE_VERSION, SLOW5_VERSION_STRING)
#define SLOW5_ASCII_ENTRY_VERSION_FORMAT    SLOW5_HDR_ENTRY(SLOW5_HDR_FILE_VERSION, SLOW5_VERSION_STRING_FORMAT)
#define SLOW5_ASCII_ENTRY_NUM_GROUPS_FORMAT SLOW5_HDR_ENTRY(SLOW5_HDR_NUM_GROUPS, SLOW5_ASCII_NUM_GROUPS_FORMAT)
#define SLOW5_ASCII_HDR_FORMAT              SLOW5_ASCII_ENTRY_VERSION_FORMAT SLOW5_ASCII_ENTRY_NUM_GROUPS_FORMAT
#define SLOW5_ASCII_TYPE_HDR_MIN            SLOW5_COLUMN_HDR_PREFIX SLOW5_COLS(SLOW5_GENERATE_TYPE_STRING_SEP, SLOW5_GENERATE_TYPE_STRING)
#define SLOW5_ASCII_COLUMN_HDR_MIN          SLOW5_COLUMN_HDR_PREFIX SLOW5_COLS(SLOW5_GENERATE_NAME_STRING_SEP, SLOW5_GENERATE_NAME_STRING)
#define SLOW5_ASCII_MISSING                 "."
#define SLOW5_ASCII_MISSING_CHAR            '.'

// Binary SLOW5 specs
#define SLOW5_BINARY_NAME                   "blow5"
#define SLOW5_BINARY_EXTENSION              "." SLOW5_BINARY_NAME
#define SLOW5_BINARY_MAGIC_NUMBER           { 'B', 'L', 'O', 'W', '5', '\1' }
#define SLOW5_BINARY_EOF                    { '5', 'W', 'O', 'L', 'B' }
#define SLOW5_BINARY_HDR_SIZE_OFFSET        (64L)

// error codes
#define SLOW5_ERR_EOF           (-1)    // EOF reached
#define SLOW5_ERR_ARG           (-2)    // bad argument
#define SLOW5_ERR_TRUNC         (-3)    // file truncated, or size of header/record in blow5 differs to actual size read
#define SLOW5_ERR_RECPARSE      (-4)    // record parsing error
#define SLOW5_ERR_IO            (-5)    // other file I/O error (check errno for details)
#define SLOW5_ERR_NOIDX         (-6)    // index not loaded
#define SLOW5_ERR_NOTFOUND      (-7)    // read id not found
#define SLOW5_ERR_OTH           (-8)    // other error (big endian, internal error, etc.)
#define SLOW5_ERR_UNK           (-9)    // file format unknown
#define SLOW5_ERR_MEM           (-10)   // memory (re)allocation error
#define SLOW5_ERR_NOAUX         (-11)   // no auxiliary map
#define SLOW5_ERR_NOFLD         (-12)   // field not found
#define SLOW5_ERR_PRESS         (-13)   // (de)compression failure
#define SLOW5_ERR_MAGIC         (-14)   // magic number invalid
#define SLOW5_ERR_VERSION       (-15)   // version incompatible
#define SLOW5_ERR_HDRPARSE      (-16)   // header parsing error
#define SLOW5_ERR_TYPE          (-17)   // error relating to slow5 data type


#ifdef __cplusplus
}
#endif

#endif
