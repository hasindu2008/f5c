/// @file htslib/ftidx.h
/// FASTT random access.
/*
   adpapted from htslib/faidx.h by Hasindu Gamaarachchi <hasindu@unsw.edu.au>
*/

/*

htslib/faidx.h:
 
   Copyright (C) 2008, 2009, 2013, 2014, 2016, 2017-2018 Genome Research Ltd.

   Author: Heng Li <lh3@sanger.ac.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#ifndef HTSLIB_FTIDX_H
#define HTSLIB_FTIDX_H

#include "htslib/hts_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

#define ASYNC 1

/** @file

  Index FASTT files and extract lines given a read ID.

  The fti file index columns for FASTT are:
    - read name
    - 0 unused
    - offset: number of bytes to skip to get to the start of the line
        from the beginning of the file
    - 0 unused
    - 0 unused

 */

struct __ftidx_t;
/// Opaque structure representing FASTT index
typedef struct __ftidx_t ftidx_t;

/// File format to be dealing with.
enum fti_format_options {
    FTI_NONE,
    FTI_FASTT, 
    FTI_FASTB //unused
};

/// Build index for a FASTT bgzip-compressed FASTT.
/**  @param  fn  FASTT file name
     @param  fnfti Name of .fti file to build.
     @param  fngzi Name of .gzi file to build (if fn is bgzip-compressed).
     @return     0 on success; or -1 on ftilure

If fnfti is NULL, ".fti" will be appended to fn to make the FTI file name.
If fngzi is NULL, ".gzi" will be appended to fn for the GZI file.  The GZI
file will only be built if fn is bgzip-compressed.
*/
int fti_build3(const char *fn, const char *fnfti, const char *fngzi) HTS_RESULT_USED;

/// Build index for a FASTT or bgzip-compressed FASTT file.
/** @param  fn  FASTT file name
    @return     0 on success; or -1 on ftilure

File "fn.fti" will be generated.  This function is equivalent to
fti_build3(fn, NULL, NULL);
*/
int fti_build(const char *fn) HTS_RESULT_USED;

/// Destroy a ftidx_t struct
void fti_destroy(ftidx_t *fti);

enum fti_load_options {
    FTI_CREATE = 0x01,
};

/// Load FASTT indexes.
/** @param  fn  File name of the FASTT file (can be compressed with bgzip).
    @param  fnfti File name of the FASTT index.
    @param  fngzi File name of the bgzip index.
    @param  flags Option flags to control index file caching and creation.
    @return Pointer to a ftidx_t struct on success, NULL on ftilure.

If fnfti is NULL, ".fti" will be appended to fn to make the FTI file name.
If fngzi is NULL, ".gzi" will be appended to fn for the bgzip index name.
The bgzip index is only needed if fn is compressed.

If (flags & FTI_CREATE) is true, the index files will be built using
fti_build3() if they are not already present.
*/
ftidx_t *fti_load3(const char *fn, const char *fnfti, const char *fngzi,
                   int flags);

/// Load index from "fn.fti".
/** @param  fn  File name of the FASTT file
    @return Pointer to a ftidx_t struct on success, NULL on ftilure.

This function is equivalent to fti_load3(fn, NULL, NULL, FTI_CREATE|FTI_CACHE);
*/
ftidx_t *fti_load(const char *fn);

/// Load FASTT indexes.
/** @param  fn  File name of the FASTT file (can be compressed with bgzip).
    @param  fnfti File name of the FASTT index.
    @param  fngzi File name of the bgzip index.
    @param  flags Option flags to control index file caching and creation.
    @param  format FASTT file format
    @return Pointer to a ftidx_t struct on success, NULL on ftilure.

If fnfti is NULL, ".fti" will be appended to fn to make the FTI file name.
If fngzi is NULL, ".gzi" will be appended to fn for the bgzip index name.
The bgzip index is only needed if fn is compressed.

If (flags & FTI_CREATE) is true, the index files will be built using
fti_build3() if they are not already present.
*/
ftidx_t *fti_load3_format(const char *fn, const char *fnfti, const char *fngzi,
                   int flags, enum fti_format_options format);

/// Load index from "fn.fti".
/** @param  fn  File name of the FASTT file
    @param  format FASTT file format
    @return Pointer to a ftidx_t struct on success, NULL on ftilure.

This function is equivalent to fti_load3_format(fn, NULL, NULL, FTI_CREATE|FTI_CACHE, format);
*/
ftidx_t *fti_load_format(const char *fn, enum fti_format_options format);

/// Fetch the sequence in a region
/** @param  fti  Pointer to the ftidx_t struct
    @param  reg  Region in the format "chr2:20,000-30,000"
    @param  len  Length of the region; -2 if seq not present, -1 general error
    @return      Pointer to the sequence; `NULL` on ftilure

The returned sequence is allocated by `malloc()` family and should be destroyed
by end users by calling `free()` on it.
*/
char *fti_fetch(const ftidx_t *fti, const char *reg, int *len);
#ifdef ASYNC
    char *fti_fetch_async(const ftidx_t *fti, const char *reg, int *len,struct aiocb *aiocb);
#endif

/// Fetch the quality string for a region for FASTQ files
/** @param  fti  Pointer to the ftidx_t struct
    @param  reg  Region in the format "chr2:20,000-30,000"
    @param  len  Length of the region; -2 if seq not present, -1 general error
    @return      Pointer to the quality string; null on ftilure

The returned quality string is allocated by `malloc()` family and should be destroyed
by end users by calling `free()` on it.
*/
char *fti_fetchqual(const ftidx_t *fti, const char *reg, int *len);

/// Fetch the number of sequences
/** @param  fti  Pointer to the ftidx_t struct
    @return      The number of sequences
*/
int ftidx_fetch_nseq(const ftidx_t *fti) HTS_DEPRECATED("Please use ftidx_nseq instead");

/// Fetch the sequence in a region
/** @param  fti  Pointer to the ftidx_t struct
    @param  c_name Region name
    @param  p_beg_i  Beginning position number (zero-based)
    @param  p_end_i  End position number (zero-based)
    @param  len  Length of the region; -2 if c_name not present, -1 general error
    @return      Pointer to the sequence; null on ftilure

The returned sequence is allocated by `malloc()` family and should be destroyed
by end users by calling `free()` on it.
*/
char *ftidx_fetch_seq(const ftidx_t *fti, const char *c_name, int p_beg_i, int p_end_i, int *len);

/// Fetch the quality string in a region for FASTQ files
/** @param  fti  Pointer to the ftidx_t struct
    @param  c_name Region name
    @param  p_beg_i  Beginning position number (zero-based)
    @param  p_end_i  End position number (zero-based)
    @param  len  Length of the region; -2 if c_name not present, -1 general error
    @return      Pointer to the sequence; null on ftilure

The returned sequence is allocated by `malloc()` family and should be destroyed
by end users by calling `free()` on it.
*/
char *ftidx_fetch_qual(const ftidx_t *fti, const char *c_name, int p_beg_i, int p_end_i, int *len);

/// Query if sequence is present
/**   @param  fti  Pointer to the ftidx_t struct
      @param  seq  Sequence name
      @return      1 if present or 0 if absent
*/
int ftidx_has_seq(const ftidx_t *fti, const char *seq);

/// Return number of sequences in fti index
int ftidx_nseq(const ftidx_t *fti);

/// Return name of i-th sequence
const char *ftidx_iseq(const ftidx_t *fti, int i);

/// Return sequence length, -1 if not present
int ftidx_seq_len(const ftidx_t *fti, const char *seq);

#ifdef __cplusplus
}
#endif

#endif
