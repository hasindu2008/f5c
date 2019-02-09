/// @file htslib/ftidx.h
/// FASTA random access.
/*
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

/*  hts_defs.h -- Miscellaneous definitions.

    Copyright (C) 2013-2015,2017 Genome Research Ltd.

    Author: John Marshall <jm18@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#ifndef HTSLIB_HTS_DEFS_H
#define HTSLIB_HTS_DEFS_H

#ifdef __clang__
#ifdef __has_attribute
#define HTS_COMPILER_HAS(attribute) __has_attribute(attribute)
#endif

#elif defined __GNUC__
#define HTS_GCC_AT_LEAST(major, minor) \
    (__GNUC__ > (major) || (__GNUC__ == (major) && __GNUC_MINOR__ >= (minor)))
#endif

#ifndef HTS_COMPILER_HAS
#define HTS_COMPILER_HAS(attribute) 0
#endif
#ifndef HTS_GCC_AT_LEAST
#define HTS_GCC_AT_LEAST(major, minor) 0
#endif

#if HTS_COMPILER_HAS(__noreturn__) || HTS_GCC_AT_LEAST(3,0)
#define HTS_NORETURN __attribute__ ((__noreturn__))
#else
#define HTS_NORETURN
#endif

// GCC introduced warn_unused_result in 3.4 but added -Wno-unused-result later
#if HTS_COMPILER_HAS(__warn_unused_result__) || HTS_GCC_AT_LEAST(4,5)
#define HTS_RESULT_USED __attribute__ ((__warn_unused_result__))
#else
#define HTS_RESULT_USED
#endif

#if HTS_COMPILER_HAS(__unused__) || HTS_GCC_AT_LEAST(3,0)
#define HTS_UNUSED __attribute__ ((__unused__))
#else
#define HTS_UNUSED
#endif

#if HTS_COMPILER_HAS(__deprecated__) || HTS_GCC_AT_LEAST(4,5)
#define HTS_DEPRECATED(message) __attribute__ ((__deprecated__ (message)))
#elif HTS_GCC_AT_LEAST(3,1)
#define HTS_DEPRECATED(message) __attribute__ ((__deprecated__))
#else
#define HTS_DEPRECATED(message)
#endif

#if HTS_COMPILER_HAS(__deprecated__) || HTS_GCC_AT_LEAST(6,4)
#define HTS_DEPRECATED_ENUM(message) __attribute__ ((__deprecated__ (message)))
#else
#define HTS_DEPRECATED_ENUM(message)
#endif

// On mingw the "printf" format type doesn't work.  It needs "gnu_printf"
// in order to check %lld and %z, otherwise it defaults to checking against
// the Microsoft library printf format options despite linking against the
// GNU posix implementation of printf.  The __MINGW_PRINTF_FORMAT macro
// expands to printf or gnu_printf as required, but obviously may not
// exist
#ifdef __MINGW_PRINTF_FORMAT
#define HTS_PRINTF_FMT __MINGW_PRINTF_FORMAT
#else
#define HTS_PRINTF_FMT printf
#endif

#if HTS_COMPILER_HAS(__format__) || HTS_GCC_AT_LEAST(3,0)
#define HTS_FORMAT(type, idx, first) __attribute__((__format__ (type, idx, first)))
#else
#define HTS_FORMAT(type, idx, first)
#endif

#endif


#ifdef __cplusplus
extern "C" {
#endif

/** @file

  Index FASTA or FASTQ files and extract subsequence.

  The fti file index columns for FASTA are:
    - chromosome name
    - chromosome length: number of bases
    - offset: number of bytes to skip to get to the first base
        from the beginning of the file, including the length
        of the sequence description string (`>chr ..\n`)
    - line length: number of bases per line (excluding `\n`)
    - binary line length: number of bytes, including `\n`

   The index for FASTQ is similar to above:
    - chromosome name
    - chromosome length: number of bases
    - sequence offset: number of bytes to skip to get to the first base
        from the beginning of the file, including the length
        of the sequence description string (`@chr ..\n`)
    - line length: number of bases per line (excluding `\n`)
    - binary line length: number of bytes, including `\n`
    - quality offset: number of bytes to skip from the beginning of the file
        to get to the first quality value in the indexed entry.

    The FASTQ version of the index uses line length and binary line length
    for both the sequence and the quality values, so they must be line
    wrapped in the same way.
 */

struct __ftidx_t;
/// Opaque structure representing FASTA index
typedef struct __ftidx_t ftidx_t;

/// File format to be dealing with.
enum fti_format_options {
    FTI_NONE,
    FTI_FASTA,
    FTI_FASTQ
};

/// Build index for a FASTA or FASTQ or bgzip-compressed FASTA or FASTQ file.
/**  @param  fn  FASTA/FASTQ file name
     @param  fnfti Name of .fti file to build.
     @param  fngzi Name of .gzi file to build (if fn is bgzip-compressed).
     @return     0 on success; or -1 on ftilure

If fnfti is NULL, ".fti" will be appended to fn to make the FTI file name.
If fngzi is NULL, ".gzi" will be appended to fn for the GZI file.  The GZI
file will only be built if fn is bgzip-compressed.
*/
int fti_build3(const char *fn, const char *fnfti, const char *fngzi) HTS_RESULT_USED;

/// Build index for a FASTA or FASTQ or bgzip-compressed FASTA or FASTQ file.
/** @param  fn  FASTA/FASTQ file name
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

/// Load FASTA indexes.
/** @param  fn  File name of the FASTA file (can be compressed with bgzip).
    @param  fnfti File name of the FASTA index.
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
/** @param  fn  File name of the FASTA file
    @return Pointer to a ftidx_t struct on success, NULL on ftilure.

This function is equivalent to fti_load3(fn, NULL, NULL, FTI_CREATE|FTI_CACHE);
*/
ftidx_t *fti_load(const char *fn);

/// Load FASTA or FASTQ indexes.
/** @param  fn  File name of the FASTA/FASTQ file (can be compressed with bgzip).
    @param  fnfti File name of the FASTA/FASTQ index.
    @param  fngzi File name of the bgzip index.
    @param  flags Option flags to control index file caching and creation.
    @param  format FASTA or FASTQ file format
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
/** @param  fn  File name of the FASTA/FASTQ file
    @param  format FASTA or FASTQ file format
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
