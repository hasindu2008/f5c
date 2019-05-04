/*  ftidx.c -- FASTT (FAST5 in TSV) random access.
	adpapted from htslib/faidx.c by Hasindu Gamaarachchi <hasindu@unsw.edu.au>
*/

/*  faidx.c -- FASTA and FASTQ random access.

    Copyright (C) 2008, 2009, 2013-2018 Genome Research Ltd.
    Portions copyright (C) 2011 Broad Institute.

    Author: Heng Li <lh3@sanger.ac.uk>

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

//#include <config.h>


//#define BGFS_HFILE 1
#define UN_BUFFERED 1

#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <errno.h>
#include <limits.h>
#include <unistd.h>
#include <assert.h>
#include <getopt.h>
#ifdef BGFS_HFILE
    #include "htslib/bgzf.h"
#endif
#include "ftidx.h"
#include "htslib/hfile.h"
#include "htslib/khash.h"
#include "htslib/kstring.h"
//#include "hts_internal.h"
#ifdef UN_BUFFERED
#include <unistd.h>
#include <fcntl.h>
#endif

#define hts_log_warning(arg, ...)                                                      \
    fprintf(stderr, "[%s::WARNING]\033[1;33m " arg "\033[0m\n", __func__,      \
            __VA_ARGS__)
#define hts_log_error(arg, ...)                                                        \
    fprintf(stderr, "[%s::ERROR]\033[1;31m " arg "\033[0m\n", __func__,        \
            __VA_ARGS__)
#define hts_log_info(arg, ...)                                                         \
    fprintf(stderr, "[%s::INFO]\033[1;34m " arg "\033[0m\n", __func__,         \
            __VA_ARGS__)
static inline int isspace_c(char c) { return isspace((unsigned char) c); }
static inline int isdigit_c(char c) { return isdigit((unsigned char) c); }


#ifndef BGFS_HFILE

typedef struct {

    int fd;

    FILE *fp;

    int is_compressed;
    int is_gzip;
} BGZF;    

    /**
     * Read one line from a BGZF file. It is faster than bgzf_getc()
     *
     * @param fp     BGZF file handler
     * @param delim  delimitor
     * @param str    string to write to; must be initialized
     * @return       length of the string; -1 on end-of-file; <= -2 on error
     */
    static inline int bgzf_getline(BGZF *fp, int delim, kstring_t *str){
        
        str->m = 20*1024*1024;
        str->s = (char *)malloc(sizeof(char)*20*1024*1024);
        str->l=getline(&(str->s),&(str->m),fp->fp);
        if(str->l==0 || str->m==0 || str->s==NULL ){
            hts_log_error("%s\n", "reading issue");
            exit(1);
        }
    }

    static inline size_t f5read(BGZF *fp, kstring_t *str, size_t num_elements){
        str->m = num_elements;
        str->s = (char *)malloc(sizeof(char)*str->m);
    #ifdef UN_BUFFERED
        size_t ret=read(fp->fd,str->s,num_elements);
    #else    
        size_t ret=fread(str->s,1,num_elements,fp->fp);
    #endif
        str->l = ret;
        if(ret!=num_elements){
            fprintf(stderr,"Reading error has occurred :%s\n",strerror(errno));
            exit(EXIT_FAILURE);
        }
        return ret;
    } 

    /**
     *  Position in uncompressed BGZF
     *
     *  @param fp           BGZF file handler; must be opened for reading
     *
     *  Returns the current offset on success and -1 on error.
     */
    long bgzf_utell(BGZF *fp){
        hts_log_error("%s\n", "Not implemented");
        exit(1);
    }

 /**
     * Close the BGZF and free all associated resources.
     *
     * @param fp    BGZF file handler
     * @return      0 on success and -1 on error
     */
    int bgzf_close(BGZF *fp){
    #ifdef UN_BUFFERED
        close(fp->fd);
    #else   
        fclose(fp->fp);
    #endif
        free(fp);
        //hts_log_error("%s\n", "Not implemented");
        //exit(1);
    }

    /**
     * Open the specified file for reading or writing.
     */
    BGZF* bgzf_open(const char* path, const char *mode){
        BGZF *fp = (BGZF *)malloc(sizeof(BGZF));
        fp->is_compressed=0;
        fp->is_gzip=0;
    #ifdef UN_BUFFERED   
        fp->fd = open(path,O_RDONLY );
        if(fp->fd<0){
            hts_log_error("File %s cannot be opened\n", path);
            exit(1);
        }
    #else 
        fp->fp = fopen(path,mode);
        if(fp->fp==NULL){
            hts_log_error("File %s cannot be opened\n", path);
            exit(1);
        }
    #endif

        return fp;
    }



    /** Return the file's compression format
     *
     * @param fp  BGZF file handle
     * @return    A small integer matching the corresponding
     *            `enum htsCompression` value:
     *   - 0 / `no_compression` if the file is uncompressed
     *   - 1 / `gzip` if the file is plain GZIP-compressed
     *   - 2 / `bgzf` if the file is BGZF-compressed
     * @since 1.4
     */
    int bgzf_compression(BGZF *fp){
        return 0;
    }

      /**
     *  Position BGZF at the uncompressed offset
     *
     *  @param fp           BGZF file handler; must be opened for reading
     *  @param uoffset      file offset in the uncompressed data
     *  @param where        SEEK_SET supported atm
     *
     *  Returns 0 on success and -1 on error.
     */
    int bgzf_useek(BGZF *fp, long uoffset, int where)  {
        #ifdef UN_BUFFERED
            long ret=lseek(fp->fd, uoffset, SEEK_SET);
            if(ret>=0){
                return 0;
            }
            else{
                return -1;
            }
        #else
            return fseek(fp->fp, uoffset, SEEK_SET);
        #endif
        //hts_log_error("%s\n", "Not implemented");
        //exit(1);  
    }


    /**
     * Tell BGZF to build index while compressing.
     *
     * @param fp          BGZF file handler; can be opened for reading or writing.
     *
     * Returns 0 on success and -1 on error.
     */
    int bgzf_index_build_init(BGZF *fp){
        hts_log_error("%s\n", "Not implemented");
        exit(1);
    }

    /// Save BGZF index
    /**
     * @param fp          BGZF file handler
     * @param bname       base name
     * @param suffix      suffix to add to bname (can be NULL)
     * @return 0 on success and -1 on error.
     */
    int bgzf_index_dump(BGZF *fp,
                        const char *bname, const char *suffix) {
        hts_log_error("%s\n", "Not implemented");
        exit(1);
                        }

/// Load BGZF index
    /**
     * @param fp          BGZF file handler
     * @param bname       base name
     * @param suffix      suffix to add to bname (can be NULL)
     * @return 0 on success and -1 on error.
     */
    int bgzf_index_load(BGZF *fp,
                        const char *bname, const char *suffix){
        hts_log_error("%s\n", "Not implemented");
        exit(1);
        }

#endif



typedef struct {
    uint32_t line_len, line_blen;
    uint64_t len;
    uint64_t seq_offset;
    uint64_t qual_offset;
} ftidx1_t;
KHASH_MAP_INIT_STR(s, ftidx1_t)

struct __ftidx_t {
    BGZF *bgzf;
    int n, m;
    char **name;
    khash_t(s) *hash;
    enum fti_format_options format;
};

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

static inline int fti_insert_index(ftidx_t *idx, const char *name, uint64_t len, uint32_t line_len, uint32_t line_blen, uint64_t seq_offset, uint64_t qual_offset)
{
    if (!name) {
        hts_log_error("%s","Malformed line");
        return -1;
    }

    char *name_key = strdup(name);
    int absent;
    khint_t k = kh_put(s, idx->hash, name_key, &absent);
    ftidx1_t *v = &kh_value(idx->hash, k);

    if (! absent) {
        hts_log_warning("Ignoring duplicate sequence \"%s\" at byte offset %" PRIu64 "", name, seq_offset);
        free(name_key);
        return 0;
    }

    if (idx->n == idx->m) {
        char **tmp;
        idx->m = idx->m? idx->m<<1 : 16;
        if (!(tmp = (char**)realloc(idx->name, sizeof(char*) * idx->m))) {
            hts_log_error("%s","Out of memory");
            return -1;
        }
        idx->name = tmp;
    }
    idx->name[idx->n++] = name_key;
    v->len = len;
    v->line_len = line_len;
    v->line_blen = line_blen;
    v->seq_offset = seq_offset;
    v->qual_offset = qual_offset;

    return 0;
}


static ftidx_t *fti_build_core(BGZF *bgzf) {
    kstring_t name = { 0, 0, NULL };
    kstring_t linebuffer = { 0, 0, NULL };
    //int c, read_done, line_num;
    ftidx_t *idx;
    uint64_t seq_offset, qual_offset;
    uint64_t seq_len, qual_len;
    uint64_t char_len, cl, line_len, ll;
    //enum read_state {OUT_READ, IN_NAME, IN_SEQ, SEQ_END, IN_QUAL} state;

    idx = (ftidx_t*)calloc(1, sizeof(ftidx_t));
    idx->hash = kh_init(s);
    idx->format = FTI_NONE;

    //state = OUT_READ, read_done = 0, line_num = 1;
    seq_offset = qual_offset = seq_len = qual_len = char_len = cl = line_len = ll = 0;

    linebuffer.l=0;

    while (bgzf_getline(bgzf, '\n', &linebuffer) >0 ) {
        if (linebuffer.s[0] == '#' || linebuffer.s[0] == '\n' || linebuffer.s[0] == '\r') { //comments and header
            //fprintf(stderr,"%s\n",linebuffer.s);
            //line_num++;
            linebuffer.l=0;
            seq_offset = bgzf_utell(bgzf);
            continue;
        }
        else{
                
                char *name=strtok(linebuffer.s,"\t");
                line_len=linebuffer.l;
                //fprintf(stderr,"%s %ld\n",name,seq_offset);
                if (fti_insert_index(idx, name, seq_len, line_len, char_len, seq_offset, qual_offset) != 0){
                    goto ftil;
                }
                seq_len = qual_len = char_len = line_len = 0;
                seq_offset = bgzf_utell(bgzf);
                //line_num++;
                linebuffer.l=0;
                //name.l=0;
        }
    }

    // while ((c = bgzf_getc(bgzf)) >= 0) {
    //     switch (state) {
    //         case OUT_READ:
    //             switch (c) {

    //                 case '#':
    //                     idx->format = FTI_FASTT;
    //                     state = IN_NAME;
    //                 break;

    //                 case '\r':
    //                     // Blank line with cr-lf ending?
    //                     if ((c = bgzf_getc(bgzf)) == '\n') {
    //                         line_num++;
    //                     } else {
    //                         hts_log_error("Format error, carriage return not followed by new line at line %d", line_num);
    //                         goto ftil;
    //                     }
    //                 break;

    //                 case '\n':
    //                     // just move onto the next line
    //                     line_num++;
    //                 break;

    //                 default: {
    //                     char s[4] = { '"', c, '"', '\0' };
    //                     hts_log_error("Format error, unexpected %s at line %d", isprint(c) ? s : "character", line_num);
    //                     goto ftil;
    //                 break;
    //                 }
    //             }
    //         break;

    //         case IN_NAME:
    //             if (read_done) {
    //                 if (fti_insert_index(idx, name.s, seq_len, line_len, char_len, seq_offset, qual_offset) != 0)
    //                     goto ftil;

    //                 read_done = 0;
    //             }

    //             name.l = 0;

    //             do {
    //                 if (!isspace(c)) {
    //                     kputc(c, &name);
    //                 } else if (name.l > 0 || c == '\n') {
    //                     break;
    //                 }
    //             } while ((c = bgzf_getc(bgzf)) >= 0);

    //             kputsn("", 0, &name);

    //             if (c < 0) {
    //                 hts_log_error("The last entry '%s' has no sequence", name.s);
    //                 goto ftil;
    //             }

    //             // read the rest of the line if necessary
    //             if (c != '\n') while ((c = bgzf_getc(bgzf)) >= 0 && c != '\n');

    //             state = OUT_READ; seq_len = qual_len = char_len = line_len = 0;
    //             seq_offset = bgzf_utell(bgzf);
    //             line_num++;
    //         break;

    //     }
    // }

    // if (read_done) {
    //     if (fti_insert_index(idx, name.s, seq_len, line_len, char_len, seq_offset, qual_offset) != 0)
    //         goto ftil;
    // } else {
    //     goto ftil;
    // }

    free(name.s);
    free(linebuffer.s);
    return idx;

ftil:
    free(name.s);
    free(linebuffer.s);
    fti_destroy(idx);
    return NULL;
}


static int fti_save(const ftidx_t *fti, hFILE *fp) {
    khint_t k;
    int i;
    char buf[96]; // Must be big enough for format below.

    for (i = 0; i < fti->n; ++i) {
        ftidx1_t x;
        k = kh_get(s, fti->hash, fti->name[i]);
        assert(k < kh_end(fti->hash));
        x = kh_value(fti->hash, k);

        if (fti->format == FTI_FASTT) {
            snprintf(buf, sizeof(buf),
                 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu32 "\t%" PRIu32 "\n",
                 x.len, x.seq_offset, x.line_blen, x.line_len);
        } else {
            snprintf(buf, sizeof(buf),
                 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu32 "\t%" PRIu32 "\t%" PRIu64 "\n",
                 x.len, x.seq_offset, x.line_blen, x.line_len, x.qual_offset);
        }

        if (hputs(fti->name[i], fp) != 0) return -1;
        if (hputs(buf, fp) != 0) return -1;
    }
    return 0;
}


static ftidx_t *fti_read(hFILE *fp, const char *fname, int format)
{
    ftidx_t *fti;
    char *buf = NULL, *p;
    ssize_t l, lnum = 1;

    fti = (ftidx_t*)calloc(1, sizeof(ftidx_t));
    if (!fti) return NULL;

    fti->hash = kh_init(s);
    if (!fti->hash) goto ftil;

    buf = (char*)calloc(0x10000, 1);
    if (!buf) goto ftil;

    while ((l = hgetln(buf, 0x10000, fp)) > 0) {
        uint32_t line_len, line_blen, n;
        uint64_t len;
        uint64_t seq_offset;
        uint64_t qual_offset = 0;

        for (p = buf; *p && !isspace_c(*p); ++p);

        if (p - buf < l) {
            *p = 0; ++p;
        }

        if (format == FTI_FASTT) {
            n = sscanf(p, "%" SCNu64 "%" SCNu64 "%" SCNu32 "%" SCNu32 "", &len, &seq_offset, &line_blen, &line_len);

            if (n != 4) {
                hts_log_error("Could not understand FASTT index %s line %zd", fname, lnum);
                goto ftil;
            }
        } else {
            n = sscanf(p, "%" SCNu64 "%" SCNu64 "%" SCNu32 "%" SCNu32 "%" SCNu64 "", &len, &seq_offset, &line_blen, &line_len, &qual_offset);

            if (n != 5) {
                if (n == 4) {
                    hts_log_error("Possibly this is a FASTT index, try using ftidx.  Problem in %s line %zd", fname, lnum);
                } else {
                    hts_log_error("Could not understand FASTB index %s line %zd", fname, lnum);
                }

                goto ftil;
            }
        }

        if (fti_insert_index(fti, buf, len, line_len, line_blen, seq_offset, qual_offset) != 0) {
            goto ftil;
        }

        if (buf[l - 1] == '\n') ++lnum;
    }

    if (l < 0) {
        hts_log_error("Error while reading %s: %s", fname, strerror(errno));
        goto ftil;
    }
    free(buf);
    return fti;

 ftil:
    free(buf);
    fti_destroy(fti);
    return NULL;
}

void fti_destroy(ftidx_t *fti)
{
    int i;
    if (!fti) return;
    for (i = 0; i < fti->n; ++i) free(fti->name[i]);
    free(fti->name);
    kh_destroy(s, fti->hash);
    if (fti->bgzf) bgzf_close(fti->bgzf);
    free(fti);
}


static int fti_build3_core(const char *fn, const char *fnfti, const char *fngzi)
{
    kstring_t fti_kstr = { 0, 0, NULL };
    kstring_t gzi_kstr = { 0, 0, NULL };
    BGZF *bgzf = NULL;
    hFILE *fp = NULL;
    ftidx_t *fti = NULL;
    int save_errno, res;
    const char *file_type;

    bgzf = bgzf_open(fn, "r");

    if ( !bgzf ) {
        hts_log_error("Failed to open the file %s", fn);
        goto ftil;
    }

    if ( bgzf->is_compressed ) {
        if (bgzf_index_build_init(bgzf) != 0) {
            hts_log_error("%s","Failed to allocate bgzf index");
            goto ftil;
        }
    }

    fti = fti_build_core(bgzf);

    if ( !fti ) {
        if (bgzf->is_compressed && bgzf->is_gzip) {
            hts_log_error("%s","Cannot index files compressed with gzip, please use bgzip");
        }
        goto ftil;
    }

    if (fti->format == FTI_FASTT) {
        file_type   = "FASTT";
    } else {
        file_type   = "FASTB";
    }

    if (!fnfti) {
        if (ksprintf(&fti_kstr, "%s.fti", fn) < 0) goto ftil;
        fnfti = fti_kstr.s;
    }

    if (!fngzi) {
        if (ksprintf(&gzi_kstr, "%s.gzi", fn) < 0) goto ftil;
        fngzi = gzi_kstr.s;
    }

    if ( bgzf->is_compressed ) {
        if (bgzf_index_dump(bgzf, fngzi, NULL) < 0) {
            hts_log_error("Failed to make bgzf index %s", fngzi);
            goto ftil;
        }
    }

    res = bgzf_close(bgzf);
    bgzf = NULL;

    if (res < 0) {
        hts_log_error("Error on closing %s : %s", fn, strerror(errno));
        goto ftil;
    }

    fp = hopen(fnfti, "wb");

    if ( !fp ) {
        hts_log_error("Failed to open %s index %s : %s", file_type, fnfti, strerror(errno));
        goto ftil;
    }

    if (fti_save(fti, fp) != 0) {
        hts_log_error("Failed to write %s index %s : %s", file_type, fnfti, strerror(errno));
        goto ftil;
    }

    if (hclose(fp) != 0) {
        hts_log_error("Failed on closing %s index %s : %s", file_type, fnfti, strerror(errno));
        goto ftil;
    }

    free(fti_kstr.s);
    free(gzi_kstr.s);
    fti_destroy(fti);
    return 0;

 ftil:
    save_errno = errno;
    free(fti_kstr.s);
    free(gzi_kstr.s);
    bgzf_close(bgzf);
    fti_destroy(fti);
    errno = save_errno;
    return -1;
}


int fti_build3(const char *fn, const char *fnfti, const char *fngzi) {
    return fti_build3_core(fn, fnfti, fngzi);
}


int fti_build(const char *fn) {
    return fti_build3(fn, NULL, NULL);
}


static ftidx_t *fti_load3_core(const char *fn, const char *fnfti, const char *fngzi,
                   int flags, int format)
{
    kstring_t fti_kstr = { 0, 0, NULL };
    kstring_t gzi_kstr = { 0, 0, NULL };
    hFILE *fp = NULL;
    ftidx_t *fti = NULL;
    int res, gzi_index_needed = 0;
    const char *file_type;

    if (format == FTI_FASTT) {
        file_type   = "FASTT";
    } else {
        file_type   = "FASTB";
    }

    if (fn == NULL)
        return NULL;

    if (fnfti == NULL) {
        if (ksprintf(&fti_kstr, "%s.fti", fn) < 0) goto ftil;
        fnfti = fti_kstr.s;
    }
    if (fngzi == NULL) {
        if (ksprintf(&gzi_kstr, "%s.gzi", fn) < 0) goto ftil;
        fngzi = gzi_kstr.s;
    }

    fp = hopen(fnfti, "rb");

    if (fp) {
        // index file present, check if a compressed index is needed
        hFILE *gz = NULL;
        BGZF *bgzf = bgzf_open(fn, "rb");

        if (bgzf == 0) {
            hts_log_error("Failed to open %s file %s", file_type, fn);
            goto ftil;
        }

        if (bgzf_compression(bgzf) == 2) { // BGZF compression
            if ((gz = hopen(fngzi, "rb")) == 0) {

                if (!(flags & FTI_CREATE) || errno != ENOENT) {
                    hts_log_error("Failed to open %s index %s: %s", file_type, fngzi, strerror(errno));
                    bgzf_close(bgzf);
                    goto ftil;
                }

                gzi_index_needed = 1;
                res = hclose(fp); // closed as going to be re-indexed

                if (res < 0) {
                    hts_log_error("Failed on closing %s index %s : %s", file_type, fnfti, strerror(errno));
                    goto ftil;
                }
            } else {
                res = hclose(gz);

                if (res < 0) {
                    hts_log_error("Failed on closing %s index %s : %s", file_type, fngzi, strerror(errno));
                    goto ftil;
                }
            }
        }

        bgzf_close(bgzf);
    }

    if (fp == 0 || gzi_index_needed) {
        if (!(flags & FTI_CREATE) || errno != ENOENT) {
            hts_log_error("Failed to open %s index %s: %s", file_type, fnfti, strerror(errno));
            goto ftil;
        }

        hts_log_info("Build %s index", file_type);

        if (fti_build3_core(fn, fnfti, fngzi) < 0) {
            goto ftil;
        }

        fp = hopen(fnfti, "rb");
        if (fp == 0) {
            hts_log_error("Failed to open %s index %s: %s", file_type, fnfti, strerror(errno));
            goto ftil;
        }
    }

    fti = fti_read(fp, fnfti, format);
    if (fti == NULL) {
        hts_log_error("Failed to read %s index %s", file_type, fnfti);
        goto ftil;
    }

    res = hclose(fp);
    fp = NULL;
    if (res < 0) {
        hts_log_error("Failed on closing %s index %s : %s", file_type, fnfti, strerror(errno));
        goto ftil;
    }

    fti->bgzf = bgzf_open(fn, "rb");
    if (fti->bgzf == 0) {
        hts_log_error("Failed to open %s file %s", file_type, fn);
        goto ftil;
    }

    if ( fti->bgzf->is_compressed==1 ) {
        if ( bgzf_index_load(fti->bgzf, fngzi, NULL) < 0 ) {
            hts_log_error("Failed to load .gzi index: %s", fngzi);
            goto ftil;
        }
    }
    free(fti_kstr.s);
    free(gzi_kstr.s);
    return fti;

 ftil:
    if (fti) fti_destroy(fti);
    if (fp) hclose_abruptly(fp);
    free(fti_kstr.s);
    free(gzi_kstr.s);
    return NULL;
}


ftidx_t *fti_load3(const char *fn, const char *fnfti, const char *fngzi,
                   int flags) {
    return fti_load3_core(fn, fnfti, fngzi, flags, FTI_FASTT);
}


ftidx_t *fti_load(const char *fn)
{
    return fti_load3(fn, NULL, NULL, FTI_CREATE);
}


ftidx_t *fti_load3_format(const char *fn, const char *fnfti, const char *fngzi,
                   int flags, enum fti_format_options format) {
    return fti_load3_core(fn, fnfti, fngzi, flags, format);
}


ftidx_t *fti_load_format(const char *fn, enum fti_format_options format) {
    return fti_load3_format(fn, NULL, NULL, FTI_CREATE, format);
}


static char *fti_retrieve(const ftidx_t *fti, const ftidx1_t *val,
                          uint64_t offset, long beg, long end, int *len) {
    char *s;
    size_t l;
    int c = 0;
    // int ret = bgzf_useek(fti->bgzf,
    //                      offset
    //                      + beg / val->line_blen * val->line_len
    //                      + beg % val->line_blen, SEEK_SET);
    int ret = bgzf_useek(fti->bgzf,
                         offset, SEEK_SET);
    if (ret < 0) {
        *len = -1;
        hts_log_error("%s","Failed to retrieve block. (Seeking in a compressed, .gzi unindexed, file?)");
        return NULL;
    }

    // l = 0;
    // s = (char*)malloc((size_t) end - beg + 2);
    // if (!s) {
    //     *len = -1;
    //     return NULL;
    // }

    kstring_t linebuffer = { 0, 0, NULL };
#ifdef BGFS_HFILE
    ret=bgzf_getline(fti->bgzf, '\n', &linebuffer);
#else
    ret=f5read(fti->bgzf, &linebuffer, val->line_len);
#endif

    // while ( l < end - beg && (c=bgzf_getc(fti->bgzf))>=0 )
    //     if (isgraph(c)) s[l++] = c;
    // if (c < 0) {
    if(ret<0){
        hts_log_error("Failed to retrieve block: %s",
            c == -1 ? "unexpected end of file" : "error reading file");
        if(linebuffer.s){
            free(linebuffer.s);
        }
        *len = -1;
        return NULL;
    }

    l=linebuffer.l;
    s=linebuffer.s;
    //s[l] = '\0';
    *len = l < INT_MAX ? l : INT_MAX;
    return s;
}


static int fti_get_val(const ftidx_t *fti, const char *str, int *len, ftidx1_t *val, long *fbeg, long *fend) {
    char *s, *ep;
    size_t i, l, k, name_end;
    khiter_t iter;
    khash_t(s) *h;
    long beg, end;

    beg = end = -1;
    h = fti->hash;
    name_end = l = strlen(str);
    s = (char*)malloc(l+1);
    if (!s) {
        *len = -1;
        return 1;
    }

    // remove space
    for (i = k = 0; i < l; ++i)
        if (!isspace_c(str[i])) s[k++] = str[i];
    s[k] = 0;
    name_end = l = k;
    // determine the sequence name
    for (i = l; i > 0; --i) if (s[i - 1] == ':') break; // look for colon from the end
    if (i > 0) name_end = i - 1;
    if (name_end < l) { // check if this is really the end
        int n_hyphen = 0;
        for (i = name_end + 1; i < l; ++i) {
            if (s[i] == '-') ++n_hyphen;
            else if (!isdigit_c(s[i]) && s[i] != ',') break;
        }
        if (i < l || n_hyphen > 1) name_end = l; // malformated region string; then take str as the name
        s[name_end] = 0;
        iter = kh_get(s, h, s);
        if (iter == kh_end(h)) { // cannot find the sequence name
            iter = kh_get(s, h, str); // try str as the name
            if (iter != kh_end(h)) {
                s[name_end] = ':';
                name_end = l;
            }
        }
    } else iter = kh_get(s, h, str);
    if(iter == kh_end(h)) {
        hts_log_warning("Reference %s not found in file, returning empty sequence", str);
        free(s);
        *len = -2;
        return 1;
    }
    *val = kh_value(h, iter);
    // parse the interval
    if (name_end < l) {
        int save_errno = errno;
        errno = 0;
        for (i = k = name_end + 1; i < l; ++i)
            if (s[i] != ',') s[k++] = s[i];
        s[k] = 0;
        if (s[name_end + 1] == '-') {
            beg = 0;
            i = name_end + 2;
        } else {
            beg = strtol(s + name_end + 1, &ep, 10);
            for (i = ep - s; i < k;) if (s[i++] == '-') break;
        }
        end = i < k? strtol(s + i, &ep, 10) : val->len;
        if (beg > 0) --beg;
        // Check for out of range numbers.  Only going to be a problem on
        // 32-bit platforms with >2Gb sequence length.
        if (errno == ERANGE && (uint64_t) val->len > LONG_MAX) {
            hts_log_error("Positions in range %s are too large for this platform", s);
            free(s);
            *len = -3;
            return 1;
        }
        errno = save_errno;
    } else beg = 0, end = val->len;
    if (beg >= val->len) beg = val->len;
    if (end >= val->len) end = val->len;
    if (beg > end) beg = end;
    free(s);

    *fbeg = beg;
    *fend = end;

    return 0;
}


char *fti_fetch(const ftidx_t *fti, const char *str, int *len)
{
    ftidx1_t val;
    long beg, end;

    if (fti_get_val(fti, str, len, &val, &beg, &end)) {
        return NULL;
    }

    // now retrieve the sequence
    return fti_retrieve(fti, &val, val.seq_offset, beg, end, len);
}


char *fti_fetchqual(const ftidx_t *fti, const char *str, int *len) {
    ftidx1_t val;
    long beg, end;

    if (fti_get_val(fti, str, len, &val, &beg, &end)) {
        return NULL;
    }

    // now retrieve the sequence
    return fti_retrieve(fti, &val, val.qual_offset, beg, end, len);
}


int ftidx_fetch_nseq(const ftidx_t *fti)
{
    return fti->n;
}

int ftidx_nseq(const ftidx_t *fti)
{
    return fti->n;
}

const char *ftidx_iseq(const ftidx_t *fti, int i)
{
    return fti->name[i];
}

int ftidx_seq_len(const ftidx_t *fti, const char *seq)
{
    khint_t k = kh_get(s, fti->hash, seq);
    if ( k == kh_end(fti->hash) ) return -1;
    return kh_val(fti->hash, k).len;
}


static int ftidx_adjust_position(const ftidx_t *fti, ftidx1_t *val, const char *c_name, int *p_beg_i, int *p_end_i, int *len) {
    khiter_t iter;

    // Adjust position
    iter = kh_get(s, fti->hash, c_name);

    if (iter == kh_end(fti->hash)) {
        *len = -2;
        hts_log_error("The sequence \"%s\" not found", c_name);
        return 1;
    }

    *val = kh_value(fti->hash, iter);

    if(*p_end_i < *p_beg_i)
        *p_beg_i = *p_end_i;

    if(*p_beg_i < 0)
        *p_beg_i = 0;
    else if(val->len <= *p_beg_i)
        *p_beg_i = val->len - 1;

    if(*p_end_i < 0)
        *p_end_i = 0;
    else if(val->len <= *p_end_i)
        *p_end_i = val->len - 1;

    return 0;
}


char *ftidx_fetch_seq(const ftidx_t *fti, const char *c_name, int p_beg_i, int p_end_i, int *len)
{
    ftidx1_t val;

    // Adjust position
    if (ftidx_adjust_position(fti, &val, c_name, &p_beg_i, &p_end_i, len)) {
        return NULL;
    }

    // Now retrieve the sequence
    return fti_retrieve(fti, &val, val.seq_offset, p_beg_i, (long) p_end_i + 1, len);
}


char *ftidx_fetch_qual(const ftidx_t *fti, const char *c_name, int p_beg_i, int p_end_i, int *len)
{
    ftidx1_t val;

    // Adjust position
    if (ftidx_adjust_position(fti, &val, c_name, &p_beg_i, &p_end_i, len)) {
        return NULL;
    }

    // Now retrieve the sequence
    return fti_retrieve(fti, &val, val.qual_offset, p_beg_i, (long) p_end_i + 1, len);
}


int ftidx_has_seq(const ftidx_t *fti, const char *seq)
{
    khiter_t iter = kh_get(s, fti->hash, seq);
    if (iter == kh_end(fti->hash)) return 0;
    return 1;
}

