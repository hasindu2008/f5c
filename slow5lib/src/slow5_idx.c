#define _XOPEN_SOURCE 700
#include <unistd.h>
#include <inttypes.h>
//#include "klib/khash.h"
#include "slow5_idx.h"
//#include "slow5.h"
#include "slow5_extra.h"
#include "slow5_misc.h"
//#include "slow5_error.h"

extern enum slow5_log_level_opt  slow5_log_level;
extern enum slow5_exit_condition_opt  slow5_exit_condition;

#define BUF_INIT_CAP (20*1024*1024)
#define SLOW5_INDEX_BUF_INIT_CAP (64) // 2^6 TODO is this too little?

static inline struct slow5_idx *slow5_idx_init_empty(void);
static int slow5_idx_build(struct slow5_idx *index, struct slow5_file *s5p);
static int slow5_idx_read(struct slow5_idx *index);

static inline struct slow5_idx *slow5_idx_init_empty(void) {

    struct slow5_idx *index = (struct slow5_idx *) calloc(1, sizeof *index);
    SLOW5_MALLOC_CHK(index);
    index->hash = kh_init(slow5_s2i);

    return index;
}

// TODO return NULL if idx_init fails
struct slow5_idx *slow5_idx_init(struct slow5_file *s5p) {

    struct slow5_idx *index = slow5_idx_init_empty();
    if (!index) {
        return NULL;
    }
    index->pathname = slow5_get_idx_path(s5p->meta.pathname);
    if (!index->pathname) {
        slow5_idx_free(index);
        return NULL;
    }

    FILE *index_fp;

    // If file doesn't exist
    if ((index_fp = fopen(index->pathname, "r")) == NULL) {
        SLOW5_INFO("Index file not found. Creating an index at '%s'.", index->pathname)
        if (slow5_idx_build(index, s5p) != 0) {
            slow5_idx_free(index);
            return NULL;
        }
        index->fp = fopen(index->pathname, "w");
        if (slow5_idx_write(index, s5p->header->version) != 0) {
            slow5_idx_free(index);
            return NULL;
        }
        fclose(index->fp);
        index->fp = NULL;
    } else {
        index->fp = index_fp;
        int err;
        if (slow5_filestamps_cmp(index->pathname, s5p->meta.pathname, &err) < 0.0) {
            SLOW5_WARNING("Index file '%s' is older than slow5 file '%s'.",
                    index->pathname, s5p->meta.pathname);
        }
        if (err == -1) {
            slow5_idx_free(index);
            return NULL;
        }
        if (slow5_idx_read(index) != 0) {
            slow5_idx_free(index);
            return NULL;
        }
        if (index->version.major != s5p->header->version.major ||
                index->version.minor != s5p->header->version.minor ||
                index->version.patch != s5p->header->version.patch) {
            SLOW5_ERROR("Index file version '" SLOW5_VERSION_STRING_FORMAT "' is different to slow5 file version '" SLOW5_VERSION_STRING_FORMAT "'. Please re-index.",
                    index->version.major, index->version.minor, index->version.patch,
                    s5p->header->version.major, s5p->header->version.minor, s5p->header->version.patch);
            slow5_idx_free(index);
            return NULL;
        }
    }

    return index;
}

/**
 * Create the index file for slow5 file.
 * Overrides if already exists.
 *
 * @param   s5p         slow5 file structure
 * @param   pathname    pathname to write index to
 * @return  -1 on error, 0 on success
 */
int slow5_idx_to(struct slow5_file *s5p, const char *pathname) {

    struct slow5_idx *index = slow5_idx_init_empty();
    if (slow5_idx_build(index, s5p) != 0) {
        slow5_idx_free(index);
        return -1;
    }

    index->fp = fopen(pathname, "w");
    if (slow5_idx_write(index, s5p->header->version) != 0) {
        slow5_idx_free(index);
        return -1;
    }

    slow5_idx_free(index);
    return 0;
}

/*
 * return 0 on success
 * return <0 on failure
 * TODO fix error handling
 */
static int slow5_idx_build(struct slow5_idx *index, struct slow5_file *s5p) {

    uint64_t curr_offset = ftello(s5p->fp);
    if (fseeko(s5p->fp, s5p->meta.start_rec_offset, SEEK_SET != 0)) {
        return -1;
    }

    uint64_t offset = 0;
    uint64_t size = 0;

    if (s5p->format == SLOW5_FORMAT_ASCII) {
        size_t cap = BUF_INIT_CAP;
        char *buf = (char *) malloc(cap * sizeof *buf);
        SLOW5_MALLOC_CHK(buf);
        ssize_t buf_len;
        char *bufp;

        offset = ftello(s5p->fp);
        while ((buf_len = getline(&buf, &cap, s5p->fp)) != -1) { // TODO this return is closer int64_t not unsigned
            bufp = buf;
            char *read_id = strdup(slow5_strsep(&bufp, SLOW5_SEP_COL));
            size = buf_len;

            if (slow5_idx_insert(index, read_id, offset, size) == -1) {
                // TODO handle error and free
                SLOW5_ERROR("Inserting '%s' to index failed", read_id);
                return -1;
            }
            offset += buf_len;
        }

        free(buf);

    } else if (s5p->format == SLOW5_FORMAT_BINARY) {

        while (1) {
            // Set start offset
            offset = ftello(s5p->fp);

            // Get record size
            slow5_rec_size_t record_size;
            size_t bytes_read;
            if ((bytes_read = fread(&record_size, 1, sizeof record_size, s5p->fp)) != sizeof record_size) {
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
                        /* slow5_errno = SLOW5_ERR_EOF; */
                        break;
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
                return slow5_errno;
            }

            size = sizeof record_size + record_size;

            uint8_t *read_comp = (uint8_t *) malloc(record_size);
            SLOW5_MALLOC_CHK(read_comp);
            if (fread(read_comp, record_size, 1, s5p->fp) != 1) {
                free(read_comp);
                return -1;
            }

            uint8_t *read_decomp = read_comp;
            if (s5p->compress) {

                int8_t full_depress=!(s5p->compress->record_press->method == SLOW5_COMPRESS_ZLIB);
                if(!full_depress){
                    size_t part_len = 256; //guess max value to contain the read ID
                    part_len = part_len > record_size ? record_size : part_len;
                    size_t n=0;
                    read_decomp = (uint8_t *) slow5_ptr_depress(s5p->compress->record_press, read_comp, part_len, &n);
                    if(read_decomp == NULL || n < sizeof(slow5_rid_len_t) ){
                        SLOW5_WARNING("%s","Partial decompression failed. Falling back to the slower method. Open an issue on GitHub.");
                        free(read_decomp);
                        full_depress = 1;
                    }
                    else {
                        slow5_rid_len_t read_id_len;
                        memcpy(&read_id_len, read_decomp, sizeof read_id_len);
                        if(n < read_id_len + sizeof(slow5_rid_len_t)){
                            SLOW5_WARNING("%s","Partial decompression failed. Read ID too long. Falling back to the slower method. Open an issue on GitHub.");
                            free(read_decomp);
                            full_depress = 1;
                        }
                    }
                }

                if(full_depress){
                    read_decomp = (uint8_t *) slow5_ptr_depress(s5p->compress->record_press, read_comp, record_size, NULL);

                    if (read_decomp == NULL) {
                        free(read_comp);
                        free(read_decomp);
                        return -1;
                    }

                }
                free(read_comp);

            }

            // Get read id length
            uint64_t cur_len = 0;
            slow5_rid_len_t read_id_len;
            memcpy(&read_id_len, read_decomp + cur_len, sizeof read_id_len);
            cur_len += sizeof read_id_len;

            // Get read id
            char *read_id = (char *) malloc((read_id_len + 1) * sizeof *read_id); // +1 for '\0'
            SLOW5_MALLOC_CHK(read_id);
            memcpy(read_id, read_decomp + cur_len, read_id_len * sizeof *read_id);
            read_id[read_id_len] = '\0';

            // Insert index record
            if (slow5_idx_insert(index, read_id, offset, size) == -1) {
                SLOW5_ERROR("Inserting '%s' to index failed", read_id);
                // TODO handle error and free
                return -1;
            }

            free(read_decomp);
        }
    }

    if (fseeko(s5p->fp, curr_offset, SEEK_SET != 0)) {
        return -1;
    }

    return 0;
}

/*
 * write an index to its file
 * returns 0 on success, <0 on error
 */
int slow5_idx_write(struct slow5_idx *index, struct slow5_version version) {

    const char magic[] = SLOW5_INDEX_MAGIC_NUMBER;
    if (fwrite(magic, sizeof *magic, sizeof magic, index->fp) != sizeof magic) {
        return SLOW5_ERR_IO;
    }

    if (fwrite(&version.major, sizeof version.major, 1, index->fp) != 1 ||
            fwrite(&version.minor, sizeof version.minor, 1, index->fp) != 1 ||
            fwrite(&version.patch, sizeof version.patch, 1, index->fp) != 1) {
        return SLOW5_ERR_IO;
    }

    uint8_t padding = SLOW5_INDEX_HEADER_SIZE_OFFSET -
            sizeof magic * sizeof *magic -
            sizeof version.major -
            sizeof version.minor -
            sizeof version.patch;
    uint8_t *zeroes = (uint8_t *) calloc(padding, sizeof *zeroes);
    SLOW5_MALLOC_CHK(zeroes);
    if (fwrite(zeroes, sizeof *zeroes, padding, index->fp) != padding) {
        return SLOW5_ERR_IO;
    }
    free(zeroes);

    for (uint64_t i = 0; i < index->num_ids; ++ i) {

        khint_t pos = kh_get(slow5_s2i, index->hash, index->ids[i]);
        if (pos == kh_end(index->hash)) {
            return SLOW5_ERR_NOTFOUND;
        }

        struct slow5_rec_idx read_index = kh_value(index->hash, pos);

        slow5_rid_len_t read_id_len = strlen(index->ids[i]);
        if (fwrite(&read_id_len, sizeof read_id_len, 1, index->fp) != 1 ||
                fwrite(index->ids[i], sizeof *index->ids[i], read_id_len, index->fp) != read_id_len ||
                fwrite(&read_index.offset, sizeof read_index.offset, 1, index->fp) != 1 ||
                fwrite(&read_index.size, sizeof read_index.size, 1, index->fp) != 1) {
            return SLOW5_ERR_IO;
        }
    }

    const char eof[] = SLOW5_INDEX_EOF;
    if (fwrite(eof, sizeof *eof, sizeof eof, index->fp) != sizeof eof) {
        return SLOW5_ERR_IO;
    }

    return 0;
}

static int slow5_idx_read(struct slow5_idx *index) {

    struct slow5_version max_supported = SLOW5_VERSION_ARRAY;
    const char magic[] = SLOW5_INDEX_MAGIC_NUMBER;
    char buf_magic[sizeof magic]; // TODO is this a vla?
    if (fread(buf_magic, sizeof *magic, sizeof magic, index->fp) != sizeof magic) {
        return SLOW5_ERR_IO;
    }
    if (memcmp(magic, buf_magic, sizeof *magic * sizeof magic) != 0) {
        return SLOW5_ERR_MAGIC;
    }

    if (fread(&index->version.major, sizeof index->version.major, 1, index->fp) != 1 ||
        fread(&index->version.minor, sizeof index->version.minor, 1, index->fp) != 1 ||
        fread(&index->version.patch, sizeof index->version.patch, 1, index->fp) != 1) {
        return SLOW5_ERR_IO;
    }

    if (slow5_is_version_compatible(index->version, max_supported) == 0){
        SLOW5_ERROR("Index file version '" SLOW5_VERSION_STRING_FORMAT "' is higher than the max slow5 version '" SLOW5_VERSION_STRING "' supported by this slow5lib! Please re-index or use a newer version of slow5lib.",
                index->version.major, index->version.minor, index->version.patch);
        return SLOW5_ERR_VERSION;
    }

    if (fseek(index->fp, SLOW5_INDEX_HEADER_SIZE_OFFSET, SEEK_SET) == -1) {
        return SLOW5_ERR_IO;
    }

    while (1) {
        slow5_rid_len_t read_id_len;
        if (fread(&read_id_len, sizeof read_id_len, 1, index->fp) != 1) {
            SLOW5_ERROR("Malformed slow5 index. Failed to read the read ID length.%s", feof(index->fp) ? " Missing index end of file marker." : "");
            if (feof(index->fp)) {
                slow5_errno = SLOW5_ERR_TRUNC;
            } else {
                slow5_errno = SLOW5_ERR_IO;
            }
            return slow5_errno;
        }
        char *read_id = (char *) malloc((read_id_len + 1) * sizeof *read_id); // +1 for '\0'
        SLOW5_MALLOC_CHK(read_id);

        size_t bytes_read;
        if ((bytes_read = fread(read_id, sizeof *read_id, read_id_len, index->fp)) != read_id_len) {
            free(read_id);
            bytes_read += sizeof read_id_len;
            const char eof[] = SLOW5_INDEX_EOF;
            if (bytes_read == sizeof eof) {
                /* check if eof marker */
                int is_eof = slow5_is_eof(index->fp, eof, sizeof eof);
                if (is_eof == -1) { /* io/mem error */
                    SLOW5_ERROR("%s", "Internal error while checking for index eof marker.");
                } else if (is_eof == -2) {
                    SLOW5_ERROR("%s", "Malformed index. End of file marker found, but end of file not reached.");
                } else if (is_eof == 1) {
                    /* slow5_errno = SLOW5_ERR_EOF; */
                    break;
                }
            } else {
                slow5_errno = SLOW5_ERR_IO;
            }
            return slow5_errno;
        }
        read_id[read_id_len] = '\0'; // Add null byte

        uint64_t offset;
        uint64_t size;

        if (fread(&offset, sizeof offset, 1, index->fp) != 1 ||
                fread(&size, sizeof size, 1, index->fp) != 1) {
            return SLOW5_ERR_IO;
        }

        if (slow5_idx_insert(index, read_id, offset, size) == -1) {
            SLOW5_ERROR("Inserting '%s' to index failed", read_id);
            // TODO handle error and free
            return -1;
        }
    }

    return 0;
}

int slow5_idx_insert(struct slow5_idx *index, char *read_id, uint64_t offset, uint64_t size) {

    int absent;
    khint_t k = kh_put(slow5_s2i, index->hash, read_id, &absent);
    if (absent == -1 || absent == 0) {
        // TODO error if read_id duplicated?
        SLOW5_ERROR("Read ID '%s' is duplicated", read_id);
        return -1;
    }

    struct slow5_rec_idx *read_index = &kh_value(index->hash, k);

    if (index->num_ids == index->cap_ids) {
        // Realloc ids array
        index->cap_ids = index->cap_ids ? index->cap_ids << 1 : 16; // TODO possibly integer overflow

        char **tmp = (char **) realloc(index->ids, index->cap_ids * sizeof *tmp);
        SLOW5_MALLOC_CHK(tmp);

        index->ids = tmp;
    }

    index->ids[index->num_ids ++] = read_id;

    read_index->offset = offset;
    read_index->size = size;

    return 0;
}

/*
 * index, read_id cannot be NULL
 * returns -1 if read_id not in the index hash map, 0 otherwise
 */
int slow5_idx_get(struct slow5_idx *index, const char *read_id, struct slow5_rec_idx *read_index) {
    int ret = 0;

    khint_t pos = kh_get(slow5_s2i, index->hash, read_id);
    if (pos == kh_end(index->hash)) {
        SLOW5_ERROR("Read ID '%s' was not found.", read_id)
        ret = -1;
    } else if (read_index) {
        *read_index = kh_value(index->hash, pos);
    }

    return ret;
}

/*
 * SLOW5_ERR_IO - issue closing index file pointer, check errno for details
 */
void slow5_idx_free(struct slow5_idx *index) {
    if (index == NULL) {
        return;
    }

    if (index->fp != NULL) {
        if (fclose(index->fp) == EOF) {
            SLOW5_ERROR("Failure when closing index file: %s", strerror(errno));
            slow5_errno = SLOW5_ERR_IO;
        }
    }

    for (uint64_t i = 0; i < index->num_ids; ++ i) {
        free(index->ids[i]);
    }
    free(index->ids);

    kh_destroy(slow5_s2i, index->hash);

    free(index->pathname);
    free(index);
}
