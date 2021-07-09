#define _XOPEN_SOURCE 700
#include <unistd.h>
#include <inttypes.h>
//#include "klib/khash.h"
#include "slow5_idx.h"
//#include "slow5.h"
#include "slow5_extra.h"
#include "slow5_misc.h"
//#include "slow5_error.h"
//TODO MALLOC_CHK for testing

extern enum slow5_log_level_opt  slow5_log_level;
extern enum slow5_exit_condition_opt  slow5_exit_condition;

#define BUF_INIT_CAP (20*1024*1024)
#define SLOW5_INDEX_BUF_INIT_CAP (64) // 2^6 TODO is this too little?

static inline struct slow5_idx *slow5_idx_init_empty(void);
static int slow5_idx_build(struct slow5_idx *index, struct slow5_file *s5p);
static void slow5_idx_read(struct slow5_idx *index);

static inline struct slow5_idx *slow5_idx_init_empty(void) {

    struct slow5_idx *index = (struct slow5_idx *) calloc(1, sizeof *index);
    SLOW5_MALLOC_CHK(index);
    index->hash = kh_init(slow5_s2i);

    return index;
}

// TODO return NULL if idx_init fails
struct slow5_idx *slow5_idx_init(struct slow5_file *s5p) {

    struct slow5_idx *index = slow5_idx_init_empty();
    index->pathname = slow5_get_idx_path(s5p->meta.pathname);

    if(index==NULL || index->pathname==NULL ){
        //TODO fix mem leak
        return NULL;
    }

    FILE *index_fp;

    // If file doesn't exist
    if ((index_fp = fopen(index->pathname, "r")) == NULL) {
        SLOW5_INFO("Index file not found. Creating an index at %s.",index->pathname)
        if (slow5_idx_build(index, s5p) != 0) {
            slow5_idx_free(index);
            return NULL;
        }
        index->fp = fopen(index->pathname, "w");
        slow5_idx_write(index);
        fclose(index->fp);
        index->fp = NULL;
    } else {
        index->fp = index_fp;
        slow5_idx_read(index);
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
    if (slow5_idx_build(index, s5p) == -1) {
        slow5_idx_free(index);
        return -1;
    }

    index->fp = fopen(pathname, "w");
    slow5_idx_write(index);

    slow5_idx_free(index);
    return 0;
}

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
            char *read_id = strdup(slow5_strsep(&bufp, SLOW5_SEP_COL)); // TODO quicker to not split the whole line just the first delim
            size = buf_len;

            slow5_idx_insert(index, read_id, offset, size);
            offset += buf_len;
        }

        free(buf);

    } else if (s5p->format == SLOW5_FORMAT_BINARY) {
        const char eof[] = SLOW5_BINARY_EOF;
        char buf_eof[sizeof eof]; // TODO is this a vla?

        if (fread(buf_eof, sizeof *eof, sizeof eof, s5p->fp) != sizeof eof) {
            return -1;
        }
        while (memcmp(eof, buf_eof, sizeof *eof * sizeof eof) != 0) {
            if (fseek(s5p->fp, - sizeof *eof * sizeof eof, SEEK_CUR) != 0) { // Seek back
                return -1;
            }

            // Set start offset
            offset = ftello(s5p->fp);

            // Get record size
            slow5_rec_size_t record_size;
            if (fread(&record_size, sizeof record_size, 1, s5p->fp) != 1) {
                return -1;
            }

            size = sizeof record_size + record_size;

            uint8_t *read_comp = (uint8_t *) malloc(record_size);
            SLOW5_MALLOC_CHK(read_comp);
            if (fread(read_comp, record_size, 1, s5p->fp) != 1) {
                free(read_comp);
                return -1;
            }

            uint8_t *read_decomp = (uint8_t *) slow5_ptr_depress(s5p->compress, read_comp, record_size, NULL);
            if (read_decomp == NULL) {
                free(read_comp);
                free(read_decomp);
                return -1;
            }
            free(read_comp);

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
            slow5_idx_insert(index, read_id, offset, size);

            free(read_decomp);

            // Read in potential eof marker
            if (fread(buf_eof, sizeof *eof, sizeof eof, s5p->fp) != sizeof eof) {
                return -1;
            }
        }

        // Ensure actually at end of file
        if (fread(buf_eof, 1, 1, s5p->fp) != 0) {
            return -1;
        }
    }

    if (fseeko(s5p->fp, curr_offset, SEEK_SET != 0)) {
        return -1;
    }

    return 0;
}

void slow5_idx_write(struct slow5_idx *index) {

    //fprintf(index->fp, SLOW5_INDEX_HEADER);

    const char magic[] = SLOW5_INDEX_MAGIC_NUMBER;
    SLOW5_ASSERT(fwrite(magic, sizeof *magic, sizeof magic, index->fp) == sizeof magic);

    struct slow5_version version = SLOW5_INDEX_VERSION;
    SLOW5_ASSERT(fwrite(&version.major, sizeof version.major, 1, index->fp) == 1);
    SLOW5_ASSERT(fwrite(&version.minor, sizeof version.minor, 1, index->fp) == 1);
    SLOW5_ASSERT(fwrite(&version.patch, sizeof version.patch, 1, index->fp) == 1);

    uint8_t padding = SLOW5_INDEX_HEADER_SIZE_OFFSET -
            sizeof magic * sizeof *magic -
            sizeof version.major -
            sizeof version.minor -
            sizeof version.patch;
    uint8_t *zeroes = (uint8_t *) calloc(padding, sizeof *zeroes);
    SLOW5_ASSERT(fwrite(zeroes, sizeof *zeroes, padding, index->fp) == padding);
    free(zeroes);

    for (uint64_t i = 0; i < index->num_ids; ++ i) {

        khint_t pos = kh_get(slow5_s2i, index->hash, index->ids[i]);
        SLOW5_ASSERT(pos != kh_end(index->hash));

        struct slow5_rec_idx read_index = kh_value(index->hash, pos);

        /*
        SLOW5_ASSERT(fprintf(index->fp, "%s" SLOW5_SEP_COL "%" PRIu64 SLOW5_SEP_COL "%" PRIu64 "\n",
                index->ids[i],
                read_index.offset,
                read_index.size) >= 0);
        */
        slow5_rid_len_t read_id_len = strlen(index->ids[i]);
        SLOW5_ASSERT(fwrite(&read_id_len, sizeof read_id_len, 1, index->fp) == 1);
        SLOW5_ASSERT(fwrite(index->ids[i], sizeof *index->ids[i], read_id_len, index->fp) == read_id_len);
        SLOW5_ASSERT(fwrite(&read_index.offset, sizeof read_index.offset, 1, index->fp) == 1);
        SLOW5_ASSERT(fwrite(&read_index.size, sizeof read_index.size, 1, index->fp) == 1);
    }

    const char eof[] = SLOW5_INDEX_EOF;
    SLOW5_ASSERT(fwrite(eof, sizeof *eof, sizeof eof, index->fp) == sizeof eof);
}

static inline int slow5_idx_is_version_compatible(struct slow5_version file_version){

    struct slow5_version supported_max_version = SLOW5_INDEX_VERSION;

    if(file_version.major > supported_max_version.major){
        return 0;
    }
    else if (file_version.minor > supported_max_version.minor){
        return 0;
    }
    else if (file_version.patch > supported_max_version.patch){
        return 0;
    }
    else{
        return 1;
    }
}

static void slow5_idx_read(struct slow5_idx *index) {

    const char magic[] = SLOW5_INDEX_MAGIC_NUMBER;
    char buf_magic[sizeof magic]; // TODO is this a vla?
    SLOW5_ASSERT(fread(buf_magic, sizeof *magic, sizeof magic, index->fp) == sizeof magic);
    SLOW5_ASSERT(memcmp(magic, buf_magic, sizeof *magic * sizeof magic) == 0);

    SLOW5_ASSERT(fread(&index->version.major, sizeof index->version.major, 1, index->fp) == 1);
    SLOW5_ASSERT(fread(&index->version.minor, sizeof index->version.minor, 1, index->fp) == 1);
    SLOW5_ASSERT(fread(&index->version.patch, sizeof index->version.patch, 1, index->fp) == 1);

    if(slow5_idx_is_version_compatible(index->version)==0){
        struct slow5_version supported_max_version = SLOW5_INDEX_VERSION;
        SLOW5_ERROR("file version (%d.%d.%d) in your slow5 index file is higher than the maximally compatible version (%d.%d.%d) by this slow5lib. Please re-index or use a newer version of slow5lib",
                index->version.major, index->version.minor, index->version.patch,
                     supported_max_version.major,  supported_max_version.minor,  supported_max_version.patch);
        SLOW5_ASSERT(0);
    }

    SLOW5_ASSERT(fseek(index->fp, SLOW5_INDEX_HEADER_SIZE_OFFSET, SEEK_SET) != -1);

    const char eof[] = SLOW5_INDEX_EOF;
    char buf_eof[sizeof eof]; // TODO is this a vla?

    SLOW5_ASSERT(fread(buf_eof, sizeof *eof, sizeof eof, index->fp) == sizeof eof);
    while (memcmp(eof, buf_eof, sizeof *eof * sizeof eof) != 0) {
        SLOW5_ASSERT(fseek(index->fp, - sizeof *eof * sizeof eof, SEEK_CUR) == 0); // Seek back

        slow5_rid_len_t read_id_len;
        SLOW5_ASSERT(fread(&read_id_len, sizeof read_id_len, 1, index->fp) == 1);
        char *read_id = (char *) malloc((read_id_len + 1) * sizeof *read_id); // +1 for '\0'
        if(read_id == NULL){
            SLOW5_ERROR("%s","Reads_id returned was NULL");
        }

        SLOW5_ASSERT(fread(read_id, sizeof *read_id, read_id_len, index->fp) == read_id_len);
        read_id[read_id_len] = '\0'; // Add null byte

        uint64_t offset;
        uint64_t size;

        SLOW5_ASSERT(fread(&offset, sizeof offset, 1, index->fp) == 1);
        SLOW5_ASSERT(fread(&size, sizeof size, 1, index->fp) == 1);

        slow5_idx_insert(index, read_id, offset, size);

        SLOW5_ASSERT(fread(buf_eof, sizeof *eof, sizeof eof, index->fp) == sizeof eof);
    }
}

void slow5_idx_insert(struct slow5_idx *index, char *read_id, uint64_t offset, uint64_t size) {

    int absent;
    khint_t k = kh_put(slow5_s2i, index->hash, read_id, &absent);
    SLOW5_ASSERT(absent != -1);
    SLOW5_ASSERT(absent != 0); // TODO error if read_id duplicated?

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
}

// -1 if read_id not in the hash map
// 0 otherwise
int slow5_idx_get(struct slow5_idx *index, const char *read_id, struct slow5_rec_idx *read_index) {
    int ret = 0;

    khint_t pos = kh_get(slow5_s2i, index->hash, read_id);
    if (pos == kh_end(index->hash)) {
        ret = -1;
    } else {
        if (read_index != NULL) {
            *read_index = kh_value(index->hash, pos);
        }
    }

    return ret;
}

void slow5_idx_free(struct slow5_idx *index) {
    //NULL_CHK(index); // TODO necessary?

    if (index->fp != NULL) {
        SLOW5_ASSERT(fclose(index->fp) == 0);
    }

    for (uint64_t i = 0; i < index->num_ids; ++ i) {
        free(index->ids[i]);
    }
    free(index->ids);

    kh_destroy(slow5_s2i, index->hash);

    free(index->pathname);
    free(index);
}
