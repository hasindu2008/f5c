#ifndef SLOW5_EXTRA_H
#define SLOW5_EXTRA_H

#include <dirent.h>
#include <slow5/slow5.h>

#ifdef __cplusplus
extern "C" {
#endif

/**************************************************************************************************
 ***  Ultra-Low-level API ******************************************************************************
 **************************************************************************************************/

/*
These functions are not to be exposed to the public
Used for slow5tools
Any change to a function prototype here means slow5tools must be fixed
*/

// slow5 file
struct slow5_file *slow5_init(FILE *fp, const char *pathname, enum slow5_fmt format);
struct slow5_file *slow5_init_empty(FILE *fp, const char *pathname, enum slow5_fmt format);

// slow5 header
struct slow5_hdr *slow5_hdr_init_empty(void);
struct slow5_hdr *slow5_hdr_init(FILE *fp, enum slow5_fmt format, slow5_press_method_t *method);
void slow5_hdr_free(struct slow5_hdr *header);

// slow5 header data
int slow5_hdr_data_init(FILE *fp, char *buf, size_t *cap, struct slow5_hdr *header, uint32_t *hdr_len);
khash_t(slow5_s2s) *slow5_hdr_get_data(uint32_t read_group, const struct slow5_hdr *header);
int64_t slow5_hdr_add_rg_data(struct slow5_hdr *header, khash_t(slow5_s2s) *new_data);
char *slow5_hdr_types_to_str(struct slow5_aux_meta *aux_meta, size_t *len);
char *slow5_hdr_attrs_to_str(struct slow5_aux_meta *aux_meta, size_t *len);
void slow5_hdr_data_free(struct slow5_hdr *header);

struct slow5_aux_meta *slow5_aux_meta_init_empty(void);
struct slow5_aux_meta *slow5_aux_meta_init(FILE *fp, char *buf, size_t *cap, uint32_t *hdr_len);
int slow5_aux_meta_add(struct slow5_aux_meta *aux_meta, const char *attr, enum slow5_aux_type type);
void slow5_aux_meta_free(struct slow5_aux_meta *aux_meta);

// slow5 record
int slow5_rec_set(struct slow5_rec *read, struct slow5_aux_meta *aux_meta, const char *attr, const void *data);
int slow5_rec_set_array(struct slow5_rec *read, struct slow5_aux_meta *aux_meta, const char *attr, const void *data, size_t len);
static inline int slow5_rec_set_string(struct slow5_rec *read, struct slow5_aux_meta *aux_meta, const char *attr, const void *data) {
    return slow5_rec_set_array(read, aux_meta, attr, data, strlen((const char *) data) + 1);
}
int slow5_rec_parse(char *read_mem, size_t read_size, const char *read_id, struct slow5_rec *read, enum slow5_fmt format, struct slow5_aux_meta *aux_meta);
void slow5_rec_aux_free(khash_t(slow5_s2a) *aux_map);

// slow5 extension parsing
enum slow5_fmt slow5_name_get_fmt(const char *name);
enum slow5_fmt slow5_path_get_fmt(const char *path);
const char *slow5_fmt_get_name(enum slow5_fmt format);
char *slow5_get_idx_path(const char *path);

// auxilary type
enum slow5_aux_type slow5_str_to_aux_type(const char *str, int *err);
int slow5_memcpy_type_from_str(uint8_t *data, const char *value, enum slow5_aux_type type);
void slow5_memcpy_type_from_null_str(uint8_t *data, enum slow5_aux_type type);
char *slow5_type_to_str(uint8_t *data, const char *type, size_t len, size_t *str_len);
char *slow5_aux_type_to_str(enum slow5_aux_type type);
char *slow5_data_to_str(uint8_t *data, enum slow5_aux_type type, uint64_t len, size_t *str_len);











/* SLOW5 Extra API */

/*
// Convert fast5 files to a slow5 file
// TODO decide return type
int8_t fast5_to_slow5(const char *pathname_from, struct slow5_file *s5p_to);
// Convert slow5 file to fast5 files
int8_t slow5_to_fast5(struct slow5_file *s5p_from, const char *pathname_to);


// Header
int8_t slow5_hdr_write(struct slow5_file *s5p_from, struct slow5_file *s5p_to);
int8_t slow5_hdr_data_write(struct slow5_file *s5p_from, struct slow5_file *s5p_to);
int8_t slow5_hdr_data_attr_write(const char *attr, struct slow5_file *s5p_from, struct slow5_file *s5p_to);

// Index
inline khash_t(slow5_s2i) *slow5_idx_init_empty(void);

// Convert fast5 dir/file to a slow5 file
// TODO decide return type
int8_t fast5dir_to_slow5(DIR *dirp_from, struct slow5_file *s5p_to);
int8_t fast5fp_to_slow5(fast5_t f5p_from, struct slow5_file *s5p_to);

// Convert slow5 file to fast5 dir/file
// TODO decide return type
int8_t slow5_to_fast5dir(struct slow5_file *s5p_from, DIR *dirp_to);
int8_t slow5_to_fast5fp(struct slow5_file *s5p_from, fast5_t *f5p_to);

// Merge 2 slow5 files to another slow5 file
int8_t slow5_merge_2(struct slow5_file *s5p_from_1, struct slow5_file *s5p_from_2, struct slow5_file *s5p_to);

// Split a slow5 file to a dir
int8_t slow5_split(struct slow5_file *s5p, DIR *dirp);


// Initiate an empty slow5 read object
struct slow5_read *slow5_read_init_empty(void);
// Initiate a slow5 read
struct slow5_read *slow5_read_init(void); // TODO change void


// Print out the SLOW5 structure contents
void slow5_hdr_data_print(const struct SLOW5Header *hdr);


// Get the format from a slow5 format name
enum slow5_format str_get_slow5_format(const char *str);
// Get the format of a slow5 pathname
enum slow5_format path_get_slow5_format(const char *pathname);
// Get the format of a slow5 FILE
enum slow5_format stream_get_slow5_format(const FILE *stream);

// Get the slow5 format name from the format
const char *slow5_format_get_str(enum slow5_format format);

// Get the slow5 version array from a version string
//const uint8_t *str_get_slow5_version(const char *str);
*/


#ifdef __cplusplus
}
#endif


#endif
