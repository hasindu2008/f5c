// Miscellaneous helper functions

#include <zlib.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdarg.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>
#include <inttypes.h>
#include <limits.h>
#include "slow5_misc.h"
#include <slow5/slow5_error.h>
#include <slow5/slow5_defs.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>

extern enum slow5_log_level_opt  slow5_log_level;
extern enum slow5_exit_condition_opt  slow5_exit_condition;

#define DBL_STRING_BUF_FIXED_CAP (64) // 2^6
#define FLT_STRING_BUF_FIXED_CAP (64) // 2^6

int slow5_uint_check(const char *str);
int slow5_int_check(const char *str);
int slow5_float_check(const char *str);

/*
 * From https://stackoverflow.com/questions/3774417/sprintf-with-automatic-memory-allocation
 * sprintf with automatic memory allocation and variable list argument
 * returns -1 on calloc error, which may be a result of vsnprintf returning a negative size from an invalid fmt ap combination
 * returns negative value if the 2nd vsnprintf fails
 */
int slow5_vasprintf(char **strp, const char *fmt, va_list ap) {
    va_list ap1;
    size_t size;
    char *buffer;

    va_copy(ap1, ap);
    size = vsnprintf(NULL, 0, fmt, ap1) + 1;
    va_end(ap1);
    buffer = (char *) calloc(1, size);
    SLOW5_MALLOC_CHK(buffer);

    if (!buffer)
        return -1;

    *strp = buffer;

    return vsnprintf(buffer, size, fmt, ap);
}

/*
 * From https://stackoverflow.com/questions/3774417/sprintf-with-automatic-memory-allocation
 * sprintf with automatic memory allocation
 * returns -1 on calloc/fmt error
 * returns negative value on other vsnprintf error
 */
int slow5_asprintf(char **strp, const char *fmt, ...) {
    int error;
    va_list ap;

    va_start(ap, fmt);
    error = slow5_vasprintf(strp, fmt, ap);
    va_end(ap);

    return error;
}


/*
 * from https://code.woboq.org/userspace/glibc/string/strsep.c.html
 * strsep source code
 * updates *stringp to after delim, sets *stringp to NULL if delim isn't found
 * returns pointer to where *stringp was before, NULL if *stringp is NULL
 * DONE
 */
char *
slow5_strsep (char **stringp, const char *delim)
{
  char *begin, *end;
  begin = *stringp;
  if (begin == NULL)
    return NULL;
  /* Find the end of the token.  */
  end = begin + strcspn (begin, delim);
  if (*end)
    {
      /* Terminate the token and set *STRINGP past NUL character.  */
      *end++ = '\0';
      *stringp = end;
    }
  else
    /* No more delimiters; this is the last token.  */
    *stringp = NULL;
  return begin;
}

int slow5_uint_check(const char *str) {

    // Check for:
    // empty string
    // first number is not 0 if more letters in string
    if (strlen(str) == 0 || (strlen(str) > 1 && str[0] == '0')) {
        return -1;
    }

    // Ensure only integers in string
    for (size_t i = 0; i < strlen(str); ++ i) {
        if (!isdigit(str[i])) {
            return -1;
        }
    }

    return 0;
}

int slow5_int_check(const char *str) {

    // Check for:
    // empty string
    // first number is not 0 if more letters in string
    if (strlen(str) == 0 || (strlen(str) > 1 && str[0] == '0')) {
        return -1;
    }

    // Ensure only integers in string
    for (size_t i = 0; i < strlen(str); ++ i) {
        if (!(isdigit(str[i]) || str[i] == '-')) {
            return -1;
        }
    }

    return 0;
}

int slow5_float_check(const char *str) {

    // Ensure no empty string
    if (strlen(str) == 0) {
        return -1;
    }

    // Ensure only integers in string, '.', '-'
    for (size_t i = 0; i < strlen(str); ++ i) {
        if (!(isdigit(str[i]) || str[i] == '.' || str[i] == '-')) {
            return -1;
        }
    }

    return 0;
}

// Atoi but to uint64_t
// and without any symbols
// and without 0 prefixing
uint64_t slow5_ato_uint64(const char *str, int *err) {
    uint64_t ret = 0;

    if (slow5_uint_check(str) == -1) {
        *err = -1;
        return ret;
    }

    unsigned long long int tmp = strtoull(str, NULL, 10);
    if (ret == ULLONG_MAX && errno == ERANGE) {
        *err = -1;
        return ret;
    }

    if (tmp > UINT64_MAX) {
        *err = -1;
        return ret;
    }

    ret = (uint64_t) tmp;
    *err = 0;
    return ret;
}

// Atoi but to uint32_t
// and without any symbols
// and without 0 prefixing
uint32_t slow5_ato_uint32(const char *str, int *err) {
    uint32_t ret = 0;
    if (slow5_uint_check(str) == -1) {
        *err = -1;
        return ret;
    }

    unsigned long int tmp = strtoul(str, NULL, 10);
    if (tmp == ULONG_MAX && errno == ERANGE) {
        *err = -1;
        return ret;
    }

    if (tmp > UINT32_MAX) {
        *err = -1;
        return ret;
    }

    ret = (uint32_t) tmp;
    *err = 0;
    return ret;
}

// Atoi but to uint16_t
// and without any symbols
// and without 0 prefixing
uint16_t slow5_ato_uint16(const char *str, int *err) {
    uint16_t ret = 0;
    if (slow5_uint_check(str) == -1) {
        *err = -1;
        return ret;
    }

    unsigned long int tmp = strtoul(str, NULL, 10);
    if (tmp == ULONG_MAX && errno == ERANGE) {
        *err = -1;
        return ret;
    }

    if (tmp > UINT16_MAX) {
        *err = -1;
        return ret;
    }

    ret = (uint16_t) tmp;
    *err = 0;
    return ret;
}

// Atoi but to uint8_t
// and without any symbols
// and without 0 prefixing
uint8_t slow5_ato_uint8(const char *str, int *err) {
    uint8_t ret = 0;
    if (slow5_uint_check(str) == -1) {
        *err = -1;
        return ret;
    }

    unsigned long int tmp = strtoul(str, NULL, 10);
    if (tmp > UINT8_MAX) {
        *err = -1;
        return ret;
    }

    ret = (uint8_t) tmp;
    *err = 0;
    return ret;
}

// Atoi but to int64_t
// and without any symbols
// and without 0 prefixing
int64_t slow5_ato_int64(const char *str, int *err) {
    int64_t ret = 0;
    if (slow5_int_check(str) == -1) {
        *err = -1;
        return ret;
    }

    long int tmp = strtol(str, NULL, 10);
    if (tmp > INT64_MAX || tmp < INT64_MIN) {
        *err = -1;
        return ret;
    }

    ret = (int64_t) tmp;
    *err = 0;
    return ret;
}

// Atoi but to int32_t
// and without any symbols
// and without 0 prefixing
int32_t slow5_ato_int32(const char *str, int *err) {
    int32_t ret = 0;
    if (slow5_int_check(str) == -1) {
        *err = -1;
        return ret;
    }

    long int tmp = strtol(str, NULL, 10);
    if (tmp > INT32_MAX || tmp < INT32_MIN) {
        *err = -1;
        return ret;
    }

    ret = (int32_t) tmp;
    *err = 0;
    return ret;
}

// Atoi but to int16_t
// and without any symbols
// and without 0 prefixing
int16_t slow5_ato_int16(const char *str, int *err) {
    int16_t ret = 0;
    if (slow5_int_check(str) == -1) {
        *err = -1;
        return ret;
    }

    long int tmp = strtol(str, NULL, 10);
    if (tmp > INT16_MAX || tmp < INT16_MIN) {
        *err = -1;
        return ret;
    }

    ret = (int16_t) tmp;
    *err = 0;
    return ret;
}

// Atoi but to int8_t
// and without any symbols
// and without 0 prefixing
int8_t slow5_ato_int8(const char *str, int *err) {
    int8_t ret = 0;
    if (slow5_int_check(str) == -1) {
        *err = -1;
        return ret;
    }

    long int tmp = strtol(str, NULL, 10);
    if (tmp > INT8_MAX || tmp < INT8_MIN) {
        *err = -1;
        return ret;
    }

    ret = (int8_t) tmp;
    *err = 0;
    return ret;
}

float slow5_strtof_check(const char *str, int *err) {
    float ret = 0;

    if (slow5_float_check(str) == -1) {
        *err = -1;
        return ret;
    }

    ret = strtof(str, NULL);
    if (errno == ERANGE && (ret == HUGE_VAL || ret == -HUGE_VAL || ret == (float) 0)) {
        *err = -1;
        return ret;
    }

    *err = 0;
    return ret;
}

double slow5_strtod_check(const char *str, int *err) {
    double ret = 0;

    if (slow5_float_check(str) == -1) {
        *err = -1;
        return ret;
    }

    ret = strtod(str, NULL);
    if (errno == ERANGE && (ret == HUGE_VAL || ret == -HUGE_VAL || ret == (double) 0)) {
        *err = -1;
        return ret;
    }

    *err = 0;
    return ret;
}

// Convert double to decimal string without trailing 0s
char *slow5_double_to_str(double x, size_t *len) {
    char *str = NULL;
    int n = slow5_asprintf(&str, "%f", x);

    for (int i = n - 1; i >= 1; -- i) {
        if (str[i] == '.') {
            n = i;
            str[n] = '\0';
            if (strcmp(str, "-0") == 0) {
                strcpy(str, "0"); /* remove minus sign */
                -- n;
            }
            break;
        } else if (str[i] != '0') {
            if (i != n - 1) {
                n = i + 1;
                str[n] = '\0';
            }
            break;
        }
    }

    if (len) {
        *len = n;
    }

    return str;
}

/*
static uint8_t get_type_size(const char *name) {
    uint8_t size = 0;

    if SLOW5_CHECK_TYPE(name, int8_t)
    else if SLOW5_CHECK_TYPE(name, uint8_t)
    else if SLOW5_CHECK_TYPE(name, int16_t)
    else if SLOW5_CHECK_TYPE(name, uint16_t)
    else if SLOW5_CHECK_TYPE(name, int32_t)
    else if SLOW5_CHECK_TYPE(name, uint32_t)
    else if SLOW5_CHECK_TYPE(name, int64_t)
    else if SLOW5_CHECK_TYPE(name, uint64_t)
    else if SLOW5_CHECK_TYPE(name, float)
    else if SLOW5_CHECK_TYPE(name, double)
    else if SLOW5_CHECK_TYPE(name, char)

    return size;
}
*/

/*
 * given filepaths a and b get the difference between their last modification times
 * return <a_time> - <b_time>
 * =0 both are created at the same time
 * >0 a is newer than b
 * <0 a is older than b
 * *err=-1 return 0 if failed
 * *err=0 return diff if success
 */
double slow5_filestamps_cmp(const char *a, const char *b, int *err) {
    struct stat a_stat;
    struct stat b_stat;
    if (stat(a, &a_stat) == -1) {
        SLOW5_ERROR("Failed to retrieve stats about file '%s'.", a);
        if (err) {
            *err = -1;
        }
        return 0;
    }
    if (stat(b, &b_stat) == -1) {
        SLOW5_ERROR("Failed to retrieve stats about file '%s'.", b);
        if (err) {
            *err = -1;
        }
        return 0;
    }

    if (err) {
        *err = 0;
    }
    return difftime(a_stat.st_mtime, b_stat.st_mtime);
}

/**
 * check if label is a valid c label
 * doesn't check if it is a reserved key word
 * label cannot be NULL
 * return 0 if valid
 * -1 if empty
 * -2 if contains a character that's not alphanumeric or a _
 * -3 if first character is a number
 */
int slow5_is_c_label(const char *label) {
    size_t len = strlen(label);

    if (len == 0) {
        return -1;
    }

    for (size_t i = 0; i < len; ++ i) {
        if (!isalnum(label[i]) && label[i] != '_') {
            return -2;
        }
    }

    if (isdigit(label[0])) {
        return -3;
    }

    return 0;
}
