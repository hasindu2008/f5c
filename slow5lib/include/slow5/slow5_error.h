/**
 * @file slow5_error.h
 * @brief SLOW5 error handling functions
 * @author Sasha Jenner (jenner.sasha@gmail.com), Hasindu Gamaarachchi (hasindu@garvan.org.au)
 * @date 27/02/2021
 */

/*
MIT License

Copyright (c) 2020 Sasha Jenner
Copyright (c) 2020 Hasindu Gamaarachchi

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


#ifndef SLOW5_ERROR_H
#define SLOW5_ERROR_H

#include <errno.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif /* _cplusplus */

/* This is for internal use only - do not use any of the following directly*/

/* Debug and verbosity */

// the level of verbosity in the log printed to the standard error
enum slow5_log_level_opt {
    SLOW5_LOG_OFF,      // nothing at all
    SLOW5_LOG_ERR,      // error messages
    SLOW5_LOG_WARN,     // warning and error messages
    SLOW5_LOG_INFO,     // information, warning and error messages
    SLOW5_LOG_VERB,     // verbose, information, warning and error messages
    SLOW5_LOG_DBUG      // debugging, verbose, information, warning and error messages
};

// determines if the library exit behaviour (if exit(EXIT_FAILURE) is called) on error
enum slow5_exit_condition_opt {
    // do not exit the programme even if an error occur
    //the user should have done error handling for each function call - typical library behaviour in C  (default at the moment)
    SLOW5_EXIT_OFF,
    // exit the programme if an error occurs (useful for a programmer who is lazy to handle errors)
    SLOW5_EXIT_ON_ERR,
    // exit the programme if a warning occurs. this is a bit silly, but incase someone wants
    SLOW5_EXIT_ON_WARN
};

/* slow5 errno */
int *slow5_errno_location(void);
#define slow5_errno (*slow5_errno_location())

#define SLOW5_DEBUG_PREFIX "[DEBUG] %s: " /* TODO function before debug */
#define SLOW5_VERBOSE_PREFIX "[INFO] %s: "
#define SLOW5_INFO_PREFIX "[%s::INFO]\033[1;34m "
#define SLOW5_WARNING_PREFIX "[%s::WARNING]\033[1;33m "
#define SLOW5_ERROR_PREFIX "[%s::ERROR]\033[1;31m "
#define SLOW5_NO_COLOUR "\033[0m"

#define SLOW5_LOG_DEBUG(msg, ...) { \
    if (slow5_log_level >= SLOW5_LOG_DBUG) { \
        fprintf(stderr, SLOW5_DEBUG_PREFIX msg \
                " At %s:%d\n", \
                __func__, __VA_ARGS__, __FILE__, __LINE__ - 1); \
    } \
}

#define SLOW5_VERBOSE(msg, ...) { \
    if (slow5_log_level >= SLOW5_LOG_VERB) { \
        fprintf(stderr, SLOW5_VERBOSE_PREFIX msg "\n", __func__, __VA_ARGS__); \
    } \
}

#define SLOW5_INFO(msg, ...) { \
    if (slow5_log_level >= SLOW5_LOG_INFO) { \
        fprintf(stderr, SLOW5_INFO_PREFIX msg SLOW5_NO_COLOUR "\n", __func__, __VA_ARGS__); \
    } \
}

#define SLOW5_WARNING(msg, ...) { \
    if (slow5_log_level >= SLOW5_LOG_WARN) { \
        fprintf(stderr, SLOW5_WARNING_PREFIX msg SLOW5_NO_COLOUR \
                " At %s:%d\n", \
                __func__, __VA_ARGS__, __FILE__, __LINE__ - 1); \
    } \
    if (slow5_exit_condition >= SLOW5_EXIT_ON_WARN){ \
        SLOW5_INFO("%s", "Exiting on warning."); \
        exit(EXIT_FAILURE); \
    } \
}

#define SLOW5_ERROR(msg, ...) { \
    if (slow5_log_level >= SLOW5_LOG_ERR) { \
        fprintf(stderr, SLOW5_ERROR_PREFIX msg SLOW5_NO_COLOUR \
                " At %s:%d\n", \
                __func__, __VA_ARGS__, __FILE__, __LINE__ - 1); \
    } \
}

#define SLOW5_ERROR_EXIT(msg, ...) { \
    SLOW5_ERROR(msg, __VA_ARGS__) \
    SLOW5_EXIT_IF_ON_ERR() \
}

#define SLOW5_EXIT_IF_ON_ERR() { \
    if (slow5_exit_condition >= SLOW5_EXIT_ON_ERR){ \
        SLOW5_ERROR("%s", "Exiting on error."); \
        exit(EXIT_FAILURE); \
    } \
}

#define SLOW5_MALLOC_CHK(ret) { \
    if ((ret) == NULL) { \
        SLOW5_MALLOC_ERROR() \
    } \
}

#define SLOW5_MALLOC_CHK_EXIT(ret) { \
    if ((ret) == NULL) { \
        SLOW5_MALLOC_ERROR_EXIT() \
    } \
}

#define SLOW5_MALLOC_ERROR() SLOW5_ERROR("Failed to allocate memory: %s", strerror(errno))
#define SLOW5_MALLOC_ERROR_EXIT() SLOW5_ERROR_EXIT("Failed to allocate memory: %s", strerror(errno))

#define SLOW5_ASSERT(ret) { \
    if ((ret) == 0){ \
        fprintf(stderr, SLOW5_ERROR_PREFIX "Assertion failed." SLOW5_NO_COLOUR \
                " At %s:%d\nExiting.\n", \
                __func__ , __FILE__, __LINE__ - 1); \
        exit(EXIT_FAILURE); \
    } \
}

#ifdef __cplusplus
}
#endif /* _cplusplus */

#endif /* slow5_error.h */
