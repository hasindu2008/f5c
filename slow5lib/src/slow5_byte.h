/**
 * @file slow5_byte.h
 * @brief SLOW5 byte handling functions
 * @author Hasindu Gamaarachchi (hasindu@garvan.org.au)
 * @date 05/05/2024
 */

/*
MIT License

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


#ifndef SLOW5_BYTE_H
#define SLOW5_BYTE_H

#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif /* _cplusplus */

/* This is for internal use only - do not use any of the following directly*/


static inline void slow5_byte_swap(void *dest, const void *src, size_t size){
    if(size==2){

        // uint16_t tmp = ((uint16_t*)src)[0];
        // tmp = (tmp << 8) | (tmp >> 8);
        // ((uint16_t*)dest)[0] = tmp;

        char tmp0 = ((char*)src)[0];
        char tmp1 = ((char*)src)[1];
        ((char*)dest)[0] = tmp1;
        ((char*)dest)[1] = tmp0;

    } else if(size==4){

        // uint32_t tmp = ((uint32_t*)src)[0];
        // tmp = (tmp >> 24) | ((tmp << 8) & 0x00FF0000) | ((tmp >> 8) & 0x0000FF00) | (tmp << 24);
        // ((uint32_t*)dest)[0] = tmp;

        char tmp0 = ((char*)src)[0];
        char tmp1 = ((char*)src)[1];
        char tmp2 = ((char*)src)[2];
        char tmp3 = ((char*)src)[3];
        ((char*)dest)[0] = tmp3;
        ((char*)dest)[1] = tmp2;
        ((char*)dest)[2] = tmp1;
        ((char*)dest)[3] = tmp0;

    } else if(size==8){

        // uint64_t tmp = ((uint64_t*)src)[0];
        // tmp = (tmp >> 56) | ((tmp << 40) & 0x00FF000000000000) | ((tmp << 24) & 0x0000FF0000000000) |
        //      ((tmp << 8) & 0x000000FF00000000) | ((tmp >> 8) & 0x00000000FF000000) |
        //      ((tmp >> 24) & 0x0000000000FF0000) | ((tmp >> 40) & 0x000000000000FF00) | (tmp << 56);
        // ((uint64_t*)dest)[0] = tmp;

        char tmp0 = ((char*)src)[0];
        char tmp1 = ((char*)src)[1];
        char tmp2 = ((char*)src)[2];
        char tmp3 = ((char*)src)[3];
        char tmp4 = ((char*)src)[4];
        char tmp5 = ((char*)src)[5];
        char tmp6 = ((char*)src)[6];
        char tmp7 = ((char*)src)[7];
        ((char*)dest)[0] = tmp7;
        ((char*)dest)[1] = tmp6;
        ((char*)dest)[2] = tmp5;
        ((char*)dest)[3] = tmp4;
        ((char*)dest)[4] = tmp3;
        ((char*)dest)[5] = tmp2;
        ((char*)dest)[6] = tmp1;
        ((char*)dest)[7] = tmp0;

    }
}


static inline int slow5_fwrite_bigend(const void *ptr, size_t size, size_t nitems, FILE *stream) {

    if(!(size==1 || size==2 || size==4 || size==8)){
        fprintf(stderr,"ERROR: slow5_fwrite_bigend is only implemented for data types of size %zu\n",size);
        return -1;
    }
    if(nitems!=1){
        fprintf(stderr,"ERROR: slow5_fwrite_bigend is only implemented for nitems=1\n");
        return -1;
    }

    if(size == 1){
        return fwrite(ptr, size, nitems, stream);
    } else {
        void *ptr_copy = malloc(size*nitems);
        if(ptr_copy==NULL){
            fprintf(stderr,"[%s::ERROR]\033[1;31m malloc failled\033[0m\n", __func__);
            fprintf(stderr,"At %s:%d\n",  __FILE__, __LINE__-2);
            return -1;
        }
        slow5_byte_swap(ptr_copy, ptr, size);
        int ret = fwrite(ptr_copy, size, nitems, stream);
        free(ptr_copy);
        return ret;
    }
}

static inline int slow5_fread_bigend(void *ptr, size_t size, size_t nitems, FILE *stream) {

    if(!(size==1 || size==2 || size==4 || size==8)){
        fprintf(stderr,"[%s::ERROR]\033[1;31m slow5_fwrite_bigend is only implemented for data types of size %zu\033[0m\n",__func__,size);
        fprintf(stderr,"At %s:%d\n", __FILE__, __LINE__-2);
        return -1;
    }
    if(nitems!=1){
        fprintf(stderr,"[%s::ERROR]\033[1;31m slow5_fwrite_bigend is only implemented for nitems=1\033[0m\n", __func__);
        fprintf(stderr,"At %s:%d\n",  __FILE__, __LINE__-2);
        return -1;
    }

    if(size == 1){
        return fread(ptr, size, nitems, stream);
    } else {
        void *ptr_copy = malloc(size*nitems);
        if(ptr_copy==NULL){
            fprintf(stderr,"[%s::ERROR]\033[1;31m malloc failled\033[0m\n", __func__);
            fprintf(stderr,"At %s:%d\n",  __FILE__, __LINE__-2);
            return -1;
        }
        size_t ret = fread(ptr_copy, size, nitems, stream);
        slow5_byte_swap(ptr, ptr_copy, size);
        free(ptr_copy);
        return ret;
    }
}

#define SLOW5_FWRITE(ptr, size, nitems, stream) \
    ( (slow5_bigend) ? (slow5_fwrite_bigend((ptr), (size), (nitems), (stream))) : (fwrite((ptr), (size), (nitems), (stream))) )

#define SLOW5_FREAD(ptr, size, nitems, stream) \
    ( (slow5_bigend) ? (slow5_fread_bigend((ptr), (size), (nitems), (stream))) : (fread((ptr), (size), (nitems), (stream))) )

static inline void *slow5_memcpy_bigend(void *dest, const void * src, size_t size) {
    if(!(size==1 || size==2 || size==4 || size==8)){
        fprintf(stderr,"[%s::ERROR]\033[1;31m slow5_fwrite_bigend is only implemented for data types of size %zu\033[0m\n",__func__,size);
        fprintf(stderr,"At %s:%d\n", __FILE__, __LINE__-2);
        exit(EXIT_FAILURE);
    }
    slow5_byte_swap((dest), (src), (size));
    return (dest);
}

#define SLOW5_MEMCPY(dest, src, size) \
    ( (slow5_bigend) ? (slow5_memcpy_bigend((dest), (src), (size))) : (memcpy((dest), (src), (size))) )


#define SLOW5_BYTE_SWAP(ptr) \
    if(slow5_bigend){ \
        int size = sizeof(*(ptr)); \
        slow5_byte_swap((ptr), (ptr), size); \
    }

#define SLOW5_BYTE_SWAP_VOID(ptr,size) \
    if(slow5_bigend){ \
        slow5_byte_swap((ptr), (ptr), size); \
    }

#define SLOW5_BYTE_SWAP_ARRAY(ptr, nitems) \
    if(slow5_bigend){ \
        int size = sizeof(*(ptr)); \
        for(size_t i=0; i<(nitems); i++){ \
            slow5_byte_swap((ptr)+i, (ptr)+i, size); \
        } \
    }

#define SLOW5_BYTE_SWAP_ARRAY_VOID(ptr, size, nitems) \
    if(slow5_bigend){ \
        if((size)!=2){ \
            fprintf(stderr,"[%s::ERROR]\033[1;31m Big Endian is not supported for the feature. Open an issue please.\033[0m\n", __func__); \
            fprintf(stderr,"At %s:%d\n",  __FILE__, __LINE__-2);  \
            exit(EXIT_FAILURE); \
        } \
        for(int64_t i=0; i<(nitems); i++){ \
            char *thisptr = ((char *)(ptr))+i*(size); \
            slow5_byte_swap(thisptr, thisptr, (size)); \
        } \
    }


#define SLOW5_ENDIAN_ERROR(MSG) \
    if(slow5_bigend){ \
        fprintf(stderr,"[%s::ERROR]\033[1;31m Big Endian is not supported for this" MSG "feature. Open an issue please.\033[0m\n", __func__); \
        fprintf(stderr,"At %s:%d\n",  __FILE__, __LINE__-2);  \
        exit(EXIT_FAILURE); \
    }\

#define SLOW5_ENDIAN_WARN \
    if(slow5_bigend){ \
        fprintf(stderr,"[%s::ERROR]\033[1;31m Big Endian is not supported for this multi-byte array feature. Don't trust any results. Open an issue please.\033[0m\n", __func__); \
        fprintf(stderr,"At %s:%d\n",  __FILE__, __LINE__);  \
    }\


#ifdef __cplusplus
}
#endif /* _cplusplus */

#endif /* slow5_byte.h */
