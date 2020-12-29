#ifndef STRINGLITE_H

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdarg.h>
#include <assert.h>
#include "error.h"

//adapted from https://github.com/lh3/minimap2/blob/master/kseq.h
#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	uint64_t l, m;
	char *s;
} kstring_t;
#endif

//from https://github.com/lh3/minimap2/blob/master/kseq.h
#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

//from https://github.com/lh3/minimap2/blob/master/format.c
static inline void str_enlarge(kstring_t *s, int l)
{
	if (s->l + l + 1 > s->m) {
		s->m = s->l + l + 1;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
}

static inline void str_init(kstring_t *s, int size){
	s->l=0;
	s->m=kroundup32(size);
	s->s=(char *)malloc(sizeof(char)* s->m);
	MALLOC_CHK(s->s);
	s->s[0]='\0';
}

static inline void str_free(kstring_t *s){
	free(s->s);
}


static inline void str_cat(kstring_t *s, char *t, int t_l){
	str_enlarge(s,t_l);
	memcpy(&s->s[s->l], t, t_l);
	s->l += t_l;
	s->s[s->l] = '\0';
}


static inline void sprintf_append(kstring_t *s, const char *fmt, ... ){

	va_list ap;
	char buffer[10000];

	va_start(ap, fmt);
	int len = vsnprintf(buffer,10000,fmt,ap);
	va_end(ap);

	assert(len>=0);
	if(len>=10000){
		WARNING("Too long string got truncated: %s",buffer);
	}

	str_cat(s,buffer,len);

}

static inline void str_copy(kstring_t *s, const char *st, const char *en)
{
	str_enlarge(s, en - st);
	memcpy(&s->s[s->l], st, en - st);
	s->l += en - st;
}


// static void sprintf_lite(kstring_t *s, const char *fmt, ...)
// {
// 	char buf[20]; // for integer to string conversion
// 	const char *p, *q;
// 	va_list ap;
// 	va_start(ap, fmt);
// 	for (q = p = fmt; *p; ++p) {
// 		if (*p == '%') {
// 			if (p > q) str_copy(s, q, p);
// 			++p;
// 			if (*p == 'd') {
// 				int32_t c, i, l = 0;
// 				uint32_t x;
// 				c = va_arg(ap, int32_t);
// 				x = c >= 0? c : -c;
// 				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
// 				if (c < 0) buf[l++] = '-';
// 				str_enlarge(s, l);
// 				for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
// 			} else if (*p == 'u') {
// 				int32_t i, l = 0;
// 				uint32_t x;
// 				x = va_arg(ap, uint32_t);
// 				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
// 				str_enlarge(s, l);
// 				for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
// 			} else if (*p == 's') {
// 				char *r = va_arg(ap, char*);
// 				str_copy(s, r, r + strlen(r));
// 			} else if (*p == 'c') {
// 				str_enlarge(s, 1);
// 				s->s[s->l++] = va_arg(ap, int);
// 			} else if (*p == 'l') {
// 				++p;
// 				if (*p == 'd') {
// 					int64_t c, i, l = 0;
// 					uint64_t x;
// 					c = va_arg(ap, int64_t);
// 					x = c >= 0? c : -c;
// 					do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
// 					if (c < 0) buf[l++] = '-';
// 					str_enlarge(s, l);
// 					for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
// 				} else if (*p == 'u') {
// 					int64_t i, l = 0;
// 					uint64_t x;
// 					x = va_arg(ap, uint64_t);
// 					do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
// 					str_enlarge(s, l);
// 					for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
// 				} else abort();
// 			}
// 			else abort();
// 			q = p + 1;
// 		}
// 	}
// 	if (p > q) str_copy(s, q, p);
// 	va_end(ap);
// 	s->s[s->l] = 0;
// }


#endif
