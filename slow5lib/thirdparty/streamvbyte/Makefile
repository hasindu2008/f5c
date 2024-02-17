# minimalist makefile
# NOTE: This is a modified version of the original Makefile.
# Only the files required for slow5lib are compiled by this makefile. The rest are commented out or removed.

.SUFFIXES:
#
.SUFFIXES: .cpp .o .c .h

PROC_AUTO:=$(shell uname -m)
PROCESSOR ?= $(PROC_AUTO)

ifeq ($(no_simd),1) #if no_simd is set
CFLAGS = -fPIC -std=c99 -O3 -Wall -Wextra -pedantic -Wshadow
else
ifeq ($(PROCESSOR), aarch64)
# for 64-bit ARM processors
CFLAGS = -fPIC -std=c99 -O3 -Wall -Wextra -pedantic -Wshadow -D__ARM_NEON__
else ifeq ($(PROCESSOR), armv7l)
# for 32-bit ARM processors
CFLAGS = -fPIC -std=c99 -O3 -Wall -Wextra -pedantic -Wshadow -mfpu=neon
else ifeq ($(PROCESSOR), x86_64)
CFLAGS = -fPIC -std=c99 -O3 -Wall -Wextra -pedantic -Wshadow -DSTREAMVBYTE_SSSE3=1 -mssse3
else
CFLAGS = -fPIC -std=c99 -O3 -Wall -Wextra -pedantic -Wshadow
endif
endif
# note that on macos, it is arm64 instead of aarch64. however, for now we just let it compile without simd

#LDFLAGS = -shared
LIBNAME=libstreamvbyte.so.0.0.1
LNLIBNAME=libstreamvbyte.so
STATICLIB=libstreamvbyte.a
# all:  $(LIBNAME) $(STATICLIB)
all: $(STATICLIB)

#SVB_HEADERS=./include/streamvbyte.h ./include/streamvbytedelta.h ./include/streamvbyte_zigzag.h
SVB_HEADERS=./include/streamvbyte.h ./include/streamvbyte_zigzag.h

#SVB_OBJECTS= streamvbyte_decode.o streamvbyte_encode.o streamvbytedelta_decode.o streamvbytedelta_encode.o streamvbyte_0124_encode.o  streamvbyte_0124_decode.o streamvbyte_zigzag.o
SVB_OBJECTS= streamvbyte_decode.o streamvbyte_encode.o streamvbyte_zigzag.o

streamvbyte_zigzag.o: ./src/streamvbyte_zigzag.c $(SVB_HEADERS)
	$(CC) $(CFLAGS) -c ./src/streamvbyte_zigzag.c -Iinclude


# streamvbytedelta_encode.o: ./src/streamvbytedelta_encode.c $(SVB_HEADERS)
# 	$(CC) $(CFLAGS) -c ./src/streamvbytedelta_encode.c -Iinclude

# streamvbytedelta_decode.o: ./src/streamvbytedelta_decode.c $(SVB_HEADERS)
# 	$(CC) $(CFLAGS) -c ./src/streamvbytedelta_decode.c -Iinclude

# streamvbyte_0124_encode.o: ./src/streamvbyte_0124_encode.c $(SVB_HEADERS)
# 	$(CC) $(CFLAGS) -c ./src/streamvbyte_0124_encode.c -Iinclude

# streamvbyte_0124_decode.o: ./src/streamvbyte_0124_decode.c $(SVB_HEADERS)
# 	$(CC) $(CFLAGS) -c ./src/streamvbyte_0124_decode.c -Iinclude

streamvbyte_decode.o: ./src/streamvbyte_decode.c $(SVB_HEADERS)
	$(CC) $(CFLAGS) -c ./src/streamvbyte_decode.c -Iinclude

streamvbyte_encode.o: ./src/streamvbyte_encode.c $(SVB_HEADERS)
	$(CC) $(CFLAGS) -c ./src/streamvbyte_encode.c -Iinclude

# $(LIBNAME): $(SVB_OBJECTS)
# 	$(CC) $(CFLAGS) -o $(LIBNAME) $(SVB_OBJECTS)  $(LDFLAGS)

$(STATICLIB): $(SVB_OBJECTS)
	ar rcs $@ $^

# $(LNLIBNAME): $(LIBNAME)
# 	ln -f -s $(LIBNAME) $(LNLIBNAME)


clean:
	rm -f unit *.o $(LIBNAME) $(LNLIBNAME) $(STATICLIB)
