# zstd compression is not available by default
# run `make zstd=1` to compile with zstd
# or uncomment the following line
#zstd=1

CC			= cc
AR			= ar
SVB			= thirdparty/streamvbyte
SVBLIB		= $(SVB)/libstreamvbyte.a
CPPFLAGS	+= -I include/ -I $(SVB)/include/
CFLAGS		+= -g -Wall -O2 -std=c99
LDFLAGS		+= -lm -lz
ifeq ($(zstd),1)
CFLAGS		+= -DSLOW5_USE_ZSTD
LDFLAGS		+= -lzstd
endif
ifeq ($(zstd_local),)
else
CFLAGS		+= -DSLOW5_USE_ZSTD
CPPFLAGS 	+= -I $(zstd_local)
endif
BUILD_DIR	= lib

STATICLIB	= $(BUILD_DIR)/libslow5.a
SHAREDLIB	= $(BUILD_DIR)/libslow5.so

OBJ = $(BUILD_DIR)/slow5.o \
		$(BUILD_DIR)/slow5_idx.o \
		$(BUILD_DIR)/slow5_misc.o \
		$(BUILD_DIR)/slow5_press.o \

PREFIX = /usr/local
VERSION = `git describe --tags`

SLOW5_H = include/slow5/slow5.h include/slow5/klib/khash.h include/slow5/klib/kvec.h include/slow5/slow5_defs.h include/slow5/slow5_error.h include/slow5/slow5_press.h

.PHONY: clean distclean test install uninstall slow5lib

#libslow5
slow5lib: $(SHAREDLIB) $(STATICLIB)

$(STATICLIB): $(OBJ) $(SVBLIB)
	cp $(SVBLIB) $@
	$(AR) rcs $@ $(OBJ)

$(SHAREDLIB): $(OBJ) $(SVBLIB)
	$(CC) $(CFLAGS) -shared $^ -o $@ $(LDFLAGS)

$(SVBLIB):
	make -C $(SVB) no_simd=$(no_simd) libstreamvbyte.a

$(BUILD_DIR)/slow5.o: src/slow5.c src/slow5_extra.h src/slow5_idx.h src/slow5_misc.h src/klib/ksort.h $(SLOW5_H)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -fpic -o $@

$(BUILD_DIR)/slow5_idx.o: src/slow5_idx.c src/slow5_idx.h src/slow5_extra.h src/slow5_misc.h $(SLOW5_H)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -fpic -o $@

$(BUILD_DIR)/slow5_misc.o: src/slow5_misc.c src/slow5_misc.h include/slow5/slow5_error.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -fpic -o $@

$(BUILD_DIR)/slow5_press.o: src/slow5_press.c include/slow5/slow5_press.h src/slow5_misc.h include/slow5/slow5_error.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -fpic -o $@

clean:
	rm -rf $(OBJ) $(STATICLIB) $(SHAREDLIB) $(SHAREDLIBV)
	make -C $(SVB) clean

# Delete all gitignored files (but not directories)
distclean: clean
	git clean -f -X
	rm -rf $(BUILD_DIR)/* autom4te.cache

test: slow5lib
	make -C test clean
	make -C test zstd=$(zstd)
	./test/test.sh

pyslow5:
	make clean
	rm -rf *.so python/pyslow5.cpp python/pyslow5.c build/lib.* build/temp.* build/bdist.* sdist pyslow5.egg-info dist
	python3 setup.py build
	cp build/lib.*/*.so  ./
	python3 < python/example.py
	python3 setup.py sdist


test-prep: slow5lib
	gcc test/make_blow5.c -Isrc src/slow5.c src/slow5_press.c -lm -lz src/slow5_idx.c src/slow5_misc.c -o test/bin/make_blow5 -g
	./test/bin/make_blow5

valgrind: slow5lib
	make -C test zstd=$(zstd)
	./test/test.sh mem

examples: slow5lib
	./examples/build.sh
