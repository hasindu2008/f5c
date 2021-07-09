CC			= gcc
AR			= ar
CPPFLAGS	+= -I include/
CFLAGS		+= -g -Wall -O2 -std=c99
LDFLAGS		+=  -lz
BUILD_DIR	= lib

OBJ_LIB = $(BUILD_DIR)/slow5.o \
		$(BUILD_DIR)/slow5_idx.o	\
		$(BUILD_DIR)/slow5_misc.o	\
		$(BUILD_DIR)/slow5_press.o \

PREFIX = /usr/local
VERSION = `git describe --tags`

SLOW5_H = include/slow5/slow5.h include/slow5/klib/khash.h include/slow5/klib/kvec.h include/slow5/slow5_defs.h include/slow5/slow5_error.h include/slow5/slow5_press.h

.PHONY: clean distclean test install uninstall slow5lib

#libslow5
slow5lib: $(BUILD_DIR)/libslow5.so $(BUILD_DIR)/libslow5.a

$(BUILD_DIR)/libslow5.so: $(OBJ_LIB)
	$(CC) $(CFLAGS) -shared $^  -o $@ $(LDFLAGS)

$(BUILD_DIR)/libslow5.a: $(OBJ_LIB)
	$(AR) rcs $@ $^

$(BUILD_DIR)/slow5.o: src/slow5.c src/slow5_extra.h src/slow5_idx.h src/slow5_misc.h src/klib/ksort.h $(SLOW5_H)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -fpic -o $@

$(BUILD_DIR)/slow5_idx.o: src/slow5_idx.c src/slow5_idx.h src/slow5_extra.h src/slow5_misc.h $(SLOW5_H)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -fpic -o $@

$(BUILD_DIR)/slow5_misc.o: src/slow5_misc.c src/slow5_misc.h include/slow5/slow5_error.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -fpic -o $@

$(BUILD_DIR)/slow5_press.o: src/slow5_press.c include/slow5/slow5_press.h src/slow5_misc.h include/slow5/slow5_error.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -fpic -o $@

clean:
	rm -rf $(BUILD_DIR)/*.o $(BUILD_DIR)/libslow5.so $(BUILD_DIR)/libslow5.a

# Delete all gitignored files (but not directories)
distclean: clean
	git clean -f -X
	rm -rf $(BUILD_DIR)/* autom4te.cache

dist: distclean
	mkdir -p slow5lib-$(VERSION)
	cp -r README.md LICENSE Makefile configure.ac config.mk.in \
		installdeps.mk src docs build configure slow5tools-$(VERSION)
	#mkdir -p slow5lib-$(VERSION)/scripts
	#cp scripts/install-hdf5.sh slow5lib-$(VERSION)/scripts
	tar -zcf slow5lib-$(VERSION)-release.tar.gz slow5lib-$(VERSION)
	rm -rf slow5lib-$(VERSION)

binary:
	mkdir -p slow5lib-$(VERSION)
	make clean
	make && mv slow5lib slow5lib-$(VERSION)/slow5lib_x86_64_linux
	cp -r README.md LICENSE docs slow5lib-$(VERSION)/
	#mkdir -p slow5lib-$(VERSION)/scripts
	#cp scripts/test.sh slow5lib-$(VERSION)/scripts
	tar -zcf slow5lib-$(VERSION)-binaries.tar.gz slow5lib-$(VERSION)
	rm -rf slow5lib-$(VERSION)

install: slow5lib
	mkdir -p $(DESTDIR)$(PREFIX)/bin
	mkdir -p $(DESTDIR)$(PREFIX)/share/man/man1
	cp -f $(BINARY) $(DESTDIR)$(PREFIX)/bin
	gzip < docs/slow5lib.1 > $(DESTDIR)$(PREFIX)/share/man/man1/slow5lib.1.gz

uninstall:
	rm -f $(DESTDIR)$(PREFIX)/bin/$(BINARY) \
		$(DESTDIR)$(PREFIX)/share/man/man1/slow5lib.1.gz

test: slow5lib
	./test/test.sh

pyslow5:
	make clean
	rm -rf *.so python/pyslow5.cpp build/lib.* build/temp.*
	python3 setup.py build
	cp build/lib.*/*.so  ./
	python3 < python/example.py

test-prep: slow5lib
	gcc test/make_blow5.c -Isrc src/slow5.c src/slow5_press.c -lz src/slow5_idx.c src/slow5_misc.c -o test/bin/make_blow5 -g
	./test/bin/make_blow5

valgrind: slow5lib
	./test/test.sh mem

examples: slow5lib
	./examples/build.sh
