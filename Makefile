CC       = gcc
CXX      = g++
CFLAGS   += -g -rdynamic -Wall -O2 -std=c++11 

-include config.mk

LDFLAGS += $(LIBS) -lpthread -lz

SOURCE_DIR = src
BUILD_DIR = build

SRC = $(wildcard $(SOURCE_DIR)/*.c)
DEPS = $(wildcard $(SOURCE_DIR)/*.h)
OBJ = $(patsubst $(SOURCE_DIR)/%.c,$(BUILD_DIR)/%.o,$(SRC))
BINARY = f5c

HDF5 ?= install
HTS ?= install

HTS_VERSION = 1.9
HDF5_VERSION = 1.10.4

ifdef ENABLE_PROFILE
    CFLAGS += -p
endif

ifdef cuda
    SRC_CUDA = $(wildcard $(SOURCE_DIR)/*.cu)
    DEPS_CUDA = $(SOURCE_DIR)/f5c.h \
		$(SOURCE_DIR)/fast5lite.h \
		$(SOURCE_DIR)/error.h \
		$(SOURCE_DIR)/f5cmisc.cuh
    OBJ_CUDA = $(patsubst $(SOURCE_DIR)/%.cu,$(BUILD_DIR)/%_cuda.o,$(SRC_CUDA))
    CC_CUDA = nvcc
    CFLAGS_CUDA = -g -O2 -std=c++11 -lineinfo $(CUDA_ARCH)
    CUDALIB += -L/usr/local/cuda/lib64/ -lcudart -lcudadevrt
    CUDALIB_STATIC += -L/usr/local/cuda/lib64/ -lcudart_static -lcudadevrt -lrt
    CPPFLAGS += -DHAVE_CUDA=1
endif

ifeq ($(HDF5), install)
    HDF5_LIB = $(BUILD_DIR)/lib/libhdf5.a
    HDF5_INC = -I$(BUILD_DIR)/include
    LDFLAGS += $(HDF5_LIB) -ldl
else
ifneq ($(HDF5), autoconf)
    HDF5_LIB =
    HDF5_SYS_LIB = $(shell pkg-config --libs hdf5)
    HDF5_INC = $(shell pkg-config --cflags-only-I hdf5)
endif
endif

ifeq ($(HTS), install)
    HTS_LIB = $(BUILD_DIR)/lib/libhts.a
    HTS_INC = -I$(BUILD_DIR)/include
    LDFLAGS += $(HTS_LIB)
else
ifneq ($(HTS), autoconf)
    HTS_LIB =
    HTS_SYS_LIB = $(shell pkg-config --libs htslib)
    HTS_INC = $(shell pkg-config --cflags-only-I htslib)
endif	
endif

CPPFLAGS += $(HDF5_INC) $(HTS_INC)
LDFLAGS += $(HTS_SYS_LIB) $(HDF5_SYS_LIB)

.PHONY: clean distclean test

ifdef cuda
$(BINARY): $(HTS_LIB) $(HDF5_LIB) $(OBJ) $(BUILD_DIR)/gpucode.o $(OBJ_CUDA)
	$(CXX) $(CFLAGS) $(OBJ) $(BUILD_DIR)/gpucode.o $(OBJ_CUDA) $(LDFLAGS) $(CUDALIB) -o $@

$(BINARY)_static: $(HTS_LIB) $(HDF5_LIB) $(OBJ) $(BUILD_DIR)/gpucode.o $(OBJ_CUDA)
	$(CXX) -static $(CFLAGS) $(OBJ) $(BUILD_DIR)/gpucode.o $(OBJ_CUDA) $(CUDALIB_STATIC) $(LDFLAGS) -lsz -laec $^ -o $@
else
$(BINARY): $(HTS_LIB) $(HDF5_LIB) $(OBJ)
	$(CXX) $(CFLAGS) $(OBJ) $(LDFLAGS) $(CUDALIB) -o $@

$(BINARY)_static: $(HTS_LIB) $(HDF5_LIB) $(OBJ)
	$(CXX) -static $(CFLAGS) $(OBJ) $(CUDALIB_STATIC) $(LDFLAGS) -lsz -laec $^ -o $@
endif

$(OBJ): $(BUILD_DIR)/%.o: $(SOURCE_DIR)/%.c $(DEPS) $(SOURCE_DIR)/config.h
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

$(SOURCE_DIR)/config.h:
	echo "/* Default config.h generated by Makefile */" >> $@
	echo "#define HAVE_HDF5_H 1" >> $@

$(BUILD_DIR)/gpucode.o: $(OBJ_CUDA)
	@mkdir -p $(BUILD_DIR)
	$(CC_CUDA) $(CFLAGS_CUDA) -dlink $^ -o $@

$(OBJ_CUDA): $(BUILD_DIR)/%_cuda.o: $(SOURCE_DIR)/%.cu $(DEPS_CUDA)
	@mkdir -p $(BUILD_DIR)
	$(CC_CUDA) -x cu $(CFLAGS_CUDA) $(CPPFLAGS) $(HDF5_INC) $(HTS_INC) -rdc=true -c $< -o $@

$(BUILD_DIR)/lib/libhts.a:
	@mkdir -p $(BUILD_DIR)
	@if command -v curl; then \
		curl -o $(BUILD_DIR)/htslib.tar.bz2 -L https://github.com/samtools/htslib/releases/download/$(HTS_VERSION)/htslib-$(HTS_VERSION).tar.bz2; \
	else \
		wget -O $(BUILD_DIR)/htslib.tar.bz2 https://github.com/samtools/htslib/releases/download/$(HTS_VERSION)/htslib-$(HTS_VERSION).tar.bz2; \
	fi
	tar -xf $(BUILD_DIR)/htslib.tar.bz2 -C $(BUILD_DIR)
	mv $(BUILD_DIR)/htslib-$(HTS_VERSION) $(BUILD_DIR)/htslib
	$(RM) $(BUILD_DIR)/htslib.tar.bz2
	cd $(BUILD_DIR)/htslib && \
	./configure --prefix=$(shell pwd)/$(BUILD_DIR) --enable-bz2=no --enable-lzma=no --with-libdeflate=no --enable-libcurl=no  --enable-gcs=no --enable-s3=no && \
	make -j8 && \
	make install

$(BUILD_DIR)/lib/libhdf5.a:
	@mkdir -p $(BUILD_DIR)
	@if command -v curl; then \
		curl -o $(BUILD_DIR)/hdf5.tar.bz2 https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-$(shell echo $(HDF5_VERSION) | awk -F. '{print $$1"."$$2}')/hdf5-$(HDF5_VERSION)/src/hdf5-$(HDF5_VERSION).tar.bz2; \
	else \
		wget -O $(BUILD_DIR)/hdf5.tar.bz2 https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-$(shell echo $(HDF5_VERSION) | awk -F. '{print $$1"."$$2}')/hdf5-$(HDF5_VERSION)/src/hdf5-$(HDF5_VERSION).tar.bz2; \
	fi
	tar -xf $(BUILD_DIR)/hdf5.tar.bz2 -C $(BUILD_DIR)
	mv $(BUILD_DIR)/hdf5-$(HDF5_VERSION) $(BUILD_DIR)/hdf5
	$(RM) $(BUILD_DIR)/hdf5.tar.bz2
	cd $(BUILD_DIR)/hdf5 && \
	./configure --prefix=$(shell pwd)/$(BUILD_DIR) && \
	make -j8 && \
	make install

clean: 
	$(RM) -r $(BINARY)* $(BUILD_DIR)/*.o $(BUILD_DIR)/*.out

# Delete all gitignored files (but not directories)
distclean: clean
	git clean -f -X; rm -rf ./autom4te.cache
	$(RM) -r $(BUILD_DIR)
	$(RM) $(SOURCE_DIR)/config.h

test: $(BINARY)
	./scripts/test.sh
