CC       = g++
CFLAGS   = -g -rdynamic -Wall -O2 -std=c++11 

include config.mk

LDFLAGS += $(LIBS) -lpthread -lz -ldl

SRC = main.c f5c.c events.c nanopolish_read_db.c model.c align.c meth.c hmm.c
OBJ = $(SRC:.c=.o)
BINARY = f5c
DEPS = config.h  error.h  f5c.h  f5cmisc.h  fast5lite.h  logsum.h  matrix.h  model.h  nanopolish_read_db.h

HDF5 ?= install
HTS ?= install

BUILD_DIR = $(shell pwd)/build

ifneq ($(cuda),) # if cuda is undefined
    DEPS_CUDA = f5c.h fast5lite.h error.h f5cmisc.cuh
    SRC_CUDA = f5c.cu align.cu 
    OBJ_CUDA = $(SRC_CUDA:.cu=_cuda.o)
    CC_CUDA = nvcc
    CFLAGS_CUDA = -g  -O2 -std=c++11 -lineinfo -DHAVE_CUDA=1 $(CUDA_ARCH)
    LDFLAGS += -L/usr/local/cuda/lib64/ -lcudart -lcudadevrt
    OBJ += gpucode.o $(OBJ_CUDA)
    CFLAGS += -DHAVE_CUDA=1
endif

ifeq ($(HDF5), install)
    HDF5_LIB = $(BUILD_DIR)/lib/libhdf5.a
    HDF5_INC = -I$(BUILD_DIR)/include
else
    HDF5_LIB =
    HDF5_INC =
    LIBS += -lhdf5
endif

ifeq ($(HTS), install)
    HTS_LIB = $(BUILD_DIR)/lib/libhts.a
    HTS_INC = -I$(BUILD_DIR)/include
else
    HTS_LIB =
    HTS_INC =
    LIBS += -lhts
endif

CFLAGS += $(HDF5_INC) $(HTS_INC)

.PHONY: clean distclean format test

$(BINARY): $(HTS_LIB) $(HDF5_LIB) $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) $(HTS_LIB) $(HDF5_LIB) $(LDFLAGS) -o $@

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $(HDF5_INC) $(HTS_INC) $< -c 

gpucode.o: $(OBJ_CUDA)
	$(CC_CUDA) $(CFLAGS_CUDA) -dlink $^ -o $@

%_cuda.o: %.cu $(DEPS_CUDA)
	$(CC_CUDA) -x cu $(CFLAGS_CUDA) $(CPPFLAGS) $(HDF5_INC) $(HTS_INC) -rdc=true -c $< -o $@

external/htslib/configure:
	git submodule update --recursive --init --remote
	cd external/htslib && autoreconf

external/hdf5/configure:
	git submodule update --recursive --init --remote
	cd external/hdf5 && autoreconf -i

$(BUILD_DIR)/lib/libhts.a: external/htslib/configure
	cd external/htslib && \
	./configure --prefix=$(BUILD_DIR) --enable-bz2=no --enable-lzma=no --with-libdeflate=no --enable-libcurl=no  --enable-gcs=no --enable-s3=no && \
	make -j8 && \
	make install

$(BUILD_DIR)/lib/libhdf5.a: external/hdf5/configure
	cd external/hdf5 && \
	./configure --prefix=$(BUILD_DIR) && \
	make -j8 && \
	make install

clean: 
	rm -rf f5c *.o *.out $(BUILD_DIR)

# Delete all gitignored files (but not directories)
distclean: clean
	git clean -f -X; rm -rf ./autom4te.cache

# Autoformat code with clang format
format:
	./scripts/autoformat.sh

test: $(BINARY)
	./scripts/test.sh

valgrind : $(BINARY)
	./scripts/test.sh valgrind

kepler :
	rsync -av *.cu *.cuh Makefile $(SRC) $(DEPS) hasindu@kepler:/storage/hasindu/f5c/ && ssh kepler 'cd /storage/hasindu/f5c/ && make cuda=1'
	rsync -av scripts/*.sh hasindu@kepler:/storage/hasindu/f5c/scripts/ 

jetson:
	rsync -av *.cu *.cuh Makefile $(SRC) $(DEPS) hasindu@jetson:~/f5c/ && ssh jetson 'cd ~/f5c/ && make cuda=1'
	rsync -av scripts/*.sh hasindu@jetson:~/f5c/scripts/

nsight:
	#nvprof  -f --kernels "align_kernel_core" --analysis-metrics -o bv.nvprof ./f5c -b test/chr22_meth_example/reads10k.bam -g test/chr22_meth_example//humangenome.fa -r test/chr22_meth_example//reads10k.fq -t 8 --print-scaling=yes -K512 --cuda-block-size=64 --debug-break=yes > /dev/null
	nvprof  -f --analysis-metrics -o bv.nvprof ./f5c -b test/chr22_meth_example/reads10k.bam -g test/chr22_meth_example//humangenome.fa -r test/chr22_meth_example//reads10k.fq -t 8 --print-scaling=yes -K512 --cuda-block-size=64 --debug-break=yes > /dev/null

benchmark_cuda:
	./f5c -b test/chr22_meth_example/reads10k.bam -g test/chr22_meth_example//humangenome.fa -r test/chr22_meth_example//reads10k.fq -t 8 --print-scaling=yes -K512 --cuda-block-size=64 > /dev/null

benchmark_cpu:
	./f5c -b test/chr22_meth_example/reads10k.bam -g test/chr22_meth_example//humangenome.fa -r test/chr22_meth_example//reads10k.fq -t 8 --print-scaling=yes -K512 --disable-cuda=yes > /dev/null


run:
	./f5c -b test/chr22_meth_example/reads10k.bam -g test/chr22_meth_example//humangenome.fa -r test/chr22_meth_example//reads10k.fq -t 8 --print-scaling=yes -K512 --cuda-block-size=64 --debug-break=yes > /dev/null
