CC       = g++
CFLAGS   = -g -rdynamic -Wall -O2 -std=c++11 
CPPFLAGS =

include config.mk

#LIB_LDFLAGS = -lghts

# INC = $(HDF5_INC) 
#$(HTS_INC)
LDFLAGS += $(LIBS) -lpthread 

#SRC = $(wildcard *.c)
SRC = main.c f5c.c events.c nanopolish_read_db.c model.c align.c meth.c hmm.c
OBJ = $(SRC:.c=.o)
BINARY = f5c
DEPS = f5c.h fast5lite.h nanopolish_read_db.h f5cmisc.h error.h

ifeq ($(cuda),) #if cuda is undefined

else
	DEPS_CUDA = f5c.h fast5lite.h error.h f5cmisc.cuh
	SRC_CUDA = f5c.cu align.cu align_single.cu
	OBJ_CUDA = $(SRC_CUDA:.cu=_cuda.o)
	CC_CUDA = nvcc
	#CFLAGS_CUDA = -g  -G -Xcompiler -rdynamic  -O2 -std=c++11
	CFLAGS_CUDA = -g  -O2 -std=c++11 -lineinfo $(CUDA_ARCH) -maxrregcount=32
	#CFLAGS_CUDA = -g  -O2 -std=c++11 -lineinfo -arch=sm_61 -maxrregcount=32
	LDFLAGS += -L/usr/local/cuda/lib64/ -lcudart -lcudadevrt
	OBJ += gpucode.o $(OBJ_CUDA)
	CFLAGS += -DHAVE_CUDA=1
endif	


.PHONY: clean distclean format test

$(BINARY) : $(OBJ) 
	$(CC) $(CFLAGS) $(OBJ) $(LDFLAGS) -o $@


%.o : %.c $(DEPS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c 
	

gpucode.o : $(OBJ_CUDA)
	$(CC_CUDA) $(CFLAGS_CUDA) -dlink $^ -o $@

%_cuda.o : %.cu $(DEPS_CUDA)
	$(CC_CUDA) -x cu $(CFLAGS_CUDA) $(CPPFLAGS) -rdc=true -c $< -o $@

	
clean: 
	rm -rf f5c *.o *.out

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

benchmark:
	./f5c -b test/chr22_meth_example/reads10k.bam -g test/chr22_meth_example//humangenome.fa -r test/chr22_meth_example//reads10k.fq -t 8 --print-scaling=yes -K512 --cuda-block-size=64 > /dev/null


run:
	./f5c -b test/chr22_meth_example/reads10k.bam -g test/chr22_meth_example//humangenome.fa -r test/chr22_meth_example//reads10k.fq -t 8 --print-scaling=yes -K512 --cuda-block-size=64 --debug-break=yes > /dev/null

