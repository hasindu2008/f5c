AUTOMAKE_OPTIONS = foreign

CC     = g++
CFLAGS   = -g -rdynamic -Wall -O2 -std=c++11 
CPPGLAGS =


HTS_INC =
HDF5_INC = /usr/include/hdf5/serial/ 
HTS_LDFLAGS = -lhts
HDF5_LDFLAGS = -lhdf5_serial

INC = $(HDF5_INC) $(HTS_INC)
LDFLAGS = $(HTS_LDFLAGS) $(HDF5_LDFLAGS) -lz -lm -lbz2 -llzma -lpthread 

#SRC = $(wildcard *.c)
SRC = main.c f5c.c events.c nanopolish_read_db.c
OBJ = $(SRC:.c=.o)
BINARY = f5c
DEPS = f5c.h fast5lite.h nanopolish_read_db.h f5cmisc.h

.PHONY: clean  

$(BINARY) : $(OBJ) 
	$(CC) $(CFLAGS) $(OBJ) $(LDFLAGS)-o $@


%.o : %.c $(DEPS)
	$(CC) $(CFLAGS) -I $(INC) $< -c 
	
	
clean: 
	rm -rf f5c *.o

