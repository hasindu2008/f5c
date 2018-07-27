CC     = g++
CFLAGS   = -g -rdynamic -Wall -O2 -std=c++11 
CPPGLAGS =

include config.mk

#LIB_LDFLAGS = -lghts

# INC = $(HDF5_INC) 
#$(HTS_INC)
LDFLAGS += $(LIBS) 

#SRC = $(wildcard *.c)
SRC = main.c f5c.c events.c nanopolish_read_db.c
OBJ = $(SRC:.c=.o)
BINARY = f5c
DEPS = f5c.h fast5lite.h nanopolish_read_db.h f5cmisc.h

.PHONY: clean  

$(BINARY) : $(OBJ) 
	$(CC) $(CFLAGS) $(OBJ) $(LDFLAGS)-o $@


%.o : %.c $(DEPS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c 
	
	
clean: 
	rm -rf f5c *.o *.out
