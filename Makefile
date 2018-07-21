CC     = g++
CFLAGS   = -g -Wall -O2 -std=c++11
CPPGLAGS =
INC = 
LDFLAGS = -lhts -lhdf5_serial -lz -lm -lbz2 -llzma -lpthread 

#SRC = $(wildcard *.c)
SRC = main.c f5c.c events.c nanopolish_read_db.c
OBJ = $(SRC:.c=.o)
BINARY = f5c
DEPS = f5c.h fast5lite.h nanopolish_read_db.h f5cmisc.h

.PHONY: clean  

$(BINARY) : $(OBJ) 
	$(CC) $(CFLAGS) $(OBJ) $(LDFLAGS) -o $@


%.o : %.c $(DEPS)
	$(CC) $(CFLAGS) $< -c 
	

clean: 
	rm -rf f5c *.o

