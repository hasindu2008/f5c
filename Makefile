CC     = g++
CFLAGS   = -g -Wall -O2 -std=c++11
LDFLAGS = -lhts -lhdf5_serial -lz -lm -lbz2 -llzma -lpthread 
INC = 

DEPS = common.h fast5lite.h nanopolish_read_db.h main.h
OBJ = main.o common.o nanopolish_read_db.o process.o
BINARY = main

$(BINARY) : $(OBJ) 
	$(CC) $(CFLAGS) $(OBJ) $(LDFLAGS) -o $@


main.o: main.c nanopolish_read_db.o  $(DEPS) 
	$(CC) $(CFLAGS) $< $(LDFLAGS) -c 

common.o : common.c $(DEPS) 
	$(CC) $(CFLAGS) $< $(LDFLAGS) -c 

nanopolish_read_db.o : nanopolish_read_db.cpp $(DEPS) 
	$(CC) $(CFLAGS) $< $(LDFLAGS) -c

process.o : process.c $(DEPS)
	$(CC) $(CFLAGS) $< $(LDFLAGS) -c

.PHONY: clean  
clean: 
	-rm main
	-rm *.o
