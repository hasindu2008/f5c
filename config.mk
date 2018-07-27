

LIBS = -lhdf5_serial -lhts -lz 
LDFLAGS = 
CPPFLAGS = 

ifeq "locallibhts-no" "locallibhts-yes"

CPPFLAGS += -I./htslib
LDFLAGS += htslib/libhts.a

endif


ifeq "locallibhdf5-no" "locallibhdf5-yes"

CPPFLAGS += -I./hdf5/include/
LDFLAGS += hdf5/lib/libhdf5.a -ldl

endif
