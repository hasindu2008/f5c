from libc.stdint cimport *

cdef extern from "pyf5c.h":

    void index_main(int argc, char** argv)
    int meth_main(int argc, char* argv[], int8_t mode);
