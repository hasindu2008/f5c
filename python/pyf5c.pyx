
# distutils: language = c++
# cython: language_level=3

from libc.stdlib cimport malloc, free
from libc.string cimport strdup
from pyf5c cimport *


def f5c_index_python(fast5_dir, fastq_file, num_threads, num_iop):

    args = ["index","-d", fast5_dir, fastq_file, "-t", num_threads, "--iop", num_iop]
    num_args=len(args)
    c_argv = <char**> malloc(sizeof(char*) * num_args)
    idx=0
    for i in range(num_args):
       string_a = str.encode(args[i])
       c_argv[i] = strdup(<char *> string_a)\

    index_main(<int> num_args,  c_argv)

    for i in range(num_args):
       free(c_argv[i])
    free(c_argv)
    return
