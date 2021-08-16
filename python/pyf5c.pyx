# distutils: language = c++
# cython: language_level=3

from libc.stdlib cimport malloc, free
from libc.string cimport strdup
from pyf5c cimport *


def f5c_index_python(fast5_dir, fastq_file, num_threads, num_iop, verbosity_level=1):

    args = ["index", fastq_file, "-t", str(num_threads), "--iop", str(num_iop), '-v', str(verbosity_level) ]
    for dir in fast5_dir:
        args.append("-d")
        args.append(dir)


    num_args = len(args)
    c_argv = <char**> malloc(sizeof(char*) * num_args)

    for i in range(len(args)):
        string_a = str.encode(args[i])
        c_argv[i] = strdup(<char *> string_a)

    index_main(<int> num_args,  c_argv)

    for i in range(num_args):
        free(c_argv[i])
    free(c_argv)
    return
