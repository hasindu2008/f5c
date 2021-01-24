import pyf5c

print("test launch")

fast5_dirs = ["test/ecoli_2kb_region/fast5_files"]
fastq_file = "test/ecoli_2kb_region/reads.fasta"
num_threads = 8
io_processes = 8

pyf5c.f5c_index_python(fast5_dirs, fastq_file, num_threads, io_processes)
print("done")
