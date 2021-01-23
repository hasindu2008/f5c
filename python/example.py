import pyf5c

print("test launch")
pyf5c.f5c_index_python("test/ecoli_2kb_region/fast5_files","test/ecoli_2kb_region/reads.fasta", "8", "8")
print("done")
