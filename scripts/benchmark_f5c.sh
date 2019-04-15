#nvprof  -f --kernels "align_kernel_core" --analysis-metrics -o bv.nvprof ./f5c -b test/chr22_meth_example/reads10k.bam -g test/chr22_meth_example//humangenome.fa -r test/chr22_meth_example//reads10k.fq -t 8 --print-scaling=yes -K512 --cuda-block-size=64 --debug-break=yes > /dev/null

nvprof  -f --analysis-metrics -o bv.nvprof ./f5c -b test/chr22_meth_example/reads10k.bam -g test/chr22_meth_example//humangenome.fa -r test/chr22_meth_example//reads10k.fq -t 8 --print-scaling=yes -K512 --cuda-block-size=64 --debug-break=yes > /dev/null

./f5c -b test/chr22_meth_example/reads10k.bam -g test/chr22_meth_example//humangenome.fa -r test/chr22_meth_example//reads10k.fq -t 8 --print-scaling=yes -K512 --cuda-block-size=64 > /dev/null

./f5c -b test/chr22_meth_example/reads10k.bam -g test/chr22_meth_example//humangenome.fa -r test/chr22_meth_example//reads10k.fq -t 8 --print-scaling=yes -K512 --disable-cuda=yes > /dev/null


./f5c -b test/chr22_meth_example/reads10k.bam -g test/chr22_meth_example//humangenome.fa -r test/chr22_meth_example//reads10k.fq -t 8 --print-scaling=yes -K512 --cuda-block-size=64 --debug-break=yes > /dev/null
