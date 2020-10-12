---
title: Miscellaneous Installation Methods
---

Note that at the moment, only the CPU version of *f5c* can be installed through *Docker* or *Conda*. If you have experience with *nvidia docker* or setting up CUDA enabled programmes on Conda, help us to set *f5c* GPU version up on those.

#### Installing Docker Image

To build a *docker* image
```
git clone https://github.com/hasindu2008/f5c && cd f5c
docker build .
```
Note down the image uuid and run f5c as :
```
docker run -v /path/to/local/data/data/:/data/ -it :image_id  ./f5c call-methylation -r /data/reads.fa -b /data/alignments.sorted.bam -g /data/ref.fa
```

#### Installing through Conda

To install through *Conda*
```
conda create -n ont_env -c bioconda -c conda-forge f5c
```

