---
title: Docker Image
---

To build a docker image
```
git clone https://github.com/hasindu2008/f5c && cd f5c
docker build .
```
Note down the image uuid and run f5c as :
```
docker run -v /path/to/local/data/data/:/data/ -it :image_id  ./f5c call-methylation -r /data/reads.fa -b /data/alignments.sorted.bam -g /data/ref.fa
```
