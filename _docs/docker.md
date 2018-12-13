---
title: Docker Image
---

To build a docker image
```sh
docker build .
```

Note down the image UUID and run f5c:
```sh
docker run -v /path/to/local/data/data/:/data/ -it :image_id ./f5c call-methylation -r /data/reads.fa -b /data/alignments.sorted.bam -g /data/ref.fa
```
