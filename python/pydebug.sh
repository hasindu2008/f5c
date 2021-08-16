#!/bin/bash

make clean
rm -rf *.so python/pyf5c.cpp build/lib.* build/temp.*
make -C slow5lib
GCC_ASAN_PRELOAD=$(gcc -print-file-name=libasan.so)
CC=g++ CFLAGS="-fsanitize=address -fno-omit-frame-pointer" python3 setup.py build
cp build/lib.*/*.so  ./
echo $GCC_ASAN_PRELOAD
LD_PRELOAD=$GCC_ASAN_PRELOAD  LD_LIBRARY_PATH=htslib/:slow5lib/lib/ python3 < python/example.py
