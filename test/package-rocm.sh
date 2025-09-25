#!/bin/bash

VERSION=`git describe --tags`
make distclean

scripts/install-hdf5.sh && ./scripts/install-hts.sh && ./scripts/install-zstd.sh
./configure --enable-localhdf5 --enable-localzstd
make rocm=1 -j ROCM_ARCH='"--offload-arch=gfx1030 --offload-arch=gfx1100 --offload-arch=gfx900 --offload-arch=gfx906 --offload-arch=gfx908 --offload-arch=gfx90a"'

mkdir lib/
patchelf --set-rpath '$ORIGIN/lib' f5c
cp ~/slorado/thirdparty/torch/libtorch/lib/libamd* lib/
cp ~/slorado/thirdparty/torch/libtorch/lib/libhsa-runtime64.so lib/
cp ~/slorado/thirdparty/torch/libtorch/lib/libnuma.so lib/
cp ~/slorado/thirdparty/torch/libtorch/lib/libtinfo.so lib/
cp ~/slorado/thirdparty/torch/libtorch/lib/libelf.so lib/
cp ~/slorado/thirdparty/torch/libtorch/lib/libdrm.so lib/
cp ~/slorado/thirdparty/torch/libtorch/lib/libdrm_amdgpu.so lib/

mkdir f5c-${VERSION}
mv f5c f5c-${VERSION}/f5c_x86_64_linux_rocm
mv lib f5c-${VERSION}/
cp -r README.md LICENSE docs f5c-${VERSION}/
mkdir -p f5c-${VERSION}/scripts && cp scripts/test.sh scripts/common.sh scripts/test.awk scripts/install-vbz.sh f5c-${VERSION}/scripts
tar -zcf f5c-${VERSION}-binaries-rocm-specific.tar.gz f5c-${VERSION}