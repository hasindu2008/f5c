#!/bin/bash

#rhel 8.8


VERSION=`git describe --tags`
make distclean
rm -r f5c-*
scripts/install-hdf5.sh && ./scripts/install-hts.sh && ./scripts/install-zstd.sh
autoreconf --install
./configure --enable-localhdf5 --enable-localzstd
make rocm=1 -j ROCM_ARCH='"--offload-arch=gfx1030 --offload-arch=gfx1100 --offload-arch=gfx900 --offload-arch=gfx906 --offload-arch=gfx908 --offload-arch=gfx90a --offload-arch=gfx942"'

mkdir lib/
patchelf --force-rpath  --set-rpath '$ORIGIN/lib' f5c
##https://github.com/BonsonW/slorado/releases/download/v0.3.0-beta/slorado-v0.3.0-beta-x86_64-rocm-linux-binaries.tar.xz
##5.7 rocm
cp ~/slorado/thirdparty/torch/libtorch/lib/libamd* lib/
cp ~/slorado/thirdparty/torch/libtorch/lib/libhsa-runtime64.so lib/
cp ~/slorado/thirdparty/torch/libtorch/lib/libnuma.so lib/
cp ~/slorado/thirdparty/torch/libtorch/lib/libtinfo.so lib/
cp ~/slorado/thirdparty/torch/libtorch/lib/libelf.so lib/
cp ~/slorado/thirdparty/torch/libtorch/lib/libdrm.so lib/
cp ~/slorado/thirdparty/torch/libtorch/lib/libdrm_amdgpu.so lib/

for f in lib/*.so*; do
   patchelf --force-rpath  --set-rpath '$ORIGIN' $f
done

mkdir f5c-${VERSION}
mv f5c f5c-${VERSION}/f5c_x86_64_linux_rocm
mv lib f5c-${VERSION}/
cp -r README.md LICENSE docs f5c-${VERSION}/
mkdir -p f5c-${VERSION}/scripts && cp scripts/test.sh scripts/common.sh scripts/test.awk scripts/install-vbz.sh f5c-${VERSION}/scripts
tar -zcf f5c-${VERSION}-binaries-rocm-specific.tar.gz f5c-${VERSION}