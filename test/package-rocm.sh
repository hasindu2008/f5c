#!/bin/bash

#rhel 8.8

die() {
    echo "$*" 1>&2
    exit 1
}

VERSION=`git describe --tags`
make distclean
rm -r f5c-*
scripts/install-hdf5.sh && ./scripts/install-hts.sh && ./scripts/install-zstd.sh || die "Failed to install dependencies"
autoreconf --install || die "Failed to autoreconf"
./configure --enable-localhdf5 --enable-localzstd || die "Failed to configure"

#--offload-arch=gfx950
MI900="--offload-arch=gfx803 --offload-arch=gfx900 --offload-arch=gfx906 --offload-arch=gfx908 --offload-arch=gfx90a --offload-arch=gfx942"
RAD10="--offload-arch=gfx1010 --offload-arch=gfx1011 --offload-arch=gfx1012 --offload-arch=gfx1030 --offload-arch=gfx1031 --offload-arch=gfx1032"
RAD11="--offload-arch=gfx1100 --offload-arch=gfx1101 --offload-arch=gfx1102"
#RAD12="--offload-arch=gfx1200 --offload-arch=gfx1201"
# --offload-arch=gfx1150
IGPU="--offload-arch=gfx1035 --offload-arch=gfx1103"
make rocm=1 -j ROCM_ARCH="$MI900 $RAD10 $RAD11 $RAD12 $IGPU" || die "Failed to make"
#make rocm=1 -j ROCM_ARCH='"--offload-arch=gfx1030 --offload-arch=gfx1100 --offload-arch=gfx900 --offload-arch=gfx906 --offload-arch=gfx908 --offload-arch=gfx90a --offload-arch=gfx942"'

mkdir lib/ || die "Failed to make lib dir"
patchelf --force-rpath  --set-rpath '$ORIGIN/lib' f5c || die "Failed to patchelf f5c"
##https://github.com/BonsonW/slorado/releases/download/v0.3.0-beta/slorado-v0.3.0-beta-x86_64-rocm-linux-binaries.tar.xz
##5.7 rocm
cp ~/slorado/thirdparty/torch/libtorch/lib/libamd* lib/ || die "Failed to copy amd libs"
cp ~/slorado/thirdparty/torch/libtorch/lib/libhsa-runtime64.so lib/ || die "Failed to copy hsa-runtime64.so"
cp ~/slorado/thirdparty/torch/libtorch/lib/libnuma.so lib/ || die "Failed to copy libnuma.so"
cp ~/slorado/thirdparty/torch/libtorch/lib/libtinfo.so lib/ || die "Failed to copy libtinfo.so"
cp ~/slorado/thirdparty/torch/libtorch/lib/libelf.so lib/ || die "Failed to copy libelf.so"
cp ~/slorado/thirdparty/torch/libtorch/lib/libdrm.so lib/ || die "Failed to copy libdrm.so"
cp ~/slorado/thirdparty/torch/libtorch/lib/libdrm_amdgpu.so lib/ || die "Failed to copy libdrm_amdgpu.so"

for f in lib/*.so*; do
   patchelf --force-rpath  --set-rpath '$ORIGIN' $f || die "Failed to patchelf $f"
done

mkdir f5c-${VERSION} || die "Failed to make f5c-${VERSION} dir"
mv f5c f5c-${VERSION}/f5c_x86_64_linux_rocm || die "Failed to move f5c"
mv lib f5c-${VERSION}/ || die "Failed to move lib"
cp -r README.md LICENSE docs f5c-${VERSION}/ || die "Failed to copy docs"
mkdir -p f5c-${VERSION}/scripts && cp scripts/test.sh scripts/common.sh scripts/test.awk scripts/install-vbz.sh f5c-${VERSION}/scripts || die "Failed to copy scripts"
tar -zcf f5c-${VERSION}-rocm-binaries-experimental.tar.gz f5c-${VERSION} || die "Failed to make tar"