#!/bin/bash

#Todo : use curl as a failsafe

# terminate script
die() {
    echo -e "${RED}$1${NC}" >&2
    echo
    exit 1
}

ZSTD_VERSION=1.3.1
test -e zstd-v${ZSTD_VERSION}.tar.gz && rm zstd-v${ZSTD_VERSION}.tar.gz
test -d zstd && rm -r zstd

wget https://github.com/facebook/zstd/archive/refs/tags/v${ZSTD_VERSION}.tar.gz -O zstd-v${ZSTD_VERSION}.tar.gz || die "Downloading zstd source failed"
tar -xzf zstd-v${ZSTD_VERSION}.tar.gz || die "extracting failed"
rm zstd-v${ZSTD_VERSION}.tar.gz
mv zstd-${ZSTD_VERSION} zstd || die "Moving zstd-v${ZSTD_VERSION} to zstd failed"
cd zstd || die "cd to zstd failed"
make -j8 || die " Building zstd failed"
echo "Successfully installed zstd to ./zstd."
