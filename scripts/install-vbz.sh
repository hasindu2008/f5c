#!/bin/bash

# Download libvbz_hdf_plugin.so to $HOME/.local/lib
# Check https://github.com/nanoporetech/vbz_compression/releases
# User may have to change the $ONT_VBZ_VERSION and $ONT_VBZ_FILE.
# User may have to change the $PLUGIN_PATH
# Make sure to export the $PLUGIN_PATH using export HDF5_PLUGIN_PATH=$PLUGIN_PATH before using it

NC='\033[0m' # No Color
RED='\033[0;31m'
GREEN='\033[0;32m'

# terminate script
die() {
    echo -e "${RED}$1${NC}" >&2
    echo
    exit 1
}

print() {
    echo -e "${GREEN}$1${NC}" >&2
}

MANUAL_LINK="https://f5c.page.link/troubleshoot"

uname -o || die "Could not determine the O/S. See ${MANUAL_LINK}"
uname -m || die "Could not determine the architecture. See ${MANUAL_LINK}"

ARCH=$(uname -m)
OS=$(uname -o)


if [ "${OS}" != "GNU/Linux"  ];
then
    die "Unhandled O/S ${OS}. See  ${MANUAL_LINK}"
fi

if [[ ${ARCH} != "x86_64"  && ${ARCH} != "aarch64" ]];
then
    die "Unhandled architecture ${ARCH}. See ${MANUAL_LINK}"
fi

ONT_VBZ_VERSION=v1.0.1
ONT_VBZ_FILE=ont-vbz-hdf-plugin-1.0.1-Linux-${ARCH}.tar.gz
WGET_LINK=https://github.com/nanoporetech/vbz_compression/releases/download/$ONT_VBZ_VERSION/$ONT_VBZ_FILE

PLUGIN_PATH=$HOME/.local/hdf5/lib/plugin
test -d $PLUGIN_PATH || mkdir -p $PLUGIN_PATH || die "Creating directory $PLUGIN_PATH failed"

ONT_VBZ_DIR=ont_vbz_plugin_$ONT_VBZ_VERSION
test -d $ONT_VBZ_DIR && rm -r "$ONT_VBZ_DIR"
mkdir $ONT_VBZ_DIR || die "Failed creating $ONT_VBZ_DIR"

wget $WGET_LINK || die "Could not download $ONT_VBZ_FILE"
tar -xzvf $ONT_VBZ_FILE -C $ONT_VBZ_DIR || die "Extracting $ONT_VBZ_FILE failed"
rm $ONT_VBZ_FILE || die "Cannot delete $ONT_VBZ_FILE"

find $ONT_VBZ_DIR -name '*.so' -exec mv -t $PLUGIN_PATH {} + || die "Could not move .so file to $PLUGIN_PATH. Check if .so exists inside $ONT_VBZ_DIR"
rm -r $ONT_VBZ_DIR || die "Cannot delete $ONT_VBZ_DIR"

export HDF5_PLUGIN_PATH=${PLUGIN_PATH} || die "exporting PLUGIN_PATH failed"

print "successfully installed the vbz plugin under ${PLUGIN_PATH}"
print "Please invoke the following on your shell"
print "export HDF5_PLUGIN_PATH=${PLUGIN_PATH} "
print "For persistence, add the above line to your ~/.bashrc"

exit 0
