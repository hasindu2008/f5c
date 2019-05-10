#!/usr/bin/env bash

autoreconf
./scripts/install-hts.sh
./scripts/install-hdf5.sh
./configure --enable-localhdf5
