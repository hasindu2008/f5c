#!/usr/bin/env bash

autoreconf
./scripts/install-hdf5.sh
./scripts/install-hts.sh
./configure --enable-localhdf5
