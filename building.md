# Building f5c

Note : Building from the Github repository requires `autoreconf` which can be installed on Ubuntu using `sudo apt-get install autoconf automake`.

Clone the git repository.
```
git clone https://github.com/hasindug/f5c && cd f5c
```
Alternatively, download the [latest release](https://github.com/hasindu2008/f5c/releases) tarball and exract. 
eg :
```
wget "https://github.com/hasindu2008/f5c/releases/download/v0.0-alpha/f5c-v0.0-alpha-release.tar.gz" && tar xvf f5c-v0.0-alpha-release.tar.gz && cd f5c-v0.0-alpha/
```

While we have tried hard to avoid the dependency hell, three dependencies (zlib, HDF5 and HTS) could not be avoided.

Currently 3 building methods are supported.
1. Locally compiled HTS library and system wide HDF5 library (recommended)
2. Locally compiled HTS and HDF5 libraries (HDF5 local compilation - takes a bit of time)
3. System wide HTS and HDF5 libraries (not recommended as HTS versions can be old)


#### Method 1 (recommended)

Dependencies : Install the HDF5 (and zlib development libraries)
```
On Debian/Ubuntu : sudo apt-get install libhdf5-dev zlib1g-dev
On Fedora/CentOS : sudo dnf/yum install hdf5-devel zlib-devel
On Arch Linux: sudo pacman -S hdf5
On OS X : brew install hdf5
```

Now build f5c
```
autoreconf              # skip if compiling a release, only required when building from github
scripts/install-hts.sh  # download and compiles htslib in the current folder
./configure
make                    # or make cuda=1 if compiling for CUDA
```

#### Method 2 (time consuming)

Dependencies : Install the zlib development libraries
```
On Debian/Ubuntu : sudo apt-get install zlib1g-dev
On Fedora/CentOS : sudo dnf/yum install zlib-devel
```

Now build f5c
```
autoreconf                      # skip if compiling a release, only required when building from github
scripts/install-hts.sh          # download and compiles htslib in the current folder
scripts/install-hdf5.sh         # download and compiles HDF5 in the current folder
./configure --enable-localhdf5
make                            # or make cuda=1 if compiling for CUDA
```

#### Method 3 (not recommended)

Dependencies : Install HDF5 and hts
```
On Debian/Ubuntu : sudo apt-get install libhdf5-dev zlib1g-dev libhts1
```

Now build f5c
```
autoreconf                      # skip if compiling a release, only required when building from github
./configure --enable-systemhts
make                            # or make cuda=1 if compiling for CUDA
```
