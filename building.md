#Building f5c

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
git clone https://github.com/hasindug/f5c
cd f5c
autoreconf #not required if a release
scripts/install-hts.sh
./configure
make
```

#### Method 2 (time consuming)

Dependencies : Install the HDF5 and zlib development libraries
```
On Debian/Ubuntu : sudo apt-get install zlib1g-dev
On Fedora/CentOS : sudo dnf/yum install zlib-devel
```

Now build f5c
```
git clone https://github.com/hasindug/f5c
cd f5c
autoreconf #optional
scripts/install-hts.sh #optional
scripts/install-hdf5.sh #optional
./configure --enable-localhdf5 #optional
make
```

#### Method 3 (not recommended)

Dependencies : Install HDF5 and hts
```
On Debian/Ubuntu : sudo apt-get install libhdf5-dev zlib1g-dev libhts1
```

Now build f5c
```
git clone https://github.com/hasindug/f5c
cd f5c
autoreconf
./configure --enable-systemhts
make
```
