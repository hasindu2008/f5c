test -e  https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && rm htslib-1.9.tar.bz2
test -d htslib && rm -r htslib
wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 || exit 1
tar -xf htslib-1.9.tar.bz2  || exit
rm htslib-1.9.tar.bz2
mv htslib-1.9 htslib || exit
cd htslib || exit 1
chmod +x ./configure || exit 1
./configure --disable-bz2 --disable-lzma || exit
make -j8 || exit
echo "Successfully installed htslib to ./htslib."
echo 'Now run ./configure again!'
