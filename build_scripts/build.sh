if [ -d build ]
then
rm -r build
fi
mkdir build && cd build
cmake ..
make
ctest
cpack -G DEB
mv addcomp-0.1.1-Linux.deb addcomp.deb
cd ..

