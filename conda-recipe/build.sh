#!/bin/bash
echo "STARTING build.sh"
echo $(pwd)

rm -rf build
mkdir build
cd build

if [ `uname` == Linux ]; then
    export LDFLAGS="-L${PREFIX}/lib"
    export CC=x86_64-conda_cos6-linux-gnu-gcc
    export CXX=x86_64-conda_cos6-linux-gnu-g++
    cmake .. \
        -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
        -DCMAKE_INSTALL_RPATH:STRING="${PREFIX}/lib"
#        -DCMAKE_NO_BUILTIN_CHRPATH:BOOL=ON
fi

if [ `uname` == Darwin ]; then
    export CC=gcc
    export CXX=g++
    cmake .. \
        -Wno-dev \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
        -DCMAKE_INSTALL_RPATH:STRING="${PREFIX}/lib"
fi


make
echo "MAKE DONE"
make install
echo "INSTALL DONE"
ln -s ${PREFIX}/bin/cudaSireconDriver ${PREFIX}/bin/cudasirecon
