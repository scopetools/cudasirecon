#!/bin/bash
echo "STARTING build.sh"
echo $(pwd)

rm -rf build
mkdir build
cd build

if [ `uname` == Linux ]; then
    export LDFLAGS="-L${PREFIX}/lib"
    # export CXXFLAGS="-L${PREFIX}/lib"
    export CC=x86_64-conda_cos6-linux-gnu-gcc
    export CXX=x86_64-conda_cos6-linux-gnu-g++
    cmake .. \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
        -DCMAKE_INSTALL_RPATH:STRING="${PREFIX}/lib"
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
make install
ln -s ${PREFIX}/bin/cudaSireconDriver ${PREFIX}/bin/cudasirecon