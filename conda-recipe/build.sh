#!/bin/bash

# NOTE: (when I wrote this note) The linux build pulled from conda packages for the 
# cudatoolkit, whereas the mac version assumes you have it installed in the local environment
# AND, for mac: the environmental variable CUDA_VERSION must be set, e.g.:
# export CUDA_VERSION=10.1 && conda build conda-recipe

if [ `uname` == Linux ]; then
    export LDFLAGS="-L${PREFIX}/lib"
    export CC=x86_64-conda_cos6-linux-gnu-gcc
    export CXX=x86_64-conda_cos6-linux-gnu-g++
fi 

if [ `uname` == Darwin ]; then
    export CC=gcc
    export CXX=g++
fi

rm -rf build
mkdir build
cd build

# build makeotf
wget http://www.fftw.org/fftw-2.1.5.tar.gz
tar -zxvf fftw-2.1.5.tar.gz
cd fftw-2.1.5
./configure --prefix=${PREFIX} --enable-type-prefix --enable-float --enable-threads
make
make install

$CC "${SRC_DIR}/otf/makeotf.c" \
    -I"${SRC_DIR}/IVE/darwin64/INCLUDE" \
    -I"${PREFIX}/include" \
    -L"${SRC_DIR}/IVE/darwin64/LIB" \
    -L"${PREFIX}/lib" \
    -limlib -lsrfftw -lsfftw \
    -o "${PREFIX}/bin/makeotf" 

cd "${SRC_DIR}/build"

if [ `uname` == Linux ]; then
    cmake .. \
        -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
        -DCMAKE_INSTALL_RPATH:STRING="${PREFIX}/lib"
#        -DCMAKE_NO_BUILTIN_CHRPATH:BOOL=ON
fi

if [ `uname` == Darwin ]; then
    install_name_tool -add_rpath @executable_path/../lib "${PREFIX}/bin/makeotf"

    echo "CUDA_VERSION: ${CUDA_VERSION}"
    CUDA_TOOLKIT_ROOT_DIR="/Developer/NVIDIA/CUDA-${CUDA_VERSION}"
    echo "CUDA_TOOLKIT_ROOT_DIR: ${CUDA_TOOLKIT_ROOT_DIR}"
    cmake .. \
        -Wno-dev \
        -D CUDA_TOOLKIT_ROOT_DIR="${CUDA_TOOLKIT_ROOT_DIR}" \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
        -DCMAKE_INSTALL_RPATH:STRING="${PREFIX}/lib"
fi

make -j4
make install
ln -s ${PREFIX}/bin/cudaSireconDriver ${PREFIX}/bin/cudasirecon
ln -s ${PREFIX}/bin/cudaSireconDriver ${PREFIX}/bin/sirecon

if [ `uname` == Darwin ]; then
    cp "${CUDA_TOOLKIT_ROOT_DIR}"/lib/libcufft.*.dylib "${PREFIX}"/lib/
fi