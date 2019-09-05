#!/bin/bash

# NOTE: (when I wrote this note) The linux build pulled from conda packages for the 
# cudatoolkit, whereas the mac version assumes you have it installed in the local environment
# AND, for mac: the environmental variable CUDA_VERSION must be set, e.g.:
# export CUDA_VERSION=10.1 && conda build conda-recipe
rm -rf build
mkdir build
cd build

if [ `uname` == Linux ]; then
    export LDFLAGS="-L${PREFIX}/lib"
    export CC=$PREFIX/bin/x86_64-conda_cos6-linux-gnu-gcc
    export CXX=$PREFIX/bin/x86_64-conda_cos6-linux-gnu-g++
    # export CC=gcc
    # export CXX=g++
    PLATFORM=linux64

    #CUDA_TOOLKIT_ROOT_DIR="/usr/local/cuda-${CUDA_VERSION}"
    #CUDA_LIB_DIR="${CUDA_TOOLKIT_ROOT_DIR}"/lib64

    cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
    -DCMAKE_INSTALL_RPATH:STRING="${PREFIX}/lib"
#    -DCUDA_TOOLKIT_ROOT_DIR="${CUDA_TOOLKIT_ROOT_DIR}" \

fi 

if [ `uname` == Darwin ]; then
    export CC=gcc
    export CXX=g++
    PLATFORM=darwin64
    CUDA_TOOLKIT_ROOT_DIR="/Developer/NVIDIA/CUDA-${CUDA_VERSION}"
    CUDA_LIB_DIR="${CUDA_TOOLKIT_ROOT_DIR}"/lib

    cmake .. \
        -DCUDA_TOOLKIT_ROOT_DIR="${CUDA_TOOLKIT_ROOT_DIR}" \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
        -DCMAKE_INSTALL_RPATH:STRING="${PREFIX}/lib"

fi


make -j4
make install
ln -s ${PREFIX}/bin/cudaSireconDriver ${PREFIX}/bin/cudasirecon
ln -s ${PREFIX}/bin/cudaSireconDriver ${PREFIX}/bin/sirecon

#if [ `uname` == Linux ]; then
#    cp "${CUDA_LIB_DIR}"/libcufft.*.so "${PREFIX}"/lib/
#fi

cd "${SRC_DIR}/build"

# install fftw
# wget http://www.fftw.org/fftw-2.1.5.tar.gz
# tar -zxvf fftw-2.1.5.tar.gz
# cd fftw-2.1.5
# ./configure --prefix=${PREFIX} --enable-type-prefix --enable-float --enable-threads
# make -j 4
# make install

# FFTW_ROOT="${SRC_DIR}/fftw2/${PLATFORM}/"

# # build makeotf
# $CC "${SRC_DIR}/otf/makeotf.c" \
#     -I"${SRC_DIR}/IVE/${PLATFORM}/INCLUDE" \
#     -I"${FFTW_ROOT}/include" \
#     -L"${SRC_DIR}/IVE/${PLATFORM}/LIB" \
#     -L"${FFTW_ROOT}/lib" \
#     -limlib -lsrfftw -lsfftw -lm \
#     -o "${PREFIX}/bin/makeotf" 

# if [ `uname` == Darwin ]; then
#     install_name_tool -add_rpath @executable_path/../lib "${PREFIX}/bin/makeotf"
# fi


if [ `uname` == Darwin ]; then
    cp "${CUDA_LIB_DIR}"/libcufft.*.dylib "${PREFIX}"/lib/
fi

echo "####################### build.sh done"