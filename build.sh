mkdir cmake_build
cd cmake_build

# This is just to remove the -std= from CXXFLAGS ... probably a better way
INCLUDE=${CONDA_PREFIX}/include
CXXFLAGS="-fvisibility-inlines-hidden -fmessage-length=0 \
          -march=nocona -mtune=haswell -ftree-vectorize \
          -fPIC -fstack-protector-strong -fno-plt -O2 \
          -ffunction-sections -pipe -isystem $INCLUDE"

cmake $CMAKE_ARGS \
    -DCMAKE_C_COMPILER=$CC \
    -DCMAKE_CXX_COMPILER=$CXX \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=/home/tjl10/dev/cudasirecon/ \
    ../src

make -j 2
make install