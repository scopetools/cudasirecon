rm -rf cmake_build

if [ ! -d "src/IVE" ] 
then
    wget src/ https://www.dropbox.com/s/2twvw0go3dr3aim/IVE.zip
    unzip -o -d src IVE.zip
    rm IVE.zip
fi


mkdir cmake_build
cd cmake_build

# This is just to remove the -std= from CXXFLAGS ... probably a better way
INCLUDE=${CONDA_PREFIX}/include
CXXFLAGS="-Wfatal-errors -fvisibility-inlines-hidden -fmessage-length=0 \
          -march=nocona -mtune=haswell -ftree-vectorize \
          -fPIC -fstack-protector-strong -fno-plt -O2 \
          -ffunction-sections -pipe -isystem $INCLUDE"

cmake $CMAKE_ARGS \
    -DCMAKE_C_COMPILER=$CC \
    -DCMAKE_CXX_COMPILER=$CXX \
    -DCMAKE_BUILD_TYPE=Release \
    ../src

make -j 4
# make install