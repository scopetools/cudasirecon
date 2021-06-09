rm -rf cmake_build

if [ ! -d "src/IVE" ] 
then
    wget src/ https://www.dropbox.com/s/2twvw0go3dr3aim/IVE.zip
    unzip -o -d src IVE.zip
    rm IVE.zip
fi


mkdir cmake_build
cd cmake_build

# This is just to remove the -std= from CXXFLAGS added by conda
# ... probably a better way
CXXFLAGS="-Wfatal-errors -Wno-deprecated-declarations -fPIC -pipe -isystem ${CONDA_PREFIX}/include"

cmake ../src $CMAKE_ARGS \
    -DCMAKE_C_COMPILER=$CC \
    -DCMAKE_CXX_COMPILER=$CXX \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_MRC=OFF
    # -DBUILD_OTF_VIEWER=ON \
    

make -j 4
# make install