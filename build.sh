rm -rf cmake_build

if [ ! -d "src/IVE" ] 
then
    wget src/ https://www.dropbox.com/s/2twvw0go3dr3aim/IVE.zip
    unzip -o -d src IVE.zip
    rm IVE.zip
fi


mkdir cmake_build
cd cmake_build

CXXFLAGS="$CXXFLAGS -Wfatal-errors -Wno-deprecated-declarations"
cmake ${CMAKE_ARGS} -DCMAKE_BUILD_TYPE=Release ../src

make -j 4
# make install