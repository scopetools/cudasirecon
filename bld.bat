mkdir cmake_build
cd cmake_build

cmake -G Ninja ^
    -DBUILD_MRC=OFF ^
    -DCMAKE_BUILD_TYPE=Release ^
    ../src

ninja
@REM ninja install
