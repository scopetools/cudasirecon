mkdir cmake_build
cd cmake_build

cmake -G Ninja ^
    @REM -DCMAKE_PREFIX_PATH:PATH="%CONDA_PREFIX%" ^
    -DCMAKE_BUILD_TYPE=Release ^
    ../src
if errorlevel 1 exit 1

ninja
if errorlevel 1 exit 1

ninja install
if errorlevel 1 exit 1