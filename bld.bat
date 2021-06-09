del /Q cmake_build

@REM if not exist src\IVE\ (
@REM     powershell -Command "Invoke-WebRequest https://www.dropbox.com/s/2twvw0go3dr3aim/IVE.zip -OutFile IVE.zip"
@REM     powershell -Command "Expand-Archive IVE.zip -DestinationPath src\"
@REM     del IVE.zip )

mkdir cmake_build
cd cmake_build

cmake -G Ninja ^
    -DBUILD_MRC=ON ^
    -DBUILD_OTF_VIEWER=OFF ^
    -DCMAKE_BUILD_TYPE=Release ^
    -DCMAKE_INSTALL_PREFIX="%CONDA_PREFIX%/Library" ^
    ../src

ninja
ninja install
