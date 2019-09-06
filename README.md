# CUDA_SIMrecon

Mats Gustafsson & Lin Shao's 3-beam SIM reconstruction software, with CUDA acceleration.

Algorithm as described in [Gustafsson et al (2008) *Biophys*. **94(12)**: 4957–4970. doi: 10.1529/biophysj.107.120345](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2397368/)


## Installation

Packages for Linux and OS X are available through conda:

```bash
$ conda install -c talley -c conda-forge cudasirecon

# you may wish to put it in a new environment ("sim" here):
$ conda create -n sim -y -c talley -c conda-forge cudasirecon
$ conda activate sim
```

## Usage

```bash
# the binary will be available as
$ cudasirecon
# (or a shorter symlink on mac/linux for convenience)
$ sirecon

# for command line help
$ cudasirecon --help

# a typical call
$ cudasirecon data.dv data-PROC.dv otf.otf -c config

# ... where otf.otf was generated using the "makeotf" command:
# note, makeotf not currently working/shipping on windows
$ makeotf /path/to/psf.dv otf.otf -angle -1.855500 -ls 0.2075 -na 1.4 -nimm 1.515 -fixorigin 3 20 -leavekz 7 11 3

# for more help on the makeotf command:
$ makeotf -help
```

### Full list of options/flags

```
$ cudasirecon --help

  --input-file arg              input file (or data folder in TIFF mode)
  --output-file arg             output file (or filename pattern in TIFF mode)
  --otf-file arg                OTF file
  --usecorr arg                 use the flat-field correction file provided
  --ndirs arg (=3)              number of directions
  --nphases arg (=5)            number of phases per direction
  --nordersout arg (=0)         number of output orders; must be <= norders
  --angle0 arg (=1.648)         angle of the first direction in radians
  --ls arg (=0.172000006)       line spacing of SIM pattern in microns
  --na arg (=1.20000005)        Detection numerical aperture
  --nimm arg (=1.33000004)      refractive index of immersion medium
  --zoomfact arg (=2)           lateral zoom factor
  --explodefact arg (=1)        artificially exploding the reciprocal-space 
                                distance between orders by this factor
  --zzoom arg (=1)              axial zoom factor
  --nofilteroverlaps [=arg(=0)] do not filter the overlaping region between 
                                bands usually used in trouble shooting
  --background arg (=0)         camera readout background
  --wiener arg (=0.00999999978) Wiener constant
  --forcemodamp arg             modamps forced to these values
  --k0angles arg                user given pattern vector k0 angles for all 
                                directions
  --otfRA [=arg(=1)]            using rotationally averaged OTF
  --otfPerAngle [=arg(=1)]      using one OTF per SIM angle
  --fastSI [=arg(=1)]           SIM data is organized in Z->Angle->Phase order;
                                default being Angle->Z->Phase
  --k0searchAll [=arg(=0)]      search for k0 at all time points
  --equalizez [=arg(=1)]        bleach correcting for z
  --equalizet [=arg(=1)]        bleach correcting for time
  --dampenOrder0 [=arg(=1)]     dampen order-0 in final assembly
  --nosuppress [=arg(=0)]       do not suppress DC singularity in final 
                                assembly (good idea for 2D/TIRF data)
  --nokz0 [=arg(=1)]            do not use kz=0 plane of the 0th order in the 
                                final assembly
  --gammaApo arg (=1)           output apodization gamma; 1.0 means triangular 
                                apo
  --saveprefiltered arg         save separated bands (half Fourier space) into 
                                a file and exit
  --savealignedraw arg          save drift-fixed raw data (half Fourier space) 
                                into a file and exit
  --saveoverlaps arg            save overlap0 and overlap1 (real-space complex 
                                data) into a file and exit
  -c [ --config ] arg           name of a file of a configuration.
  --2lenses [=arg(=1)]          I5S data
  --writeTitle [=arg(=1)]       Write command line to image header (may cause 
                                issues with bioformats)
  -h [ --help ]                 produce help message
```

### Config file

The config file can specify any flags/options listed above, and a typical 3D sim config file may look like this:

```bash
nimm=1.515
background=90
wiener=0.001
# angles of illumination in radians
k0angles=-0.804300,-1.8555,0.238800
ls=0.2035
ndirs=3
nphases=5
na=1.42
otfRA=1
dampenOrder0=1
```

## Requirements

* Currently only accepts images as .dv or .mrc format.  If you need to convert TIFF files (or any other format you can get into a numpy array) to DV/MRC format you can install the [mrc python package](https://github.com/tlambert03/mrc) with `pip install mrc`.  Then use something like:
```python
import tifffile
import mrc

fname = '/path/to/file.tif'
im = tifffile.imread(fname)
# it will be critical that the pixel calibrations 'dx','dy', 'dz' are correct for sim reconstruction
mrc.imwrite(fname.replace('.tif', '.dv'), im,
            metadata={'dx': 0.08, 'dy': 0.08, 'dz': 0.25, 'wave': [525,0,0,0,0]})
```
  
* Requires CUDA-capable NVIDIA GPU and driver.

The program has been compiled against different versions of the CUDA toolkit. The required CUDA libraries are bundled in the conda distributions so you don't need to install the CUDA toolkit separately. If desired, you can pick which version of CUDA you'd like based on your needs, but please note that different versions of the CUDA toolkit have different GPU driver requirements.  Not all versions are available on all platforms.  To see what versions are available on your platform, type `conda search -c talley cudasirecon`.


| build  | min CUDA driver | Install With |
| ------------- | ------------  | -----------  |
| 10.1  | ≥ 418.39     | `conda install cudasirecon=*=cu10.1`  |
| 10.0  | ≥ 410.48     | `conda install cudasirecon=*=cu10.0`  |
|  9.2  | ≥ 396.26    | `conda install cudasirecon=*=cu9.2`  |
|  9.0  | ≥ 384.81    | `conda install cudasirecon=*=cu9.0`  |

If your CUDA Driver version is too low for the version of cudasirecon that you have installed, you may get an error that looks like: `!!Error occurred: cudaSetDevice failed`

If you run into trouble, feel free to [open an issue](https://github.com/tlambert03/CUDA_SIMrecon/issues) and describe your setup.

# Compiling from source

Building the binary from source can be somewhat tricky (hence the conda packages), but if you'd like to build from scratch, here are some notes for each platform.  You can also look in the `conda-recipe` folder for tips.

## All platforms

The program requires the IVE/Priism libraries, which are not distributed with this source code and must be acquired seperately from UCSF.  Place them in a folder called `IVE` at the top level of the source folder (same folder as the cudaSirecon folder).  It should minimally have the following files and folders (example shown for linux, use `.a` or `.lib` as necessary for osx or windows)

```
CUDA_SIMrecon
├── ...
└── IVE/
    ├── darwin64/
    │   ├── INCLUDE/
    │   └── LIB/
    ├── linux64/
    │   ├── INCLUDE/
    │   │   ├── IM.h
    │   │   ├── IMInclude.h
    │   │   └── IWApiConstants.h
    │   └── LIB/
    │       ├── libimlib.a
    │       └── libive.a
    └── win64/
        ├── INCLUDE/
        └── LIB/
```

## Building on Linux
This has only been tested on Ubuntu 16.04.  You will need to have the [NVIDIA CUDA toolkit](https://developer.nvidia.com/cuda-toolkit) installed (I have used versions 8.0 - 10.1 ... so that part doesn't really matter unless you are targetting a GPU with a specific CUDA compute capability)... See note below about optionally installing with conda.

I use [conda](https://docs.conda.io/en/latest/miniconda.html) for the remaining dependencies and build as follows:

```bash
$ conda create -n simbuild -c conda-forge -y gcc_linux-64=5.4.0 gxx_linux-64=5.4.0 cmake liblapack boost-cpp xorg-libx11
$ conda activate simbuild

# optional: if you want to install the CUDA toolkit through conda rather than the NVIDIA website,
# you need to use the dev versions that have the nvcc compiler.
# conda install -c conda-forge conda cudatookit-dev=10.0

# create a build directory inside of CUDA_SIMrecon
$ mkdir build
$ cd build
# run cmake, optionally directing it to the CUDA toolkit version you have
$ export LDFLAGS="-L${CONDA_PREFIX}/lib"
$ export CC="${CONDA_PREFIX}/bin/x86_64-conda_cos6-linux-gnu-gcc"
$ export CXX="${CONDA_PREFIX}/bin/x86_64-conda_cos6-linux-gnu-g++"
$ export CUDA_VERSION=10.1  # for example
$ cmake .. \
    -DCUDA_TOOLKIT_ROOT_DIR="/usr/local/cuda-${CUDA_VERSION}"
    -DCMAKE_BUILD_TYPE=Release

# if there were no errors, build it
$ make
```
If all went well, the executable should be at `./build/cudaSirecon/cudaSireconDriver`


## Building on Mac

I build with AppleClang 8.1.0.8020042.  Later versions may not be compatible with the NVIDIA compiler.  If during compilation you get an error like `nvcc fatal : The version ('xxxxx') of the host compiler ('Apple clang') is not supported`, then you need to download and install [Command Line Tool for 8.3.2](https://developer.apple.com/download/more/), then run `sudo xcode-select --switch /Library/Developer/CommandLineTools`.  You will also need to have the [NVIDIA CUDA toolkit](https://developer.nvidia.com/cuda-toolkit) installed (I have used versions 8.0 - 10.1 ... so that part doesn't really matter unless you are targetting a GPU with a specific CUDA compute capability). See note below about optionally installing that with conda.

I use [conda](https://docs.conda.io/en/latest/miniconda.html) for the remaining dependencies and build as follows:

```bash
$ conda create -n simbuild -c conda-forge -y cmake liblapack boost-cpp xorg-libx11
$ conda activate simbuild

# as of this writing, versions were
# cmake=3.15.1
# liblapack=3.8.0
# boost-cpp=1.70.0
# xorg-libx11=1.6.8

# optional: if you want to install the CUDA toolkit through conda rather than the NVIDIA website,
# you need to use the dev versions that have the nvcc compiler.
# conda install -c conda-forge conda cudatookit-dev=10.0

# create a build directory inside of CUDA_SIMrecon
$ mkdir build
$ cd build
# run cmake, optionally directing it to the CUDA toolkit version you have
$ cmake .. \
    -DCUDA_TOOLKIT_ROOT_DIR=/Developer/NVIDIA/CUDA-10.1 \
    -DCMAKE_BUILD_TYPE=Release

# if there were no errors, build it
$ make
```

If all went well, the executable should be at `./build/cudaSirecon/cudaSireconDriver`

## Building on Windows

I have built using [Build Tools for Visual Studio](https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2017) versions 2013, 2015, and 2017 (have not tried 2019).  You can install the full Visual Studio Community edition, but the build tools are sufficient.  You will also need to have the [NVIDIA CUDA toolkit](https://developer.nvidia.com/cuda-toolkit) installed (I have used versions 8.0 - 10.1 ... so that part doesn't really matter unless you are targetting a GPU with a specific CUDA compute capability).

I use [conda](https://docs.conda.io/en/latest/miniconda.html) for the remaining dependencies, but it's important that you do this all in a command prompt where you have activated the VS build tools that you have installed... for instance: `x64 Native Tools Command Prompt for VS 2017` or `Visual C++ 2015 x64 Native Build Tools Command Prompt` should be available in the start menu if you have installed VSBuild Tools 2017 or 2015 respectively.

#### Boost

Unforunately, I have not been able to get autolinking to work with the boost libraries in conda, so I resort to manually downloading and compiling boost and placing it at `C:\boost`.  It must be compiled with the same version of visual studio that you are building cudasirecon with:

[Download boost](https://www.boost.org/users/download/) (I am currently using v1.71.0), then prepare the boost libraries as described [here](https://www.boost.org/doc/libs/1_71_0/more/getting_started/windows.html#prepare-to-use-a-boost-library-binary).  Briefly, `cd` into the boost folder you downloaded and unzipped, then run `bootstrap` followed by `.\b2`  (if you get errors, you *may* need to specify the toolset with `bootstrap msvc`).  The compiled libraries will be put into the `./stage` directory.  Then move the whole folder to `C:\` such that you have `C:\boost\stage\`

Finally, 
```bash
> conda create -n simbuild -c conda-forge -y ninja cmake openblas
> conda activate simbuild

# as of this writing, versions were
# cmake=3.15.3
# openblas=0.3.7
# ninja=1.9.0

# create a build directory inside of CUDA_SIMrecon
> mkdir build
> cd build

# run cmake, optionally directing it to the CUDA toolkit version you have
> set CUDA_VERSION=10.1
> set CUDA_TOOLKIT_ROOT_DIR=C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v%CUDA_VERSION%
> cmake .. -G "Ninja" -Wno-dev ^
    -DCMAKE_BUILD_TYPE=Release ^
    -DCUDA_TOOLKIT_ROOT_DIR="%CUDA_TOOLKIT_ROOT_DIR%" 

# if there were no errors, build it
> ninja
```

If all went well, the executable should be at `./build/cudaSirecon/cudaSireconDriver.exe`


## Building `makeotf`


To build the `makeotf` program you also need precompiled fftw-2.x (not 3) libs in a folder called fftw2 in the base directory, with subfolders for each platform that you want to build for (darwin64, linux64, win64), as demonstrated above for the IVE libraries.  See notes below about installing and compiling for mac/linux

The `CMakeLists.txt` file in the root of the project will try to build `makeotf` by default, but you can also build it all manually as follows:

```bash
export SRC_DIR=/path/to/CUDA_SIMrecon  # fix this for your system
export FFTW_ROOT="${SRC_DIR}/fftw2"
export PLATFORM=linux64  # or darwin64 for mac
export IVE_ROOT="${SRC_DIR}/IVE/${PLATFORM}"

# download, compile, and install fftw
wget http://www.fftw.org/fftw-2.1.5.tar.gz
tar -zxvf fftw-2.1.5.tar.gz
cd fftw-2.1.5
./configure --prefix=$FFTW_ROOT --enable-type-prefix --enable-float --enable-threads
make -j 4
make install

# then build makeotf
gcc "${SRC_DIR}/otf/makeotf.c" \
    -I"${IVE_ROOT}/INCLUDE" \
    -I"${FFTW_ROOT}/include" \
    -L"${IVE_ROOT}/LIB" \
    -L"${FFTW_ROOT}/lib" \
    -limlib -lsrfftw -lsfftw -lm \
    -o "${SRC_DIR}/otf/makeotf"
```