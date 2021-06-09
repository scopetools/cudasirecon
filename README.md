# cudasirecon

Mats Gustafsson & Lin Shao's 3-beam SIM reconstruction software, with CUDA acceleration.

Algorithm as described in [Gustafsson et al (2008) *Biophys*. **94(12)**: 4957–4970. doi: 10.1529/biophysj.107.120345](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2397368/)


## Installation

Packages for Linux and Windows *coming soon* on conda-forge:

```bash
$ conda install -c conda-forge cudasirecon

# you may wish to put it in a new environment ("sim" here):
$ conda create -n sim -y -c conda-forge cudasirecon
$ conda activate sim
```

## Usage

```bash
# the binary will be available as
$ cudasirecon

# for command line help
$ cudasirecon --help

# a typical call for tiff files
$ cudasirecon some/folder file_pattern otf.tif -c config

# a typical call for dv/mrc files
$ cudasirecon data.dv data-PROC.dv otf.dv -c config

# ... where otf.otf was generated using the "makeotf" command:
# note, makeotf not currently working/shipping on windows
$ makeotf /path/to/psf.dv otf.dv -angle -1.855500 -ls 0.2075 -na 1.4 -nimm 1.515 -fixorigin 3 20 -leavekz 7 11 3

# for more help on the makeotf command:
$ makeotf -help
```

### Full list of options/flags

```txt
$ cudasirecon --help

  --input-file arg                  input file (or data folder in TIFF mode)
  --output-file arg                 output file (or filename pattern in TIFF 
                                    mode)
  --otf-file arg                    OTF file
  --usecorr arg                     use the flat-field correction file provided
  --ndirs arg (=3)                  number of directions
  --nphases arg (=5)                number of phases per direction
  --nordersout arg (=0)             number of output orders; must be <= norders
  --angle0 arg (=1.648)             angle of the first direction in radians
  --ls arg (=0.172000006)           line spacing of SIM pattern in microns
  --na arg (=1.20000005)            Detection numerical aperture
  --nimm arg (=1.33000004)          refractive index of immersion medium
  --zoomfact arg (=2)               lateral zoom factor
  --explodefact arg (=1)            artificially exploding the reciprocal-space
                                    distance between orders by this factor
  --zzoom arg (=1)                  axial zoom factor
  --nofilteroverlaps [=arg(=0)]     do not filter the overlaping region between
                                    bands usually used in trouble shooting
  --background arg (=0)             camera readout background
  --wiener arg (=0.00999999978)     Wiener constant
  --forcemodamp arg                 modamps forced to these values
  --k0angles arg                    user given pattern vector k0 angles for all
                                    directions
  --otfRA [=arg(=1)]                using rotationally averaged OTF
  --otfPerAngle [=arg(=1)]          using one OTF per SIM angle
  --fastSI [=arg(=1)]               SIM data is organized in Z->Angle->Phase 
                                    order; default being Angle->Z->Phase
  --k0searchAll [=arg(=0)]          search for k0 at all time points
  --norescale [=arg(=0)]            bleach correcting for z
  --equalizez [=arg(=1)]            bleach correcting for z
  --equalizet [=arg(=1)]            bleach correcting for time
  --dampenOrder0 [=arg(=1)]         dampen order-0 in final assembly
  --nosuppress [=arg(=0)]           do not suppress DC singularity in final 
                                    assembly (good idea for 2D/TIRF data)
  --nokz0 [=arg(=1)]                do not use kz=0 plane of the 0th order in 
                                    the final assembly
  --gammaApo arg (=1)               output apodization gamma; 1.0 means 
                                    triangular apo
  --saveprefiltered arg             save separated bands (half Fourier space) 
                                    into a file and exit
  --savealignedraw arg              save drift-fixed raw data (half Fourier 
                                    space) into a file and exit
  --saveoverlaps arg                save overlap0 and overlap1 (real-space 
                                    complex data) into a file and exit
  -c [ --config ] arg               name of a file of a configuration.
  --2lenses [=arg(=1)]              I5S data
  --bessel [=arg(=1)]               bessel-SIM data
  --besselExWave arg (=0.488000005) Bessel SIM excitation wavelength in microns
  --besselNA arg (=0.143999994)     Bessel SIM excitation NA)
  --deskew arg (=0)                 Deskew angle; if not 0.0 then perform 
                                    deskewing before processing
  --deskewshift arg (=0)            If deskewed, the output image's extra shift
                                    in X (positive->left)
  --noRecon                         No reconstruction will be performed; useful
                                    when combined with --deskew
  --cropXY arg (=0)                 Crop the X-Y dimension to this number; 0 
                                    means no cropping
  --xyres arg (=0.100000001)        x-y pixel size (only used for TIFF files)
  --zres arg (=0.200000003)         z pixel size (only used for TIFF files)
  --zresPSF arg (=0.150000006)      z pixel size used in PSF TIFF files)
  --wavelength arg (=530)           emission wavelength (only used for TIFF 
                                    files)
  --writeTitle [=arg(=1)]           Write command line to image header (may 
                                    cause issues with bioformats)
  -h [ --help ]                     produce help message
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

  
## GPU requirements

This software requires a CUDA-compatible NVIDIA GPU.

The libraries available on conda-forge have been compiled against different
versions of the CUDA toolkit.  The required CUDA libraries are bundled in the
conda distributions so you don't need to install the CUDA toolkit separately.
If desired, you can pick which version of CUDA you'd like based on your needs,
but please note that different versions of the CUDA toolkit have different GPU
driver requirements:

To specify a specific cudatoolkit version, install as follows (for instance, to use
`cudatoolkit=10.2`):

```sh
# NOTE: conda-forge coming soon... not available yet.
conda install -c conda-forge cudasirecon cudatoolkit=10.2
```

| CUDA  | Linux driver | Win driver |
| ----- | ------------ | ---------- |
| 10.2  | ≥ 440.33     | ≥ 441.22   |
| 11.0  | ≥ 450.36.06  | ≥ 451.22   |
| 11.1  | ≥ 455.23     | ≥ 456.38   |
| 11.2  | ≥ 460.27.03  | ≥ 460.82   |


If your CUDA Driver version is too low for the version of cudasirecon that you have installed, you may get an error that looks like: `!!Error occurred: cudaSetDevice failed`

If you run into trouble, feel free to [open an issue](https://github.com/scopetools/cudasirecon/issues) and describe your setup.


## Multichannel reconstruction

`cudasirecon` does not currently accept multi-channel files.  So it is necessary
to temporarily pull out each channel into a new file prior to reconstruction.
The provided `recon.py` script is an example of how to use the
[`mrc`](https://github.com/tlambert03/mrc) package to extract individual
channels from a .dv file, reconstruct them, and merge them back (and clean up
the intermediate file).  It is used as follows (note, `mrc`, `numpy`, and
`cudasirecon` must be in your path):

```
python recon.py /path/to/raw_data.dv
```

`recon.py` will also accept any key value pairs that `cudasirecon` also accepts
(to see that full list, type `cudasirecon -h` at the command prompt).  For
instance, to override just a couple of the reconstruction parameters, you could
do something like this:

```
python recon.py /path/to/raw_data.dv wiener 0.001 background 150
```

There are a couple of hard-coded filepaths in `recon.py`.  Specifically, it
currently expects to find the OTFs and config files in the same directory as the
recon.py script.  You can change that by putting in an absolute directory to
some other folder for the variables at the top of the file:

```python
# path to your otf directory.  Defaults to the same as the recon.py file
OTF_DIR = os.path.abspath(os.path.dirname(__file__))
# path to your config directory.  Defaults to the same as the recon.py file
CONFIG_DIR = os.path.abspath(os.path.dirname(__file__))
```

Note also, that the config and OTF files must be named with the emission
wavelength in the filenames. For example: `config528.txt` and `otf528.otf`


## Compiling locally from source

#### on linux:
```sh
conda env create -f environment-linux.yml
conda activate simbuild
./build.sh  # see build.sh for cmake details

# test it:
./test_data/run.sh
```

#### on windows:
```sh
conda env create -f environment-windows.yml
conda activate simbuild
bld.bat  # see bld.bat for cmake details
```

### Building with MRC support

The IVE/Priism libraries (for MRC/DV support), are not distributed with this
source code and must be acquired seperately from UCSF.  Place them in a folder
called `IVE` at the top level of the source folder (same folder as the
cudaSirecon folder).  It should minimally have the following files and folders
(example shown for linux, use `.a` or `.lib` as necessary for osx or windows)

then build with `cmake -DBUILD_MRC=ON ....`

<details>

<summary>required IVE folder structure</summary>

```
cudasirecon
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

</details>

