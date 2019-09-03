# CUDA_SIMrecon
Lin Shao's 3-beam SIM reconstruction software

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
