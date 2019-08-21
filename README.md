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
$ cudaSireconDriver
# (or a shorter symlink for convenience)
$ sirecon

# for command line help
$ sirecon --help

# a typical call
$ sirecon data.dv data-PROC.dv otf.otf -c config

# ... where otf.otf was generated using the "makeotf" command:
$ makeotf /path/to/psf.dv otf.otf -angle -1.855500 -ls 0.2075 -na 1.4 -nimm 1.515 -fixorigin 3 20 -leavekz 7 11 3

# for more help on the makeotf command:
$ makeotf -help
```

A typical 3D sim config file may look like this:

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

* Currently only accepts images as .dv or .mrc format.
* Requires CUDA-capable NVIDIA GPU and driver.
