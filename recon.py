import os
from contextlib import contextmanager
from subprocess import run

import mrc
import numpy as np

# NOTE: this assumes that you have config files and otf files with the emission wavelength
# names in them... such as `config528.txt` and `otf528.otf`.  You may place those files
# anywhere, but you should update the _DIR variables below

# path to your otf directory.  Defaults to the same as the recon.py file
OTF_DIR = os.path.abspath(os.path.dirname(__file__))
# path to your config directory.  Defaults to the same as the recon.py file
CONFIG_DIR = os.path.abspath(os.path.dirname(__file__))
# name of the reconstruction app
APP = "cudasirecon"


# you could of course do something more complicated with the config/otf finding
def get_otf_and_config(wave, otf_dir=OTF_DIR, config_dir=CONFIG_DIR):
    otf = os.path.join(otf_dir, f"otf{wave}.otf")
    config = os.path.join(config_dir, f"config{wave}.txt")
    if not os.path.exists(otf):
        print(f"Could not find otf file: {otf}")
        return None
    if not os.path.exists(config):
        print(f"Could not find config file: {config}")
        return None
    return otf, config


def merge_files(file_list, outfile, delete=False):
    data = [mrc.imread(fname) for fname in file_list]
    waves = [0, 0, 0, 0, 0]
    for i, item in enumerate(data):
        waves[i] = item.Mrc.hdr.wave[0]
    hdr = data[0].Mrc.hdr
    m = mrc.Mrc2(outfile, mode="w")
    array = np.stack(data, -3)
    m.initHdrForArr(array)
    mrc.copyHdrInfo(m.hdr, hdr)
    m.hdr.NumWaves = len(file_list)
    m.hdr.wave = waves
    m.writeHeader()
    m.writeStack(array)
    m.close()
    if delete:
        try:
            [os.remove(f) for f in file_list]
        except Exception:
            pass
    return outfile


def write_single_channel(array, outpath, hdr, wave):
    """writes a new single-channel file from array data, copying information from hdr"""
    print(f"Extracting channel {wave}...")
    m = mrc.Mrc2(outpath, mode="w")
    m.initHdrForArr(array)
    mrc.copyHdrInfo(m.hdr, hdr)
    m.hdr.NumWaves = 1
    m.hdr.wave = [wave, 0, 0, 0, 0]
    m.writeHeader()
    m.writeStack(array)
    m.close()


@contextmanager
def file_splitter(fname):
    """Context manager that takes care of splitting the file and making sure that
    duplicated data gets cleaned up when done.
    """
    im = mrc.imread(fname)
    header = im.Mrc.header
    namesplit = os.path.splitext(fname)
    try:
        splitfiles = []
        if header.NumWaves == 1:
            # if it's a single channel file, we don't need to split
            wave = header.wave[0]
            otf_config = get_otf_and_config(wave)
            if otf_config is not None:
                splitfiles = [(fname, wave, *otf_config)]
        else:
            # otherwise break out individual wavelenghts
            for c in range(header.NumWaves):
                wave = header.wave[c]
                out = f"{namesplit[0]}-{wave}{namesplit[1]}"
                # assumes channel is the 3rd to last dimension
                data = np.take(im, c, -3)
                otf_config = get_otf_and_config(wave)
                if otf_config is not None:
                    write_single_channel(data, out, header, wave)
                    splitfiles.append((out, wave, *otf_config))
        yield splitfiles
    finally:
        # clean up the split files
        if header.NumWaves > 1:
            print("\nCleaning up extracted channels...")
            for file, b, c, d in splitfiles:
                os.remove(file)


def reconstruct(file, otf, config=None, outfile=None, **kwargs):
    """perform the reconstruction"""
    if outfile is None:
        namesplit = os.path.splitext(file)
        outfile = f"{namesplit[0]}_PROC{namesplit[1]}"
    cmd = [APP, file, outfile, otf]
    if config:
        cmd.extend(["-c", config])
    if kwargs:
        for key, value in kwargs.items():
            cmd.extend([f"--{key}", str(value)])
    print(f"Running command: {' '.join(cmd)}")
    run(cmd)
    return outfile


def reconstructMulti(inFile, outfile=None, **kwargs):
    """Splits multi-channel file into individual channels
    then reconstructs each channel and merges the results
    """
    if outfile is None:
        namesplit = os.path.splitext(inFile)
        outfile = f"{namesplit[0]}_PROC{namesplit[1]}"

    with file_splitter(inFile) as splitfiles:
        processed = []
        for file, wave, otf, config in splitfiles:
            try:
                processed.append(reconstruct(file, otf, config, **kwargs))
            except Exception as e:
                print(f"Cannot reconstruct file {file} due to error {e}")
        if len(processed) > 1:
            print("\nMerging multi-channel reconstructions...")
            final_file = merge_files(processed, outfile, delete=True)
        else:
            final_file = processed[0]
            if final_file != outfile:
                os.rename(final_file, outfile)
    return final_file


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("infile", type=argparse.FileType("r"))
    args, extras = parser.parse_known_args()
    extras = {k: v for k, v in zip(extras[::2], extras[1::2])}
    reconstructMulti(os.path.abspath(args.infile.name), **extras)
