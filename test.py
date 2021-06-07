from ctypes import (
    CDLL,
    POINTER,
    Structure,
    c_bool,
    c_char_p,
    c_float,
    c_int,
    c_long,
    c_short,
    c_uint,
    c_ushort,
    c_void_p,
)
from pathlib import Path

import numpy as np
import tifffile as tf

BASE = Path("/home/tjl10/dev/cudasirecon")


SO = str(BASE / "cmake_build/cudaSirecon/libpycudasirecon.so")
CONFIG = str(BASE / "test_data/config-tiff")
img = tf.imread(str(BASE / "test_data/raw.tif"))
print("img: ", img.shape, img.dtype)
lib = CDLL(SO)

_init = lib.SR_new
_init.restype = c_void_p
_init.argtypes = [
    c_int,
    c_int,
    c_int,
    c_char_p,
]

_set_raw = lib.SR_setRaw
_set_raw.restype = None
_set_raw.argtypes = [
    c_void_p,
    np.ctypeslib.ndpointer(c_float, flags="C_CONTIGUOUS"),
    c_int,
    c_int,
    c_int,
]

_get_result = lib.SR_getResult
_get_result.restype = None
_get_result.argtypes = [
    c_void_p,
    np.ctypeslib.ndpointer(c_float, flags="C_CONTIGUOUS"),
]


class ReconParams(Structure):
    _fields_ = [
        ("k0startangle", c_float),
        ("linespacing", c_float),
        ("na", c_float),
        ("nimm", c_float),
        ("ndirs", c_int),
        ("nphases", c_int),
        ("norders_output", c_int),
        ("norders", c_int),
        ("phaseSteps", c_long),
        ("bTwolens", c_int),
        ("bFastSIM", c_int),
        ("bBessel", c_int),
        ("BesselNA", c_float),
        ("BesselLambdaEx", c_float),
        ("deskewAngle", c_float),
        ("extraShift", c_int),
        ("bNoRecon", c_bool),
        ("cropXYto", c_uint),
        ("bWriteTitle", c_int),
        ("zoomfact", c_float),
        ("z_zoom", c_int),
        ("nzPadTo", c_int),
        ("explodefact", c_float),
        ("bFilteroverlaps", c_int),
        ("recalcarrays", c_int),
        ("napodize", c_int),
        ("bSearchforvector", c_int),
        ("bUseTime0k0", c_int),
        ("apodizeoutput", c_int),
        ("apoGamma", c_float),
        ("bSuppress_singularities", c_int),
        ("suppression_radius", c_int),
        ("bDampenOrder0", c_int),
        ("bFitallphases", c_int),
        ("do_rescale", c_int),
        ("equalizez", c_int),
        ("equalizet", c_int),
        ("bNoKz0", c_int),
        ("wiener", c_float),
        ("wienerInr", c_float),
        # there are more ...
    ]

    def __str__(self):
        return (
            f"<ReconParams ndirs={self.ndirs} nphases={self.nphases} "
            f"linespacing={self.linespacing:0.3f} ... >"
        )


class ImageParams(Structure):
    _fields_ = [
        ("nx", c_int),  # image's width after deskewing or same as "nx_raw"
        ("nx_raw", c_int),  # raw image's width before deskewing
        ("ny", c_int),
        ("nz", c_int),
        ("nz0", c_int),
        ("nwaves", c_short),
        ("wave0", c_short),
        ("wave1", c_short),
        ("wave2", c_short),
        ("wave3", c_short),
        ("wave4", c_short),
        ("ntimes", c_short),
        ("curTimeIdx", c_ushort),
        ("dxy", c_float),
        ("dz", c_float),
        ("dz_raw", c_float),
        ("inscale", c_float),
    ]

    def __str__(self):
        return (
            f"<ImageParams nz={self.nz} ny={self.ny} nx={self.nx} nt={self.ntimes} "
            f"dxy={self.dxy:0.3f} dz={self.dz:0.3f} ... >"
        )


_get_recon_params = lib.SR_getReconParams
_get_recon_params.argtypes = [c_void_p]
_get_recon_params.restype = c_void_p

_get_image_params = lib.SR_getImageParams
_get_image_params.argtypes = [c_void_p]
_get_image_params.restype = c_void_p


class SimReconstructor:
    def __init__(self, image: np.ndarray, config=CONFIG) -> None:
        nz, ny, nx = image.shape
        self.image = image
        self.obj = _init(nx, ny, nz, config.encode())
        print(self.get_recon_params())
        print(self.get_image_params())
        self.set_raw(image)
        self.get_result()

    def set_raw(self, array):
        nz, ny, nx = array.shape
        _set_raw(self.obj, array, nx, ny, nz)

    def get_result(self):
        *_, ny, nx = self.image.shape
        nz = 9
        _result = np.empty((nz, ny * 2, nx * 2), dtype=np.float32)
        _get_result(self.obj, _result)
        return _result

    def get_recon_params(self) -> ReconParams:
        return ReconParams.from_address(_get_recon_params(self.obj))

    def get_image_params(self) -> ImageParams:
        return ImageParams.from_address(_get_image_params(self.obj))


sr = SimReconstructor(img)
