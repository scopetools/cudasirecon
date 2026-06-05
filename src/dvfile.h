#ifndef DVFILE_H
#define DVFILE_H

#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

// Data Types
#define IW_AS_IS -1
#define IW_BYTE 0
#define IW_SHORT 1
#define IW_FLOAT 2
#define IW_COMPLEX_SHORT 3
#define IW_COMPLEX 4
#define IW_EMTOM 5
#define IW_USHORT 6
#define IW_LONG 7

// image sequence definitions
#define ZTW_SEQUENCE 0 /* "non-interleaved", by defintion       */
#define WZT_SEQUENCE 1 /* "interleaved", from R3D and others    */
#define ZWT_SEQUENCE 2 /* new sequence. Unsupported as of 11/97 */

// Possible values for the file type
#define IM_NORMAL_IMAGES 0
#define IM_TILT_SERIES 1
#define IM_STEREO_TILT_SERIES 2
#define IM_AVERAGED_IMAGES 3
#define IM_AVERAGED_STEREO_PAIRS 4
#define IM_EM_TILT_SERIES 5
#define IM_MULTIPOSITION 20
#define IM_PUPIL_FUNCTION 8000

static const std::unordered_map<int, size_t> pixelTypeSizes = {
    {IW_BYTE, sizeof(uint8_t)},      {IW_SHORT, sizeof(int16_t)},
    {IW_FLOAT, sizeof(float)},       {IW_COMPLEX_SHORT, 2 * sizeof(int16_t)},
    {IW_COMPLEX, 2 * sizeof(float)}, {IW_EMTOM, sizeof(int16_t)},
    {IW_USHORT, sizeof(uint16_t)},   {IW_LONG, sizeof(int32_t)}};

typedef struct IW_MRC_Header {
  int32_t nx, ny, nz;        // nz=num_planes*num_waves*num_times
  int32_t mode;              // data type
  int32_t nxst, nyst, nzst;  // index of the first col/row/section
  int32_t mx, my, mz;        // number of intervals in x/y/z
  float xlen, ylen, zlen;    // pixel spacing for x/y/z
  float alpha, beta, gamma;  // cell angles
  int32_t mapc, mapr, maps;  // column/row/section axis
  float amin, amax, amean;   // min/max/mean intensity
  int32_t ispg, inbsym;      // space group number, number of bytes in extended header
  int16_t nDVID, nblank;     // ID value, unused
  int32_t ntst;              // starting time index (used for time series data)
  char ibyte[24];            // 24 bytes of blank space
  int16_t nint, nreal;       // number of integers/floats in extended header per section
  int16_t nres, nzfact;      // number of sub-resolution data sets, reduction quotient for z axis
  float min2, max2, min3, max3, min4, max4;  // min/max intensity for 2nd, 3rd, 4th wavelengths
  int16_t file_type, lens, n1, n2, v1, v2;   // file type, lens ID, n1, n2, v1, v2
  float min5, max5;                          // min/max intensity for 5th wavelength
  int16_t num_times;                         // number of time points
  int16_t interleaved;                       // (0 = ZTW, 1 = WZT, 2 = ZWT)
  float tilt_x, tilt_y, tilt_z;              // x/y/z axis tilt angles
  int16_t num_waves, iwav1, iwav2, iwav3, iwav4, iwav5;  // number & values of wavelengths
  float zorig, xorig, yorig;                             // z/x/y origin
  int32_t nlab;                                          // number of titles
  char label[800];

  std::string sequence_order() const {
    switch (interleaved) {
      case ZTW_SEQUENCE: return "ZTW";
      case WZT_SEQUENCE: return "WZT";
      case ZWT_SEQUENCE: return "ZWT";
      default: return "ZTW";
    }
  }

  int num_planes() const { return nz / (num_waves ? num_waves : 1) / (num_times ? num_times : 1); }

  std::string image_type() const {
    switch (file_type) {
      case IM_NORMAL_IMAGES: return "NORMAL";
      case 100: return "NORMAL";
      case IM_TILT_SERIES: return "TILT_SERIES";
      case IM_STEREO_TILT_SERIES: return "STEREO_TILT_SERIES";
      case IM_AVERAGED_IMAGES: return "AVERAGED_IMAGES";
      case IM_AVERAGED_STEREO_PAIRS: return "AVERAGED_STEREO_PAIRS";
      case IM_EM_TILT_SERIES: return "EM_TILT_SERIES";
      case IM_MULTIPOSITION: return "MULTIPOSITION";
      case IM_PUPIL_FUNCTION: return "PUPIL_FUNCTION";
      default: return "UNKNOWN";
    }
  }

  void print() {
    std::cout << "Header:" << std::endl;
    std::cout << "  Dimensions: " << ny << "x" << nx << "x" << num_planes() << std::endl;
    std::cout << "  Number of wavelengths: " << num_waves << std::endl;
    std::cout << "  Number of time points: " << num_times << std::endl;
    std::cout << "  Pixel size: " << mode << std::endl;
    std::cout << "  Pixel spacing: " << xlen << "x" << ylen << "x" << zlen << std::endl;
    std::cout << "  mxyz: " << mx << "x" << my << "x" << mz << std::endl;
    std::cout << "  Min/Max/Mean: " << amin << "/" << amax << "/" << amean << std::endl;
    std::cout << "  Image type: " << image_type() << std::endl;
    std::cout << "  Sequence order: " << sequence_order() << std::endl;
  }

} IW_MRC_HEADER, *IW_MRC_HEADER_PTR;

class DVFile {
 private:
  std::unique_ptr<std::fstream> _file;
  std::string _path;
  bool _big_endian;
  IW_MRC_Header hdr;
  bool closed = true;
  // IVE "ConversionFlag": when true (the default), integer pixel types are
  // converted to/from float on read/write (the caller works in float).  Toggled
  // per stream by IMAlCon.
  bool _convert = true;

  // Private default constructor
  DVFile() = default;

 public:
  // Constructor for opening an existing file
  DVFile(const std::string& path) {
    _path = path;
    _file = std::make_unique<std::fstream>(path, std::ios::in | std::ios::out | std::ios::binary);
    if (!_file->is_open()) {
      throw std::runtime_error("Failed to open file");
    }

    // Determine byte order
    _file->seekg(24 * 4);
    char dvid[2];
    _file->read(dvid, 2);
    if (dvid[0] == (char)0xA0 && dvid[1] == (char)0xC0) {
      _big_endian = false;
    } else if (dvid[0] == (char)0xC0 && dvid[1] == (char)0xA0) {
      _big_endian = true;
      // This reader does not byte-swap; a big-endian file would be silently
      // mis-read (header fields and pixels).  All modern DeltaVision data is
      // little-endian, so reject big-endian explicitly rather than corrupt it.
      throw std::runtime_error(
          path + " is big-endian, which this reader does not support.");
    } else {
      throw std::runtime_error(path + " is not a recognized DV file.");
    }

    // Read header
    _file->seekg(0);
    _file->read(reinterpret_cast<char*>(&hdr), sizeof(IW_MRC_Header));
    // Position the read/write pointer at the first pixel section, i.e. just past
    // the extended header (inbsym bytes).  Sequential reads (bare IMRdSec, used to
    // load OTF and PSF files) rely on the stream starting at the data; IVE leaves
    // the file positioned there after IMOpen/IMRdHdr.  Without this, files that
    // carry an extended header (real .otf/.dv from the microscope) are read
    // starting inside the extended header -> garbage/zeros.
    _seekToData();
    closed = false;
  }

  static std::unique_ptr<DVFile> createNew(const std::string& path) {
    std::unique_ptr<DVFile> dvfile(new DVFile());  // Use the private default constructor
    dvfile->_path = path;
    dvfile->_file = std::make_unique<std::fstream>(
        path, std::ios::in | std::ios::out | std::ios::binary | std::ios::trunc);

    if (!dvfile->_file->is_open()) {
      throw std::runtime_error("Failed to create file");
    }

    dvfile->closed = false;
    return dvfile;
  }

  // Method to read the extended header
  void readExtendedHeader() {
    if (hdr.inbsym <= 0) {
      throw std::runtime_error("No extended header present.");
    }

    // Calculate the size of the extended header
    int num_sections = hdr.nz;
    size_t extended_header_size = (hdr.nint + hdr.nreal) * num_sections * sizeof(int32_t);

    // Allocate buffer for the extended header
    size_t bufferSize = extended_header_size / sizeof(int32_t);
    std::vector<int32_t> extended_header(bufferSize);

    // Read the extended header
    _file->seekg(sizeof(IW_MRC_Header));
    _file->read(reinterpret_cast<char*>(extended_header.data()), extended_header_size);

    if (!_file->good()) {
      throw std::runtime_error("Failed to read the extended header.");
    }

    // Optionally, process the extended header based on your specific format
  }

  /**
   * @brief Calculate the offset of a section in the file.
   *
   * 0 (C/C++ macro is ZTW_SEQUENCE)
   * This is the normal arrangement for processed images. Sometimes referred to as
   * "non-interleaved". is = iz + nz * (it + nt * iw) Reversing that gives iz = is modulo nz, iw =
   * is / (nz * nt), and it = (is / nz) modulo nt. 1 (C/C++ macro is WZT_SEQUENCE) This is the
   * typical arrangement for images acquired from a microscope. Sometimes referred to as
   * "interleaved" is = iw + nw * (iz + nz * it). Reversing that gives iz = (is / nw) modulo nz, iw
   * = is modulo nw, and it = is / ( nw * nz). 2 (C/C++ macro is ZWT_SEQUENCE) Although not widely
   * used, ZWT will find uses with certain processing algorithms. is = iz + nz * (iw + nw * it).
   * Reversing that gives iz = is modulo nz, iw = (is / nz) modulo nw, and it = is / (nz * nw).
   */
  int sectionOffset(int iz, int iw, int it) {
    _validateZWT(iz, iw, it);
    int is;
    int nz = hdr.num_planes();
    int nt = hdr.num_times;
    int nw = hdr.num_waves;
    switch (hdr.interleaved) {
      case WZT_SEQUENCE: is = iw + nw * (iz + nz * it); break;
      case ZWT_SEQUENCE: is = iz + nz * (iw + nw * it); break;
      case ZTW_SEQUENCE:
      default: is = iz + nz * (it + nt * iw); break;
    }
    return is;
  }

  // Set the read/write file pointer to the beginning of a section
  void setCurrentZWT(int z, int w, int t) {
    int offset = sizeof(IW_MRC_Header) + hdr.inbsym + sectionOffset(z, w, t) * _sectionSize();
    _file->seekg(offset);
    _file->seekp(offset);
  }

  void setConvert(bool convert) { _convert = convert; }

  // Read the per-section extended-header values for section (z,w,t).  The DV
  // extended header stores, for each section in ZWT order, `nint` int32 values
  // followed by `nreal` float32 values (the standard DV float fields are
  // photosensor, timestamp, stageX/Y/Z, min/max/mean, expTime, ...).  `ival` and
  // `rval` must have room for `nint` ints and `nreal` floats respectively.  The
  // current get-pointer is preserved so this can be called between sequential
  // section reads.
  void readExtHdrZWT(int z, int w, int t, int* ival, float* rval) {
    int ni = hdr.nint, nr = hdr.nreal;
    if (hdr.inbsym <= 0 || (ni <= 0 && nr <= 0)) return;  // no extended header
    int is = sectionOffset(z, w, t);
    std::streamoff rec = static_cast<std::streamoff>(ni + nr) * 4;
    if ((static_cast<std::streamoff>(is) + 1) * rec > hdr.inbsym) return;  // out of range
    std::streamoff off = static_cast<std::streamoff>(sizeof(IW_MRC_Header)) +
                         static_cast<std::streamoff>(is) * rec;
    _file->clear();  // sequential reads may have left eofbit set
    std::streampos saved = _file->tellg();
    _file->seekg(off);
    if (ni > 0) _file->read(reinterpret_cast<char*>(ival), ni * 4);
    if (nr > 0) _file->read(reinterpret_cast<char*>(rval), nr * 4);
    _file->clear();
    if (saved >= 0) _file->seekg(saved);
  }

  void readSec(void* array) {
    if (closed) {
      throw std::runtime_error("Cannot read from closed file. Please reopen with .open()");
    }
    if (!_convert || hdr.mode == IW_FLOAT || hdr.mode == IW_COMPLEX) {
      // already the working type (or conversion disabled): read raw bytes
      _file->read(reinterpret_cast<char*>(array), _sectionSize());
      return;
    }
    // integer storage with ConversionFlag on: read native type, convert to float
    size_t npix = static_cast<size_t>(hdr.ny) * hdr.nx;
    float* out = reinterpret_cast<float*>(array);
    switch (hdr.mode) {
      case IW_BYTE:          _readAs<uint8_t>(out, npix); break;
      case IW_SHORT:         _readAs<int16_t>(out, npix); break;
      case IW_EMTOM:         _readAs<int16_t>(out, npix); break;
      case IW_USHORT:        _readAs<uint16_t>(out, npix); break;
      case IW_LONG:          _readAs<int32_t>(out, npix); break;
      case IW_COMPLEX_SHORT: _readAs<int16_t>(out, 2 * npix); break;
      default:
        throw std::runtime_error("Unsupported pixel mode for conversion: " +
                                 std::to_string(hdr.mode));
    }
  }

  void readSec(void* array, int t, int w, int z) {
    setCurrentZWT(z, w, t);
    readSec(array);
  }

  // call setCurrentZWT to set the z/w/t plane before calling this
  void writeSection(const void* array) {
    if (closed) {
      throw std::runtime_error("Cannot write to closed file. Please reopen with .open()");
    }
    if (!_convert || hdr.mode == IW_FLOAT || hdr.mode == IW_COMPLEX) {
      _file->write(reinterpret_cast<const char*>(array), _sectionSize());
      return;
    }
    // float working data with ConversionFlag on: convert to the stored int type
    size_t npix = static_cast<size_t>(hdr.ny) * hdr.nx;
    const float* in = reinterpret_cast<const float*>(array);
    switch (hdr.mode) {
      case IW_BYTE:          _writeAs<uint8_t>(in, npix); break;
      case IW_SHORT:         _writeAs<int16_t>(in, npix); break;
      case IW_EMTOM:         _writeAs<int16_t>(in, npix); break;
      case IW_USHORT:        _writeAs<uint16_t>(in, npix); break;
      case IW_LONG:          _writeAs<int32_t>(in, npix); break;
      case IW_COMPLEX_SHORT: _writeAs<int16_t>(in, 2 * npix); break;
      default:
        throw std::runtime_error("Unsupported pixel mode for conversion: " +
                                 std::to_string(hdr.mode));
    }
  }

  size_t getPixelSize() { return pixelTypeSizes.at(hdr.mode); }

  void open() {
    if (closed) {
      _file->open(_path, std::ios::binary);
      if (!_file->is_open()) {
        throw std::runtime_error("Failed to open file");
      }
      // Reopened at offset 0; reposition at the first data section so sequential
      // (bare) reads/writes start at the pixel data, not the header.
      _seekToData();
      closed = false;
    }
  }

  void close() {
    if (!closed) {
      _file->close();
      closed = true;
    }
  }

  ~DVFile() { close(); }

  std::string getPath() const { return _path; }

  IW_MRC_Header getHeader() const { return hdr; }

  void putHeader(const IW_MRC_Header& header) {
    if (closed) {
      throw std::runtime_error("Cannot write to closed file. Please reopen with .open()");
    }
    _file->seekp(0);
    _file->write(reinterpret_cast<const char*>(&header), sizeof(IW_MRC_Header));
    hdr = header;
    // Position writes at the first data section (past any extended header) so
    // subsequent sequential IMWrSec calls land on the pixel data.
    _seekToData();
  }

  bool isClosed() const { return closed; }

 private:
  // Byte offset of the first pixel section: main header + extended header.
  std::streamoff _dataOffset() const {
    return static_cast<std::streamoff>(sizeof(IW_MRC_Header)) +
           (hdr.inbsym > 0 ? hdr.inbsym : 0);
  }

  // Position both get/put pointers at the first pixel section (past the ext hdr).
  void _seekToData() {
    _file->seekg(_dataOffset());
    _file->seekp(_dataOffset());
  }

  // Return the size of a frame/section in bytes
  size_t _sectionSize() { return hdr.ny * hdr.nx * getPixelSize(); }

  // Read `n` elements of stored type T and convert to float into `out`.
  template <typename T>
  void _readAs(float* out, size_t n) {
    std::vector<T> tmp(n);
    _file->read(reinterpret_cast<char*>(tmp.data()), n * sizeof(T));
    for (size_t i = 0; i < n; ++i) out[i] = static_cast<float>(tmp[i]);
  }

  // Convert `n` floats in `in` to stored type T and write them.  Integer targets
  // are rounded to nearest (matching IVE), not truncated toward zero.
  template <typename T>
  void _writeAs(const float* in, size_t n) {
    std::vector<T> tmp(n);
    for (size_t i = 0; i < n; ++i)
      tmp[i] = static_cast<T>(std::lround(in[i]));
    _file->write(reinterpret_cast<const char*>(tmp.data()), n * sizeof(T));
  }

  void _validateZWT(int z, int w, int t) {
    if (t >= hdr.num_times) {
      throw std::runtime_error("Time index out of range");
    }
    if (w >= hdr.num_waves) {
      throw std::runtime_error("Wavelength index out of range");
    }
    if (z >= hdr.num_planes()) {
      throw std::runtime_error("Section index out of range");
    }
  }
};

//////////////////////////////////////////////////////////////////////////////
// IVE API
//////////////////////////////////////////////////////////////////////////////

// stream -> DVFile
static std::map<int, std::unique_ptr<DVFile>> dvfile_map;

inline DVFile& getDVFile(int istream) {
  auto it = dvfile_map.find(istream);
  if (it == dvfile_map.end()) {
    throw std::runtime_error("Stream not found: " + std::to_string(istream));
  }
  return *(it->second);
}

/**
 * @brief Open an image file and attach it to a stream.
 *
 * @param istream The input stream to be used for the operation.
 * @param name The name of the file to be opened.
 * @param attrib The file mode.
 *  - "ro" Opens an existing file or image window for reading only.
 *  - "new" Creates a file or image window and opens it for reading and writing.
 * @return int 0 if successful and 1 if not.
 */
inline int IMOpen(int istream, const char* name, const char* attrib) {
  // Check if the stream identifier is already in use and close it if necessary
  if (dvfile_map.find(istream) != dvfile_map.end()) {
    dvfile_map[istream]->close();
    dvfile_map.erase(istream);
    std::cerr << "Warning: Reusing stream identifier " << istream << ". Previous stream closed."
              << std::endl;
  }

  if (std::string(attrib) == "ro") {
    try {
      dvfile_map[istream] = std::make_unique<DVFile>(name);
    } catch (const std::exception& e) {
      std::cerr << "Error: " << e.what() << std::endl;
      return 1;  // Return non-zero to indicate failure
    }
  } else if (std::string(attrib) == "new") {
    try {
      dvfile_map[istream] = DVFile::createNew(name);
    } catch (const std::exception& e) {
      std::cerr << "Error: " << e.what() << std::endl;
      return 1;  // Return non-zero to indicate failure
    }
  } else {
    std::cerr << "Unknown file mode: " << attrib << std::endl;
    return 1;  // Return non-zero to indicate failure
  }
  return 0;  // Return 0 to indicate success
}

inline void IMClose(int istream) {
  // call destructor of DVFile
  dvfile_map.erase(istream);
}

inline void IMGetHdr(int istream, IW_MRC_HEADER* header) {
  *header = getDVFile(istream).getHeader();
}

inline void IMRdHdr(int istream, int ixyz[3], int mxyz[3], int* imode, float* min, float* max,
                    float* mean) {
  IW_MRC_HEADER header;
  IMGetHdr(istream, &header);
  ixyz[0] = header.nx;
  ixyz[1] = header.ny;
  ixyz[2] = header.nz;
  mxyz[0] = header.mx;
  mxyz[1] = header.my;
  mxyz[2] = header.mz;
  *imode = header.mode;
  *min = header.amin;
  *max = header.amax;
  *mean = header.amean;
}

/**
 * @brief Set the image conversion mode during read/write operations from image storage.
 *
 * By default in IVE, images read from storage are converted to 4-byte float, and
 * images written are converted to the stream's stored data type. The default is
 * ConversionFlag=TRUE.  We honor that: with the flag on, integer pixel types are
 * converted to/from float on read/write; with it off, raw bytes are used.
 *
 * @param istream The input stream to be used for the operation.
 * @param flag Nonzero => ConversionFlag ON (convert); zero => raw I/O.
 */
inline void IMAlCon(int istream, int flag) {
  getDVFile(istream).setConvert(flag != 0);
}

/**
 * @brief Change the image titles.
 *
 * @param istream The input stream to be used for the operation.
 * @param titles The titles to be changed.  Contains at least NumTitles title strings, each of
 * which must contain exactly 80 characters.
 * @param num_titles The number of titles to be changed.
 */
inline void IMAlLab(int istream, const char* labels, int nl) {
  // The DV header holds up to 10 title records of 80 chars each (label[800]).
  // `labels` is nl contiguous 80-char records.
  DVFile& dvfile = getDVFile(istream);
  IW_MRC_HEADER header = dvfile.getHeader();
  if (nl < 0) nl = 0;
  if (nl > 10) nl = 10;
  std::memset(header.label, ' ', sizeof(header.label));
  if (nl > 0) std::memcpy(header.label, labels, static_cast<size_t>(nl) * 80);
  header.nlab = nl;
  dvfile.putHeader(header);
}

/**
 * @brief  Enable or disable printing to standard output ("stdout").
 *
 * Certain IM functions will print information to stdout if Format=TRUE, which is the default. To
 * disable printing, set flag=FALSE.
 *
 * @param flag The flag indicating the type of operation to be performed.
 */
inline void IMAlPrt(int flag) {
  if (flag == 1) {
    std::cerr << "Warning: IMAlPrt is not implemented." << std::endl;
  }
}

/**
 * @brief Position the read/write point at a particular Z, W, T section.
 *
 * If the stream points to a scratch window, IMPosnZWT can only change the destination wavelength.
 *
 * @param istream The input stream to be used for the operation.
 * @param iz The Z section number.
 * @param iw The wavelength number.
 * @param it The time-point number.
 * @return int 0 if successful and 1 if not.
 */
inline int IMPosnZWT(int istream, int iz, int iw, int it) {
  DVFile& dvfile = getDVFile(istream);

  try {
    dvfile.setCurrentZWT(iz, iw, it);
    return 0;
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}

/**
 * @brief Read the next section.
 *
 * Reads the next section into ImgBuffer and advances the file pointer to the
 * section after that. The results are undefined if ImgBuffer does not have at
 * least nx * ny elements or the file pointer does not point to the beginning of
 * a section.
 * In most cases, ImgBuffer will contain floating-point data. When
 * image conversion is off, however, the data type of ImgBuffer should
 * correspond to whatever data type is actually stored. For example, if the
 * image data are stored as 16-bit integers, then ImgBuffer should point to a
 * 16-bit buffer. See IMAlCon.
 *
 * @param istream The input stream to be used for the operation.
 * @param ImgBuffer The image buffer to store the data.
 */
inline void IMRdSec(int istream, void* ImgBuffer) {
  try {
    getDVFile(istream).readSec(ImgBuffer);
  } catch (const std::runtime_error& e) {
    std::cerr << "Error reading section: " << e.what() << std::endl;
    // Handle the error appropriately, e.g., by rethrowing or returning an error code
    throw;  // Rethrow the exception to propagate it further
  }
}

inline void IMWrSec(int istream, const void* array) {
  try {
    getDVFile(istream).writeSection(array);
  } catch (const std::runtime_error& e) {
    std::cerr << "Error writing section: " << e.what() << std::endl;
    // Handle the error appropriately, e.g., by rethrowing or returning an error code
    throw;  // Rethrow the exception to propagate it further
  }
}

/**
 * @brief Put an entire header into a stream.
 *
 * Header should point to a memory location that contains a complete image header.
 *
 * @param istream The input stream to be used for the operation.
 * @param header The header to be saved in the stream.
 */
inline void IMPutHdr(int istream, const IW_MRC_HEADER* header) {
  getDVFile(istream).putHeader(*header);
}

/**
 * @brief Write the image header to the storage device.
 *
 * Write image header associated with StreamNum to the storage device. Use this
 * function to save the results of all IMAl functions. Header modifications are
 * not saved until IMWrHdr is used! The contents of Title will be saved in the
 * header according to the method indicated by ntflag. The minimum, maximum, and
 * mean intensity - Min, Max, and Mean, respectively - of wavelength 0 are also
 * saved in the header every time IMWrHdr is used.
 *
 * @param istream The input stream to be used for the operation.
 * @param title The title to be saved in the header.
 * @param ntflag The method to be used to save the title.
 *  - 0: use Title as the only title
 *  - 1: add Title to the end of the list
 */
inline void IMWrHdr(int istream, const char title[80], int ntflag, float dmin, float dmax,
                    float dmean) {
  DVFile& dvfile = getDVFile(istream);
  IW_MRC_HEADER header = dvfile.getHeader();
  header.amin = dmin;
  header.amax = dmax;
  header.amean = dmean;
  int nlab = header.nlab;
  if (nlab < 0) nlab = 0;
  if (nlab > 10) nlab = 10;
  if (ntflag == 0) {
    // use `title` as the only title (one 80-char record)
    std::memset(header.label, ' ', sizeof(header.label));
    std::memcpy(header.label, title, 80);
    header.nlab = 1;
  } else if (ntflag == 1) {
    // append `title` as a new 80-char record, if there is room (max 10)
    if (nlab < 10) {
      std::memcpy(header.label + static_cast<size_t>(nlab) * 80, title, 80);
      header.nlab = nlab + 1;
    }
  } else {
    throw std::runtime_error("Invalid ntflag: " + std::to_string(ntflag));
  }

  dvfile.putHeader(header);
}

/**
 * @brief Return extended header values for a particular Z section, wavelength,
 * and time-point.
 *
 * The integer and floating-point values for the requested Z section (ZSecNum),
 * wavelength (WaveNum), and time-point (TimeNum) are returned in IntValues and
 * FloatValues, respectively.
 *
 */
inline void IMRtExHdrZWT(int istream, int iz, int iw, int it, int ival[], float rval[]) {
  // `ival`/`rval` must hold at least nint ints and nreal floats (the per-section
  // extended-header sizes from the file header); the caller is responsible for
  // sizing them accordingly.
  getDVFile(istream).readExtHdrZWT(iz, iw, it, ival, rval);
}

#endif  // DVFILE_H