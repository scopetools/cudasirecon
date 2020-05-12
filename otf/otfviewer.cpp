//! Requires C++11

#define cimg_use_tiff
#define cimg_use_cpp11 1
#include "CImg.h"
using namespace cimg_library;

#include <iostream>
#include <string>
#include <IMInclude.h>

#if defined(_MSC_VER) && (_MSC_VER >= 1900)
// for IMLIB with vs >2015
// https://stackoverflow.com/questions/30412951/unresolved-external-symbol-imp-fprintf-and-imp-iob-func-sdl2
FILE _iob[] = {*stdin, *stdout, *stderr};

extern "C" FILE *  __iob_func(void)
{
    return _iob;
}

#endif

int main(int argc, char *argv[])
{
  if (argc < 2) {
    std::cout << "\nUsage: otfviewer OTF_file_name [gamma_value]\n";
    return 0;
  }

  TIFFSetWarningHandler(NULL);
  TIFFSetErrorHandler(NULL);
  TIFF *tf = TIFFOpen(argv[1], "r");
  CImg<> atiff;
  if (tf) {
    TIFFClose(tf);
    atiff.assign(argv[1]);
  }
  else {
    std::cout << "Not a TIFF file; now try MRC\n";
    int istream_no = 1;
    // Suppress IM header printout; somehow only in this mode would the rest of
    // IM calls go well on Windows.
    IMAlPrt(0);

    if (IMOpen(istream_no, argv[1], "ro")) {
      std::cerr << "File " << argv[1] << " cannot be opened.\n";
      return -1;
    }
    int ixyz[3], mxyz[3], nxyzst[3], pixeltype;
    float mmm[3];
    IW_MRC_HEADER header;
    
    IMRdHdr(istream_no, ixyz, mxyz, &pixeltype, mmm, mmm+1, mmm+2);
    IMGetHdr(istream_no, &header);

    atiff.assign(header.nx * 2, header.ny, header.nz);
    for (int z=0; z<header.nz; z++)
      IMRdSec(istream_no, atiff.data(0, 0, z));
    IMClose(istream_no);
  }

  float gamma = 0.5;
  if (argc > 2)
    // gamma = std::stod(std::string(argv[2]));
    gamma = std::stod(argv[2]);

  CImg<> amplitudes(atiff.width()/2, atiff.height(), atiff.depth());
  cimg_forXYZ(amplitudes, x, y, z) {
    float r = atiff(2*x, y, z);
    float i = atiff(2*x+1, y, z);
    amplitudes(x, y, z) = pow(std::sqrt(r*r + i*i), gamma);
  }

  amplitudes.display(argv[1], false);

  // Somehow doubled XY image looks the same:
  // CImg<> large = amplitudes.resize_doubleXY();
  // std::cout << amplitudes.width() << "  " << large.width() << std::endl;
  // large.display(argv[1], false);
  return 0;
}
