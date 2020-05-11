//! Requires C++11

#define cimg_use_tiff
#define cimg_use_cpp11 1
#include "CImg.h"
using namespace cimg_library;

#include <iostream>
#include <cmath>

int main(int argc, char *argv[])
{

  try {
    CImg<> atiff;
    if (argc > 1)
      atiff.assign(argv[1]);
    else {
      std::cerr << "Need at least a file name)\n";
      return 0;
    }

    float gamma = 0.5;
    if (argc > 2)
      gamma = std::stod(std::string(argv[2]));

    CImg<> amplitudes(atiff.width()/2, atiff.height(), atiff.depth());
    cimg_forXYZ(amplitudes, x, y, z) {
      float r = atiff(2*x, y, z);
      float i = atiff(2*x+1, y, z);
      amplitudes(x, y, z) = pow(std::sqrt(r*r + i*i), gamma);
      
    }
    amplitudes.display();
  }
  catch (CImgIOException &e) {
    std::cerr << "file open error\n";
    return -1;
  }
}
