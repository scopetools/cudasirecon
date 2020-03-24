#ifndef SIM_RECONSTRUCTOR_HPP
#define SIM_RECONSTRUCTOR_HPP

#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <boost/tokenizer.hpp>


//! All calculation starts with a SIM_Reconstructor object
/*!
  As seen in cudaSireconDriver.cpp, SIM reconstruction is done in the following steps:
  1. Instantiate a SIM_Reconstructor object by passing along the command line;
  2. For all the time points, do:
  2a. Load raw data and preprocess
  2b. call member processOneVolume() to reconstruct the current time point
  2c. write reconstruction result of current time point
  3. Close file(s) and exit
 */
class SIM_Reconstructor { 

public:
  /*! The constructor does the following:
   * Command line is parsed by this ctor using Boost program_options library;
   * Set all parameters' values;
   * Open raw data and OTF files for reading and output file for writing (except for in
   * TIFF mode output file is not created until the result of current time is getting saved);
   * Allocate all memory buffers.
   */
  SIM_Reconstructor(int argc, char **argv);

  ~SIM_Reconstructor();

  //! Load raw data from a file and do bleach correction (what "rescale" refers to)
  /*!
   * Calls private method loadImageData()
   * As is in other occasions, 'waveIdx' is to indicate color channel but mostly unused.
   * In the future, may add multi-color capability
   */
  void loadAndRescaleImage(int timeIdx, int waveIdx);

  //! The main processing occurs inside this function:
  /*!
    1. Re-initialize all key device buffers under m_reconData;
    2. Fine-tune k0 vectors, modulation amplitudes and phases for all directions and orders
    3. For all directions, pre-filter separated bands, inverse FFT, and assemble the bands
   */
  void processOneVolume();

  //! Off-load processed result to host and save it to disk
  void writeResult(int timeIdx, int waveIdx);

  int getNTimes() { return m_imgParams.ntimes; };
  void setCurTimeIdx(int it) { m_imgParams.curTimeIdx = it; };
  
  ReconParams & getReconParams() {return m_myParams;};
  ImageParams & getImageParams() {return m_imgParams;};

  void closeFiles(); //! only needed for closing MRC output file

  //! Names of input TIFF files (usually a time series whose file names all match a pattern)
  std::vector< std::string > m_all_matching_files;
  CImg<> m_otf_tiff;


private:
  std::string m_config_file;
  ReconParams m_myParams;
  ImageParams m_imgParams;
  ReconData m_reconData;
  DriftParams m_driftParams;
  int m_zoffset;
  po::options_description m_progopts;
  po::variables_map m_varsmap;
  IW_MRC_HEADER m_in_out_header;

  int m_argc;
  char ** m_argv;

  int setupProgramOptions(); //! setup command line options using Boost library
  int setParams();  //! assign parameters after parsing command line
  void setup(CImg<> &inTIFF);
  //! Read header, set up OTF and separation matrix, allocate device buffers, initialize parameters like k0guess etc.
  /*!
    Calls setup_part2() to:
    1. OTF is loaded and transferred to CUDA device by calling getOTFs();
    2. Separation (or un-mixing) matrix is allocated on CUDA device and populated;
    3. Raw image and other buffers, including overlap0/1, bigbuffer, outbuffer, are allocated on CUDA device
       by calling allocateImageBuffers();
    4. Allocate/initialize a bunch of reconstruction related parameters , including:
     a. inscale
     b. k0
     c. k0_time0
     d. k0guess
     e. sum_dir0_phase0
     f. amp
   */
  void setup();
  void openFiles();

  //! Load raw data from disk, do flat-fielding, and transfer data to GPU
  /*!
    This is called by loadAndRescaleImage();
    'zoffset' is used if z-padding is used (almost never)
    'iw' is color channel indicator; rarely used
   */
  void loadImageData(int it, int iw, int zoffset);
};

#endif
