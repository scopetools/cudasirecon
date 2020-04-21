#######################################################################
# Build instruction for CUDA SIMrecon project
#
# 1. Prerequisites
#  1.1 Cmake (> 2.8.11)
#  1.2 CUDA SDK (>=7.0 on Mac OS X, otherwise >=5.0)
#  1.3a. If Linux or Mac OS X, GCC (version suitable to the CUDA version being used; read CUDA manuals to find out)
#  1.3b. If Windows, Visual C++ (>= 2012)
#  1.4 Boost libraries (> 1.48)
#  1.5 If Linux or Mac, the development (not just the runtime) libraries of LAPACK, BLAS, and X11. For Mac, you can obtain these libraries and other  useful stuff by using tools like macports, homebrew, or fink.
#  1.6 If Linux or Mac, Priism/IVE libraries (download from http://msg.ucsf.edu/IVE/Download/index.html; follow the installation instructions and put a line ". install_folder/Priism_setup.sh" in your .bashrc file, where install_folder is where the Priism package is installed and retart a BASH shell. To test if the environment is set up right, run "echo $IVE_BASE" and see if it prints out the installation folder name.
#
# 2. Make a build dir (assumed to be "build" under the main folder where this 
#    CMakeLists.txt is located) and cd into it
#
# 3. On Linux or Mac at a shell prompt, type:
#  $ cmake -D CMAKE_BUILD_TYPE=Release ..
#  On Windows, one extra flag is needed. Type the following in a 
#  VS build environment:
#  > cmake -D CMAKE_BUILD_TYPE=Release -G "NMake Makefiles" ..
#  Make sure there isn't any error or warning messages. 
#  Always delete everything in the build dir before re-run cmake
#
#  3.1 The nvcc flag specified here only works for CUDA compute-capability 6.1. If you need to include other compute-capabilities, add '-D CUDA_NVCC_FLAGS="-gencode arch=compute_xx,code=sm_xx" ' to the "cmake ..." command.
#
# 4. Type "make" on Linux/Mac or "nmake" on Window to build the executable
#
# 5. If building is successful, an executable cudaSireconDriveris generated in
#  build/cudaSirecon. Run the test to make sure it works:
#  $ ../testData/job.sh
#
#######################################################################


cmake_minimum_required (VERSION 2.8)

project (cudaSIMRecon)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

ENABLE_TESTING()

file(READ "version.txt" ver)

string(REGEX MATCH "VERSION_MAJOR ([0-9]*)" _ ${ver})
set(ver_major ${CMAKE_MATCH_1})

string(REGEX MATCH "VERSION_MINOR ([0-9]*)" _ ${ver})
set(ver_minor ${CMAKE_MATCH_1})

string(REGEX MATCH "VERSION_PATCH ([0-9]*)" _ ${ver})
set(ver_patch ${CMAKE_MATCH_1})

message("version: ${ver_major}.${ver_minor}.${ver_patch}")
set(PROJECT_VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})


######################################################################
#
# Set permissions before adding subdirectories
#
######################################################################
set(SCI_GROUP_WRITE GROUP_WRITE)
set(SCI_WORLD_FILE_PERMS WORLD_READ)
set(SCI_WORLD_PROGRAM_PERMS WORLD_READ WORLD_EXECUTE)

find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if (APPLE)
# Apple's default clang complier does not support OpenMP
#  set(CMAKE_C_COMPILER gcc)
#  set(CMAKE_CXX_COMPILER g++)
#  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -undefined dynamic_lookup -flat_namespace")
#  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -undefined dynamic_lookup -flat_namespace")
endif()


# Building the documentation
option(BUILD_DOCUMENTATION "Use Doxygen to create the HTML based API documentation" OFF)
if(BUILD_DOCUMENTATION)
  FIND_PACKAGE(Doxygen)
  if (NOT DOXYGEN_FOUND)
    message(FATAL_ERROR 
      "Doxygen is needed to build the documentation. Please
      install it correctly")
    endif()
    #-- Configure the Template Doxyfile for our specific project
    configure_file(Doxyfile.in
      ${PROJECT_BINARY_DIR}/Doxyfile  @ONLY IMMEDIATE)
    #-- Add a custom target to run Doxygen when ever the project is built
    add_custom_target (Docs ALL 
      COMMAND ${DOXYGEN_EXECUTABLE}
      ${PROJECT_BINARY_DIR}/Doxyfile
      SOURCES
      ${PROJECT_BINARY_DIR}/Doxyfile)
endif()

find_package(CUDA)

find_package(TIFF)

if(WIN32)
  set (Boost_USE_STATIC_LIBS    ON)
  set (Boost_USE_MULTITHREADED  ON)
  set (Boost_USE_STATIC_RUNTIME OFF)

  if ( ${MSVC_VERSION} GREATER 1800)
    find_library(LEG_STDIO legacy_stdio_definitions)
    message(STATUS "adding legacy stdio definitions for >VS14.0 " ${LEG_STDIO})
  endif()

endif(WIN32)

find_package( Boost REQUIRED COMPONENTS program_options filesystem system regex)
message(STATUS "Boost_INCLUDE_DIRS " ${Boost_INCLUDE_DIRS})
message(STATUS "Boost_LIBRARY_DIRS " ${Boost_LIBRARY_DIRS})
message(STATUS "Boost_LIBRARIES " ${Boost_LIBRARIES})

if(WIN32 OR APPLE)
  link_directories ( ${Boost_LIBRARY_DIRS} )
  include_directories(  ${Boost_INCLUDE_DIRS} )
endif()

if(APPLE)
  find_package(X11)
  link_directories ( ${X11_INCLUDE_DIR}/../lib )
endif()

set(CUDA_NVCC_FLAGS
  -D_USE_MATH_DEFINES
  -gencode=arch=compute_50,code=sm_50; 
  -gencode=arch=compute_52,code=sm_52; 
  -gencode=arch=compute_60,code=sm_60; 
  -gencode=arch=compute_61,code=sm_61; 
  ${CUDA_NVCC_FLAGS} 
)

if(${CUDA_VERSION_MAJOR} GREATER 8)
  set(CUDA_NVCC_FLAGS
    ${CUDA_NVCC_FLAGS};
    -gencode=arch=compute_70,code=sm_70;
  )
endif()

if((${CUDA_VERSION_MAJOR} GREATER 9) AND (NOT APPLE))
  set(CUDA_NVCC_FLAGS
    ${CUDA_NVCC_FLAGS};
    -gencode=arch=compute_75,code=sm_75;
    -gencode=arch=compute_75,code=compute_75
  )
endif()

message(STATUS "CUDA_NVCC_FLAGS: " ${CUDA_NVCC_FLAGS})

add_subdirectory(Buffers)
add_subdirectory(cudaSirecon)
#if(NOT WIN32)
add_subdirectory(otf)
#endif()
