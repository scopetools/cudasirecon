#######################################################################
# Build instruction for CUDA SIMrecon project
#
# 1. Use conda to set up dependencies:
# conda env create -f environment.yml  # for linux
# conda env create -f environment_win.yml  # for windows
# conda activate simbuild
#
# 2. CUDA:
# The linux environment comes with the cudatoolkit-dev package that includes the nvcc compiler
# For windows, you will need to install the CUDA toolkit manually
#
# 3. MRC/DV support (optional)
# If you want to build with MRC support add the IVE libraries to the src directory
# with the following paths:
#   For Linux:
#     src/IVE/linux64/INCLUDE
#     src/IVE/linux64/LIB
#   For Windows:
#     src/IVE/win64/INCLUDE
#     src/IVE/win64/LIB
#
# If you DON'T want to build with MRC or do not have access to the IVE libraries
# please run cmake with the `-DBUILD_MRC=OFF` option
#
# 4. run build.sh (linux) or bld.bat (windows)
# inspect those scripts for more details
#
# 5. If building is successful, an executable cudaSireconDriver is generated in
#  cmake_build/cudaSirecon. Run the test to make sure it works:
#  $ test_data/run.sh
#
#######################################################################

option(BUILD_MRC "Build cudasirecon with MRC support" ON)

if (BUILD_MRC)
  message(STATUS "** Building WITH mrc support **")
  cmake_minimum_required (VERSION 3.14)
else()
  message(STATUS "** Building WITHOUT mrc support **")
  cmake_minimum_required (VERSION 3.19)
endif()

project (cudaSIMRecon CXX CUDA)

# set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

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

add_definitions(-Dcimg_display=0)

find_package(OpenMP)
message(STATUS "OpenMP FOUND: " ${OPENMP_FOUND})
if(OPENMP_FOUND)
  # nvcc and openmp with CImg aren't working at the moment
  # https://github.com/dtschump/CImg/issues/286
  add_definitions(-Dcimg_use_openmp=0)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


find_package(CUDA)
find_package(TIFF)
find_package(Boost COMPONENTS program_options filesystem system REQUIRED)
if(Boost_FOUND)
  link_directories ( ${Boost_LIBRARY_DIRS} )
  include_directories(${Boost_INCLUDE_DIRS})
endif()


if(WIN32)
  if ( ${MSVC_VERSION} GREATER 1800)
    find_library(LEG_STDIO legacy_stdio_definitions)
    message(STATUS "adding legacy stdio definitions for >VS14.0 " ${LEG_STDIO})
  endif()
endif(WIN32)


set(CUDA_NVCC_FLAGS
  -D_USE_MATH_DEFINES
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

if(${CUDA_VERSION_MAJOR} GREATER 9)
  set(CUDA_NVCC_FLAGS
    ${CUDA_NVCC_FLAGS};
    -gencode=arch=compute_75,code=sm_75;
    -gencode=arch=compute_75,code=compute_75
  )
endif()
message(STATUS "CUDA_NVCC_FLAGS: " ${CUDA_NVCC_FLAGS})

# # Set permissions before adding subdirectories
# set(SCI_GROUP_WRITE GROUP_WRITE)
# set(SCI_WORLD_FILE_PERMS WORLD_READ)
# set(SCI_WORLD_PROGRAM_PERMS WORLD_READ WORLD_EXECUTE)


add_subdirectory(Buffers)
add_subdirectory(cudaSirecon)
add_subdirectory(otf)

# Building the documentation
# option(BUILD_DOCUMENTATION "Use Doxygen to create the HTML based API documentation" OFF)
# if(BUILD_DOCUMENTATION)
#   FIND_PACKAGE(Doxygen)
#   if (NOT DOXYGEN_FOUND)
#     message(FATAL_ERROR 
#       "Doxygen is needed to build the documentation. Please
#       install it correctly")
#     endif()
#     #-- Configure the Template Doxyfile for our specific project
#     # configure_file(Doxyfile.in
#     #   ${PROJECT_BINARY_DIR}/Doxyfile  @ONLY IMMEDIATE)
#     #-- Add a custom target to run Doxygen when ever the project is built
#     add_custom_target (Docs ALL 
#       COMMAND ${DOXYGEN_EXECUTABLE}
#       ${PROJECT_BINARY_DIR}/Doxyfile
#       SOURCES
#       ${PROJECT_BINARY_DIR}/Doxyfile)
# endif()