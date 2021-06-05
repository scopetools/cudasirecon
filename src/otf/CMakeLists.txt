if(WIN32)
  set(PLATFORM win64)
else()
  set(PLATFORM linux64)
endif(WIN32)

set(PRIISM_LIB_PATH "${CMAKE_SOURCE_DIR}/IVE/${PLATFORM}/LIB")
set(PRIISM_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/IVE/${PLATFORM}/INCLUDE")

include_directories(${PRIISM_INCLUDE_PATH})
include_directories("${CMAKE_SOURCE_DIR}/cudaSirecon")  # for CImg.h
link_directories(${PRIISM_LIB_PATH})

find_library(IMLIB NAMES imlib libimlib PATHS ${PRIISM_LIB_PATH})
find_library(IVELIB NAMES ive libive PATHS ${PRIISM_LIB_PATH})


if (WIN32)
  set (FFTW3_INCLUDE_DIR $ENV{CONDA_PREFIX}/Library/lib)
  set (FFTW3F_LIBRARY fftw3f)
  find_package(FFTW3 COMPONENTS single threads REQUIRED)
else()
  set (FFTW3_INCLUDE_DIR $ENV{CONDA_PREFIX}/include)
  find_library(FFTW3F_LIBRARY fftw3f REQUIRED)
  find_library(FFTW3F_THREADS_LIBRARY fftw3f_threads REQUIRED)
endif()

add_executable(makeotf radialft.cpp)
target_link_libraries(makeotf ${FFTW3F_LIBRARY} ${TIFF_LIBRARIES})

if (IMLIB)
  target_link_libraries(makeotf ${IMLIB} ${IVELIB})
endif ()


if (WIN32)
  if (${MSVC_VERSION} GREATER 1800)
  	target_link_libraries(makeotf ${LEG_STDIO})
  endif()
else ()
  target_link_libraries(makeotf ${FFTW3F_THREADS_LIBRARY} X11)
endif ()

install(TARGETS makeotf RUNTIME DESTINATION bin)
