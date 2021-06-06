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


add_executable(makeotf radialft.cpp)
if (IMLIB)
  add_definitions(-DMRC)
  target_link_libraries(makeotf ${IMLIB} ${IVELIB})
endif ()


find_library (FFTW3F_LIBRARY fftw3f REQUIRED)
if (WIN32)
  find_package(FFTW3 COMPONENTS single threads REQUIRED)
  if (${MSVC_VERSION} GREATER 1800)
    target_link_libraries(makeotf ${LEG_STDIO})
  endif()
else()
  find_library(FFTW3F_THREADS_LIBRARY fftw3f_threads REQUIRED)
  target_link_libraries(makeotf ${FFTW3F_THREADS_LIBRARY} X11)
endif()

message(STATUS "FFTW3F_LIBRARY " ${FFTW3F_LIBRARY})
target_link_libraries(makeotf ${FFTW3F_LIBRARY} ${TIFF_LIBRARIES})


install(TARGETS makeotf RUNTIME DESTINATION bin)
