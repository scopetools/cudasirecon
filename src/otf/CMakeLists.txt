set(CMAKE_MACOSX_RPATH 1)

if(WIN32)
  set(PLATFORM win64)
elseif(APPLE)
  set(PLATFORM darwin64)
else()
  set(PLATFORM linux64)
endif(WIN32)

set(PRIISM_LIB_PATH "${CMAKE_SOURCE_DIR}/IVE/${PLATFORM}/LIB")
set(PRIISM_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/IVE/${PLATFORM}/INCLUDE")

include_directories(${PRIISM_INCLUDE_PATH})
include_directories("${CMAKE_SOURCE_DIR}/cudaSirecon")  # for CImg.h
link_directories(${PRIISM_LIB_PATH})
find_library(IMLIB imlib PATHS ${PRIISM_LIB_PATH} REQUIRED)
find_library(IVELIB ive PATHS ${PRIISM_LIB_PATH} REQUIRED)

# FFTW3
# set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/findFFTW")
# if (WIN32)
#   find_package(FFTW3 REQUIRED COMPONENTS single)
#   include_directories( ${TIFF_INCLUDE_DIR} )
#   include_directories( ${FFTW3_INCLUDE_DIR} )

# else ()
#   find_package(FFTW3 REQUIRED COMPONENTS single threads)
# endif()

if (WIN32)
  set (FFTW3_INCLUDE_DIR $ENV{CONDA_PREFIX}/Library/lib)
  set (FFTW3F_LIBRARY fftw3f)
  find_package(FFTW3 COMPONENTS single threads REQUIRED)
else()
  set (FFTW3_INCLUDE_DIR $ENV{CONDA_PREFIX}/include)
  find_library(FFTW3F_LIBRARY fftw3f REQUIRED)
  find_library(FFTW3F_THREADS_LIBRARY fftw3f_threads REQUIRED)
endif()



# add_executable(makeotf makeotf.c)
add_executable(makeotf radialft.cpp)

if (WIN32)
  if ( ${MSVC_VERSION} GREATER 1800)
	target_link_libraries(makeotf
      ${IMLIB}
      ${IVELIB}
	  ${FFTW3F_LIBRARY}
	  ${TIFF_LIBRARIES}
      ${LEG_STDIO}
	  )
  else()
	target_link_libraries(makeotf
      ${IMLIB}
      ${IVELIB}
	  ${FFTW3F_LIBRARY}
	  ${TIFF_LIBRARIES}
	  )
  endif()
else ()
  target_link_libraries(makeotf
    ${IMLIB}
    ${IVELIB}
	${FFTW3F_LIBRARY}
	${FFTW3F_THREADS_LIBRARY}
	${TIFF_LIBRARIES}
	X11
	)
endif(WIN32)

install(TARGETS makeotf RUNTIME DESTINATION bin)
