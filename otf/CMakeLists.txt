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
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/findFFTW")
if (WIN32)
  find_package(FFTW REQUIRED COMPONENTS FLOAT_LIB)
  include_directories( ${TIFF_INCLUDE_DIR} )
  include_directories( ${FFTW_INCLUDE_DIR} )
  message(STATUS "FFTW3 library" ${FFTW_FLOAT_LIB})
else ()
  find_package(FFTW3 REQUIRED COMPONENTS single threads)
endif()


if(WIN32)
endif(WIN32)

add_executable(makeotf radialft.cpp)
add_executable(otfviewer otfviewer.cpp)

if (WIN32)
  if ( ${MSVC_VERSION} GREATER 1800)
  target_link_libraries(makeotf
      ${IMLIB}
      ${IVELIB}
      ${FFTW_FLOAT_LIB}
      ${TIFF_LIBRARIES}
      ${LEG_STDIO}
  )
  target_link_libraries(otfviewer
      ${IMLIB}
      ${IVELIB}
      ${TIFF_LIBRARIES}
      ${LEG_STDIO}
  )
  else()
  target_link_libraries(makeotf
      ${IMLIB}
      ${IVELIB}
      ${FFTW_FLOAT_LIB}
      ${TIFF_LIBRARIES}
  )
  target_link_libraries(otfviewer
    ${IMLIB}
    ${IVELIB}
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
  target_link_libraries(otfviewer
    ${IMLIB}
    ${IVELIB}
	${TIFF_LIBRARIES}
	X11
	)
endif(WIN32)

install(TARGETS makeotf RUNTIME DESTINATION bin)
install(TARGETS otfviewer RUNTIME DESTINATION bin)
