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

# # FFTW2
# include_directories("${CMAKE_SOURCE_DIR}/fftw2/${PLATFORM}/include")
# find_library(FFTW sfftw
#              PATHS "${CMAKE_SOURCE_DIR}/fftw2/${PLATFORM}/lib" REQUIRED)
# find_library(RFFTW srfftw
#              PATHS "${CMAKE_SOURCE_DIR}/fftw2/${PLATFORM}/lib" REQUIRED)

# FFTW3
find_package(FFTW3 COMPONENTS single threads)

include_directories(${FFTW3_INCLUDE_DIR})

# add_executable(makeotf makeotf.c)
add_executable(makeotf radialft.cpp)

target_link_libraries(makeotf
                      ${IMLIB}
                      ${IVELIB}
					  ${FFTW3F_LIBRARY}
					  ${FFTW3F_THREADS_LIBRARY}
					  ${TIFF_LIBRARIES}
					  X11
					  )
                      # ${FFTW}
                      # ${RFFTW})

# if(UNIX AND NOT APPLE)
#   target_link_libraries(makeotf m)
# endif()

install(TARGETS makeotf RUNTIME DESTINATION bin)
