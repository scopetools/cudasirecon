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
link_directories(${PRIISM_LIB_PATH})
find_library(IMLIB imlib PATHS ${PRIISM_LIB_PATH} REQUIRED)
find_library(IVELIB ive PATHS ${PRIISM_LIB_PATH} REQUIRED)

# FFTW2
include_directories("${CMAKE_SOURCE_DIR}/fftw2/${PLATFORM}/include")
find_library(FFTW sfftw
             PATHS "${CMAKE_SOURCE_DIR}/fftw2/${PLATFORM}/lib" REQUIRED)
find_library(RFFTW srfftw
             PATHS "${CMAKE_SOURCE_DIR}/fftw2/${PLATFORM}/lib" REQUIRED)

add_executable(makeotf makeotf.c)

target_link_libraries(makeotf
                      ${IMLIB}
                      ${IVELIB}
                      ${FFTW}
                      ${RFFTW})

if(UNIX AND NOT APPLE)
  target_link_libraries(makeotf m)
endif()

install(TARGETS makeotf RUNTIME DESTINATION bin)
