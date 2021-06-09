include_directories(${CMAKE_CURRENT_SOURCE_DIR}
                    "${CMAKE_CURRENT_SOURCE_DIR}/../Buffers")

find_package(LAPACK REQUIRED)

add_library(cudasirecon SHARED cudaSirecon.cpp boostfs.cpp gpuFunctionsImpl.cu pyInterface.cpp)
add_executable(cudaSireconDriver cudaSireconDriver.cpp)

target_include_directories(cudaSireconDriver PUBLIC "${PROJECT_BINARY_DIR}")


cuda_add_cufft_to_target(cudasirecon)
cuda_add_cufft_to_target(cudaSireconDriver)

set_property(TARGET cudasirecon PROPERTY MSVC_RUNTIME_LIBRARY "")
set_target_properties(cudaSireconDriver
  PROPERTIES
  MSVC_RUNTIME_LIBRARY ""
  OUTPUT_NAME cudasirecon
)

message(STATUS "IMLIB: " ${IMLIB})
message(STATUS "IVELIB: " ${IVELIB})
message(STATUS "LEG_STDIO: " ${LEG_STDIO})
message(STATUS "LAPACK_FOUND: " ${LAPACK_FOUND})

target_link_libraries(
  cudasirecon
  Buffer
  ${TIFF_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${BOOST_LIBRARIES}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
  ${IMLIB}
  ${IVELIB})

add_dependencies(cudaSireconDriver cudasirecon Buffer)
target_link_libraries(cudaSireconDriver cudasirecon)


if(WIN32)
  if(${MSVC_VERSION} GREATER 1800)
    target_link_libraries(cudaSireconDriver ${LEG_STDIO})
  endif()
else()
  target_link_libraries(cudasirecon rt)
endif(WIN32)

# added for make install to work in conda
set(HEADERS cudaSirecon.h cudaSireconImpl.h gpuFunctions.h)

install(
  TARGETS cudaSireconDriver cudasirecon
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

install(FILES ${HEADERS} DESTINATION include)
