cmake_minimum_required(VERSION 3.0)

project(dadafits C)

set(CMAKE_C_FLAGS_RELEASE "-O3 -march=native") 

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

find_package (psrdada REQUIRED)
find_package (cfitsio REQUIRED)
find_package (CUDA REQUIRED)

# expose some variables to the source code
set (dadafits_VERSION_MAJOR 0)
set (dadafits_VERSION_MINOR 1)
configure_file ("src/config.h.in" "${PROJECT_BINARY_DIR}/config.h")
include_directories ("${PROJECT_BINARY_DIR}")

include_directories ("${CFITSIO_INCLUDE_DIR}")

add_executable(dadafits
    src/main.c
    src/downsample.c
    src/sb_util.c
    src/fits_io.c
    src/manipulate.c
    src/dadafits_internal.h
)
add_executable(fits_dump
    src/fits_dump.c
)
target_link_libraries(dadafits ${PSRDADA_LIBRARIES} ${CUDA_LIBRARIES} ${CFITSIO_LIBRARIES} -lm)
target_link_libraries(fits_dump ${PSRDADA_LIBRARIES} ${CUDA_LIBRARIES} ${CFITSIO_LIBRARIES})

install(TARGETS dadafits RUNTIME DESTINATION bin)
install(TARGETS fits_dump RUNTIME DESTINATION bin)

