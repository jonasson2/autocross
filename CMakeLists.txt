cmake_minimum_required(VERSION 3.10)

# Use gfortran as the compiler
set(CMAKE_Fortran_COMPILER gfortran)

# Project name and Fortran as the language
project(FortranPrograms Fortran)

# Determine the shared flag based on platform
if(WIN32)
    set(SHARED_FLAG "-D_WINDOWS")
elseif(NOT APPLE)
    set(SHARED_FLAG "-fPIC")
else()
    set(SHARED_FLAG "")  # No shared flag needed for macOS
endif()

# Compiler flags for gfortran
set(CMAKE_Fortran_FLAGS "-O2 -g ${SHARED_FLAG}")
set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g ${SHARED_FLAG}")
set(CMAKE_Fortran_FLAGS_RELEASE "-O2 ${SHARED_FLAG}")
set(CMAKE_Fortran_EXTENSIONS OFF)

# Set output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/modules)

# Set the Fortran module directory for object files
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_SOURCE_DIR}/modules)

# Source files
set(COMMON_SOURCES src/common_modules.f90)
set(PEARSONT3_SOURCES src/p3_modules.f90 src/p3_subroutine.f90 ${COMMON_SOURCES})
set(REDFITX_SOURCES src/rx_modules.f90 src/rx_subroutines.f90 ${COMMON_SOURCES})

# Main program files
set(PEARSONT3_MAIN src/pearsont3.f90)
set(REDFITX_MAIN src/redfitx.f90)

# Determine suffixes based on the operating system
if(WIN32)
    set(EXE_SUFFIX ".exe")
    set(SHARED_SUFFIX ".dll")
else()
    set(EXE_SUFFIX "")
    set(SHARED_SUFFIX ".so")
endif()

# Create the pearsont3 executable
add_executable(pearsont3 ${PEARSONT3_MAIN} ${PEARSONT3_SOURCES})
set_target_properties(pearsont3 PROPERTIES
    OUTPUT_NAME pearsont3
    SUFFIX "${EXE_SUFFIX}"
)

# Create the redfitx executable
add_executable(redfitx ${REDFITX_MAIN} ${REDFITX_SOURCES})
set_target_properties(redfitx PROPERTIES
  OUTPUT_NAME redfitx
  SUFFIX "${EXE_SUFFIX}"
)

# Create the redfitx dynamic library
add_library(redfitx_shared SHARED ${REDFITX_SOURCES})
set_target_properties(redfitx_shared PROPERTIES
  OUTPUT_NAME redfitx
  SUFFIX "${SHARED_SUFFIX}"
)

target_link_libraries(redfitx PRIVATE redfitx_shared)

# Create the pearsont3 dynamic library
add_library(pearsont3_shared SHARED ${PEARSONT3_SOURCES})
set_target_properties(pearsont3_shared PROPERTIES
    OUTPUT_NAME pearsont3
    SUFFIX "${SHARED_SUFFIX}"
)

# Ensure all executables and libraries are linked with the Fortran runtime
target_link_libraries(pearsont3 pearsont3_shared)
