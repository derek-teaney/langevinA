# my_project/CMakeLists.txt
cmake_minimum_required(VERSION 3.15)
project(ModelG VERSION 1.0 LANGUAGES CXX)

# Talk directly to the OS to get the compilers
execute_process(
    COMMAND pkg-config --variable=ccompiler PETSc
    OUTPUT_VARIABLE PETSC_CC
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)
execute_process(
    COMMAND pkg-config --variable=cxxcompiler PETSc
    OUTPUT_VARIABLE PETSC_CXX
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)

# If the shell successfully returned the compilers, force CMake to use them
if(PETSC_CC AND PETSC_CXX)
    set(CMAKE_C_COMPILER ${PETSC_CC})
    set(CMAKE_CXX_COMPILER ${PETSC_CXX})
    message(STATUS "Shell override: Set C Compiler to ${CMAKE_C_COMPILER}")
    message(STATUS "Shell override: Set CXX Compiler to ${CMAKE_CXX_COMPILER}")
endif()

# Enforce C++17 and optimization
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()


# Find all external packages via pkg-config
find_package(PkgConfig REQUIRED)
pkg_check_modules(PETSC REQUIRED PETSc)
pkg_check_modules(HDF5 REQUIRED hdf5)
pkg_check_modules(FFTW REQUIRED fftw3)

# Generate gitversion.h into the build directory
execute_process(
    COMMAND git rev-parse HEAD
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
set(GIT_VERSION_FILE "${CMAKE_CURRENT_BINARY_DIR}/gitversion.h")
file(WRITE ${GIT_VERSION_FILE}
    "#ifndef GITVERSION_H\n#define GITVERSION_H\nconst char *gitversion = \"${GIT_HASH}\";\n#endif\n"
)

# Ensure our targets can find the generated gitversion.h
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Hand off the actual building to the src/ directory
add_subdirectory(src)
