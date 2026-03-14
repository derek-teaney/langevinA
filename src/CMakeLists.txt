# Create your small core library
add_library(modelg
    measurer.cxx
    measurer_output.cxx
    Stepper.cxx
    NoiseGenerator.cxx
    initialize.cxx
)

# Attach all headers and external libraries to the core library
# PUBLIC means anything that links to modelg also gets these dependencies
target_include_directories(modelg PUBLIC
    ${CMAKE_SOURCE_DIR}/src # Your local .h files
    ${PETSC_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
    ${FFTW_INCLUDE_DIRS}
)

target_link_options(modelg PUBLIC
    ${PETSC_LDFLAGS}
    ${HDF5_LDFLAGS}
    ${FFTW_LDFLAGS}
)

target_link_libraries(modelg PUBLIC
    ${PETSC_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${FFTW_LIBRARIES}
)


# Build SuperPions and link it to the core library
add_executable(SuperPions.exe SuperPions.cxx)
target_link_libraries(SuperPions.exe PRIVATE modelg)
