# Create your small core library
add_library(simulation_core STATIC
    measurer.cxx
    measurer_output.cxx
    Stepper.cxx
    NoiseGenerator.cxx
    initialize.cxx
)

# Attach all headers and external libraries to the core library
# PUBLIC means anything that links to simulation_core also gets these dependencies
target_include_directories(simulation_core PUBLIC
    ${CMAKE_SOURCE_DIR}/src # Your local .h files
    ${PETSC_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
    ${FFTW_INCLUDE_DIRS}
)

target_link_options(simulation_core PUBLIC
    ${PETSC_LDFLAGS}
    ${HDF5_LDFLAGS}
    ${FFTW_LDFLAGS}
)

target_link_libraries(simulation_core PUBLIC
    ${PETSC_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${FFTW_LIBRARIES}
)


# Build SuperPions and link it to the core library
add_executable(SuperPions.exe SuperPions.cxx)
target_link_libraries(SuperPions.exe PRIVATE simulation_core)
