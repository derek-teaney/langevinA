# SuperPions Simulation

This repository contains the C/C++ source code for the Model G Simulation

## Dependencies

To build and run this simulation, you will need the following installed on your system:
* **CMake** 
* **PETSc**
* **HDF5**
* **FFTW3**

## Standard Build Instructions (Local Machine)

If you are building this project on a standard local machine or laptop with the dependencies installed, you can use the standard CMake workflow:

```bash
# 1. Create a build directory
mkdir build
cd build

# 2. Configure the project
cmake -DCMAKE_BUILD_TYPE=Release ..

# 3. Compile the code
make
```

## Building on Perlmutter

You need to setup the perlmuttter environment
```bash
source platforms/perlmutter/setupprlm.sh
```
More detailed instructions are in Perlmutter.md
