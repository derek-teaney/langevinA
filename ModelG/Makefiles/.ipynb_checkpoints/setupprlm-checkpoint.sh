module load spack
module load spack-config
PETSC_SPACK_ENV=$CFS/m3722/opt/prlm/petsc-cpu-int64
spack env activate $PETSC_SPACK_ENV
spack load petsc
module load cray-fftw
