# Load the e4s stack
====================

module load spack

# Switch the spack stack
====================
spack env activate gcc

spack load petsc

# FFTW
# We need to load fftw
====================

module load cray-fftw
