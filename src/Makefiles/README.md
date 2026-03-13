## Compilation on Perlmutter ##

There is a script which sets up the environment called

setupprlm.sh 

Just do 

```sh
source setupprlm.sh
```

Then the  the code should compile with 

make -f Makefile.prlm 

The script consists of 
```sh
module load spack
spack env activate gcc
spack load petsc
module load cray-fftw
```

The ingrediants are:
-Spack This is a scientific package manager https://docs.nersc.gov/development/build-tools/spack/
-Switch the spack environment to gcc
-Load petsc
-Load the nersc fftw



