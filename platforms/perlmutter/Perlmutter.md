# Script
========

From the ModelG directory one executes
```bash
ModelG> source ./Makefiles/setupprlm.sh
ModelG> make -f ./Makefiles/Makefile.prlm
```
You can look at Makefiles/setupprlm.sh to see what it does. 

# How it works
=============

```bash
# The system wide spack instance (installed with git is in $SPACK_ROOT)
module load spack
# This sets the system wide defaults which are in $SPACK_SYSTEM_CONFIG_PATH
module load spack-config
module load cray-fftw
# Activate the spack environment in our common file system area
spack env activate $CFS/m3722/opt/prlm/petsc-cpu-int64
# print out the envirnoment to make sure we know what we are working with
spack env status
spack load petsc
# From the ModelG directory issue the build command
make -f Makefiles/Makefile.prlm
```

# Installing Petsc In A Home Directory Spack package manager 
==========================================================

Suppose we want our own Petsc with whatever

```bash
module load spack
module load spack-config
```
The second one sets up  the spack compilers at nersc. 
Now we can do 

```bash
spack env create mypetsc
spack env activate mypetsc
# Print out the environment to make sure we know what we are working with
spack env status
# Add the software we want
spack add petsc +int64 +hdf5
spack install
```
Ok so we created the petsc-int64 environmnet with the software we need

Then, imagine we create a new login shell.  The process of using mypetsc environment  is

```bash
module load spack
module load spack-config
module load cray-fftw
spack env activate mypetsc
spack load petsc
```

With these steps we are ready to use the petsc universe

# Petsc Common Environment On Perlmutter
========================================

Derek set up PETSc to work in  $CFS/m3722/opt/ so everyone doesn't need to spin up there own version of PETSc. Though this wasn't too bad

Here we first created the shared spack environment in the shared directory. Note the '-d'

```bash
spack create -d ${shared_directory}/petsc-int64
spack acivate ${shared_directory}/petsc-int64
spack add petsc+hdf5+int64
```
Before installation, in the environment directory we editted the "spack.yaml" file

```yaml
spack:
  # add package specs to the `specs` list
  specs:
  - petsc+hdf5+int64
  config:
    install_tree:
      root: $env/opt
  view: true
  concretizer:
    unify: true
```
We added the lines
```yaml
  config:
    # This tells Spack to install packages into a subfolder named 'opt'
    # right inside your project directory.
    # $env is the director of the environment
    install_tree:
      root: $env/opt
```
This guarantees that the software is built in the common directory and not in my home directory.

Then we install in the usual way
```bash
# This may not be needed, and just makes sure that after editting spack.yaml everything is ok
spack concretize -f
spack install
```

# Spack View
============

The directory where the binaries are stored is 

```bash
Updating view at /global/cfs/cdirs/m3722/opt/prlm/petsc-cpu-int64/.spack-env/view
```

Since petsc is installing its own hdf5 we find 

```bash
(base) dteaney@perlmutter:login35:~> which h5dump
/global/cfs/cdirs/m3722/opt/prlm/petsc-cpu-int64/.spack-env/view/bin/h5dump
```

Similarly we find
```bash
pkg-config --variable=pcfiledir PETSc
/global/cfs/cdirs/m3722/opt/prlm/petsc-cpu-int64/.spack-env/view/lib/pkgconfig
```
One can look at the file and it is kind of helpful, e.g. 
```bash
(base) dteaney@perlmutter:login35:~> cat `pkg-config --variable=pcfiledir PETSc`/PETSc.pc
prefix=/global/cfs/cdirs/m3722/opt/prlm/petsc-cpu-int64/opt/linux-zen3/petsc-3.24.1-zjjlfkajgqfrhgmefahoq7trcjuwecuy
ccompiler=/opt/cray/pe/mpich/8.1.30/ofi/gnu/12.3/bin/mpicc
# etc, etc 
```
