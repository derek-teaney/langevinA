Project
=======

* My project is m3722
* UserId=dteaney(69248) GroupId=dteaney(69248)

Login Notes
===========

* Execute the script  ./cori

ssh -i ~/.ssh/nersc dteaney@cori.nersc.gov

Enter password+OTP i.e. no spaces

* Once per day need to execute ./sshproxy.sh -u dteaney

* Then you can login via ssh 

File System
===========

* The home directory has /global/homes/d/dteaney/

* The CFS (Common File System) is $CFS/m3722

/global/cfs/cdirs/m3722> 

  This is for all of our group

* The scratch directory is for working I/O for running jobs. This is $SCRATCH = dteaney@cori07:/global/cscratch1/sd/dteaney

PETSc
=====

* Don't use pkg-config.  They have done all this 

   module load cray-petsc

* After calling this command one compiles a PETSc program with the
compiler wrapper

   CC mypetsc.cxx

GSL
===

* Use the module load gsl

   modulue load gsl

* Use the environment variable to compile

   CC mygsl.cxx $GSL

* module show gsl

Compiling the o4model
=====================

* First type

module load gsl
module load petsc


Starting a debug session
========================

* To start a debug session we allocate a session 


Node
====

* One server, this should normally be 1. Then the number of cores.  We normally
take one node with many tasks.

Batch System
=============

* We submit a job with sbatch  myjob.sh. The first lines of this script is

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=8:20:00
#SBATCH --constraint=haswell
#SBATCH --qos=regular
#SBATCH --account=m3722

* The main command is like "mpiexec" is 

srun -n 32  /global/u2/d/dteaney/o4/src/o4model.exe 

*  You can monintor the running with 

squeue -u dteaney

scontrol show job  39487350

*  Don't particularly understand -c option for srun

*  Which is "better for me"  the haswell or KNL systems.

*  Why does the output of 

   NumNodes=1 NumCPUs=64 NumTasks=1 CPUs/Task=1 ReqB:S:C:T=0:0:*:*
   TRES=cpu=64,mem=118G,node=1,billing=64
   Socks/Node=* NtasksPerN:B:S:C=0:0:*:* CoreSpec=*
   MinCPUsNode=1 MinMemoryNode=118G MinTmpDiskNode=0
   Features=haswell DelayBoot=2-00:00:00

