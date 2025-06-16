#!/usr/bin/env python3
import os
import glob
import subprocess
import json
import random
import dstack
import datetime
import uuid
import math


random.seed()

data = {
    # lattice dimension
    "NX": 32,

    # Time stepping
    "finaltime": 10,
    "initialtime": 0,
    "deltat": 0.24,

    # Action
    "mass0": -4.70052,
    "dmassdt": 0,

    "lambda": 4.,
    "H": 0.003,
    "chi": 5.,
    "gamma": 1.,
    "diffusion": 0.3333333,

    # initial condition"
    "evolverType": "PV2HBSplit23",
    "seed": 122335456,
    "restart": False,
    "outputfiletag": "grun",
    "saveFrequency": 3,
    "thermalization_time": 0.0,

    # for quenched initial conditions
    "quench_mode": False, 
    "quench_mode_mass0": -4.70052,

    # For running multi-events
    "eventmode": False,
    "nevents": 1,
    "diffusiononly": False,

    # parameters for superfluid events
    "superfluidmode": False, 
    "f2_constant" : 1, 
}

#
def checkinputs():
    if data["mass0"] > 0 :
        raise SystemExit('The parameters mass0 should be negative')
    if data["dmassdt"] > 0 :
        raise SystemExit('The parameters dmassdt should be negative')
    if data["diffusiononly"]:
        raise SystemExit('Do not run in diffusiononly mode without asking derek')
    if data["chi"] != 5.:
        raise SystemExit('Chi should be five')
    if data["evolverType"] != "PV2HBSplit23":
        raise SystemExit('The evovlerType is not "PV2HBSplit23"')


# dump the data into a .json file
def datatojson():
    with open(data["outputfiletag"] + '.json', 'w') as outfile:
        json.dump(data, outfile, indent=4)

# Canonicalize the names for a given set of parameters
def get_kzfilename(tag):
    name = "%s_N%03d_m%08d_h%06d_tkz%06d" % (tag, data["NX"], round(
        100000*data["mass0"]), round(1000000*data["H"]), round(1./data["dmassdt"]))
    return name

# Canonicalize the names for a given set of parameters
def get_qkzfilename(tag):
    name = "%s_N%03d_m%08d_h%06d_q" % (tag, data["NX"], round(
        100000*data["mass0"]), round(1000000*data["H"]))
    return name

def getdefault_filename(tag):
    tag = data["outputfiletag"]
    name = "%s_N%03d_m%08d_h%06d_c%05d" % (tag, data["NX"], round(
        100000*data["mass0"]), round(1000000*data["H"]), round(100*data["chi"]))
    return name


# Find the program looking in the environment variable for the path
def find_program(program_name="SuperPions.exe") :
    path = os.environ.get('MODELGPATH') 
    if path is None:
        print("Unable to find the path MODELGPATH") 
    abspath = os.path.join(path, program_name)
    if os.path.exists(abspath) :
        print("Found the executable {}".format(abspath))
    else:
        print("Unable to find the executable {}".format(abspath))
    return abspath


#########################################################################
# Runs on perlmutter
#########################################################################


def prlmrun(time=2, debug=False, dry_run=True, moreopts=["-log_view"], seed=None, nnodes=1, nodeid=False):
    prgm = find_program()

    # Create a run directory "name"  if does not exist, and cd to it
    dstack.pushd(data["outputfiletag"], mkdir=True)

    # If nodeid is True then append a random 8 digit hex number
    # to the tag labelling the run. This is so that independent runs using the
    # same inputfile, with different seeds, can be run in the same directory
    if nodeid: 
        oldtag = data["outputfiletag"]
        runid = "ffffffff"
        if not dry_run:
            runid = str(uuid.uuid4())[:8]
        data["outputfiletag"] = data["outputfiletag"] + "_{}".format(runid)

    tag = data["outputfiletag"]

    #
    checkinputs()

    # Set the seed and write the inputfile to tag.json
    if seed is None:
        data["seed"] = random.randint(1, 2000000000)
    else:
        data["seed"] = seed
    datatojson()

    #
    # Prepare the shell script
    #
    filenamesh = tag + '.sh'

    fh = open(filenamesh, 'w')

    tasks = int(nnodes*128)
    cpuspertask = int(2*128/(tasks/nnodes))
    print("#!/bin/bash", file=fh)
    if debug:
        print("#SBATCH -A m3722", file=fh)
        print("#SBATCH -C cpu", file=fh)
        print("#SBATCH --qos debug", file=fh)
        print("#SBATCH -t 00:30:00", file=fh)
        print("#SBATCH -N {}".format(nnodes), file=fh)
        print("#SBATCH --ntasks={}".format(tasks), file=fh)
        print("#SBATCH --cpus-per-task={}".format(cpuspertask), file=fh)
    else:
        print("#SBATCH -A m3722", file=fh)
        print("#SBATCH -C cpu", file=fh)
        print("#SBATCH -q regular", file=fh)
        print("#SBATCH -t {}".format(int(math.ceil(time*60.))), file=fh)
        print("#SBATCH -N {}".format(nnodes), file=fh)
        print("#SBATCH --ntasks={}".format(tasks), file=fh)
        print("#SBATCH --cpus-per-task={}".format(cpuspertask), file=fh)

    # Set up the shell environment
    print("", file=fh)
    print("export HDF5_DISABLE_VERSION_CHECK=2", file=fh)
    print("", file=fh)
    print("#run the application:", file=fh)
    print('date  "+%%x %%T" > %s_time.out' %
          (data["outputfiletag"]), file=fh)
    # Write the command that actually runds the program
    print("srun -n %d --cpu_bind=cores -c %d %s -input %s " %
          (tasks, cpuspertask, prgm, data["outputfiletag"]+'.json'), end=' ', file=fh)
    # This additional options are  added to the srun command
    for opt in moreopts:
        print(opt, end=' ', file=fh)
    print(file=fh)

    # # Do any post processing of the run
    # programpy = find_program(program_name="x2k.py") 
    # print("python {} {}.json".format(programpy,data["outputfiletag"]), file=fh)
    # print(file=fh)

    print('date  "+%%x %%T" >> %s_time.out' %
          (data["outputfiletag"]), file=fh)
    fh.close()

    # Submit the shell script
    if not dry_run:
        subprocess.run(['sbatch', filenamesh])

    # There was a side effect that the outputfiletag got modified
    # This should be undone for transparency
    if nodeid: 
        data["outputfiletag"] = oldtag
    # return to the root directory
    dstack.popd()


#########################################################################
# Runs on seawulf  with time in batch time. One should set dry_run=False to
# actually run the code
#########################################################################
GLOBAL_PETSCPKG_PATH_SEAWULF = "${PKG_CONFIG_PATH}:/gpfs/home/adrflorio/petsc/arch-linux2-c-debug/lib/pkgconfig/"

def seawulfrun(time="00:02:00", debug=False, shared=False, dry_run=True, moreopts=[]):
    nprocesses = 24
    filenamesh = data["outputfiletag"] + '.sh'
    with open(filenamesh, 'w') as fh:
        print("#!/bin/bash", file=fh)
        if debug:
            print("#SBATCH -p debug-{}core".format(nprocesses), file=fh)
            print("#SBATCH --time=00:10:00", file=fh)
            print("#SBATCH --nodes=1", file=fh)
            print("#SBATCH --ntasks-per-node={}".format(nprocesses), file=fh)
        else:
            print("#SBATCH -p long-{}core".format(nprocesses), file=fh)
            print("#SBATCH --time={}".format(time), file=fh)
            print("#SBATCH --nodes=1", file=fh)
            print("#SBATCH --ntasks-per-node={}".format(nprocesses), file=fh)

        print("", file=fh)
        print("module load shared", file=fh)
        print("module load gcc-stack", file=fh)
        print("module load hdf5/1.10.5-parallel", file=fh)
        print("module load fftw3", file=fh)
        print("module load cmake", file=fh)
        print("module load gsl", file=fh)
        print("export PKG_CONFIG_PATH={}".format(
            GLOBAL_PETSCPKG_PATH_SEAWULF), file=fh)
        print("export MV2_ENABLE_AFFINITY=0", file=fh)
        print("", file=fh)
        print("#run the application:", file=fh)

        print('date  "+%%x %%T" > %s_time.out' %
              (data["outputfiletag"]), file=fh)
        # get the program
        path = os.path.abspath(os.path.dirname(__file__))
        prgm = path + "/SuperPions.exe"
        # set the seed and the inputfile
        data["seed"] = random.randint(1, 2000000000)

        # Write the data to an .json
        datatojson()

        # write the command that actually runds the program
        basename = "./" + os.path.basename(data["outputfiletag"])
        print("mpirun -n {} {} -input {} ".format(nprocesses,
              prgm, basename + '.json'), end=' ', file=fh)
        for opt in moreopts:
            print(opt, end=' ', file=fh)
        print(file=fh)
        print('date  "+%%x %%T" >> %s_time.out' %
              (data["outputfiletag"]), file=fh)

    if not dry_run:
        subprocess.run(['sbatch', filenamesh])

# runs the actual command current value of data  with mpiexec

########################################################################
# runs the program with current value of data  and mpiexec on local
# mac.
########################################################################
def run(program_name="SuperPions.exe", moreopts=[], dry_run=True, time=0, seed=None, ncpus="2", log_view=True, mpiexec="mpiexec"):

    prgm = find_program(program_name)
    tag = data["outputfiletag"]

    # Go to the directory
    dstack.pushd(tag, mkdir=True)

    # set the seed and the inputfile
    if seed is None:
        data["seed"] = random.randint(1, 2000000000)
    else:
        data["seed"] = seed

    datatojson()

    # Execute the program
    opts = [mpiexec, "-n", ncpus, prgm,
            "-input",  tag + '.json']
    if log_view:
        opts.append('-log_view')
    opts.extend(moreopts)
    print(opts)
    if not dry_run:
        subprocess.run(opts)

    # Go back to the working directory
    dstack.popd()


if __name__ == "__main__":
    print(getdefault_filename("foo"))
