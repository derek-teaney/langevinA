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

def get_qinput():
    qinput = {
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

        "init_amp" : 1, 
        "standing_waves" : False, 
        "init_dim" : 1
    }
    return qinput

# Definitely worth doing
def checkinputs(qinput):
    if qinput["chi"] != 1.:
        raise SystemExit('Chi should be one')
    if qinput["evolverType"] != "PV2HBSplit23":
        raise SystemExit('The evovlerType is not "PV2HBSplit23"')


# Dump the qinput into the .json file using the current value of qinput.
# This is the json file that is ultimately passed on the command line 
# Through -input tag.json
def qinputtojson(qinput):
    with open(qinput["outputfiletag"] + '.json', 'w') as outfile:
        json.dump(qinput, outfile, indent=4)

# Returns a default name for the run based on the value of qinput and the tag
# Usage 
#
# qin["outputfilename"] = get_qfilename("foo", qin)
#
# This sets outputiletag to foo for example: 
def get_qfilename(qinput, tag):
    name = "%s_N%05d_h%06d" % (tag, qinput["NX"], round(1000000*qinput["H"]))
    return name


# Findst the program looking in the environment variable for the path
def find_program(program_name) :
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
def prlmrun(program, qinput, time=2, append_nodeid=True, debug=False, dry_run=True, moreopts=[], seed=None, nnodes=1):

    # Create a run directory with name "qinput["outputfiletag"]"  if does not
    # exist, and cd to it to do the run
    dstack.pushd(qinput["outputfiletag"], mkdir=True)

    # Append a random 8 digit hex number to the tag labelling the run. This is
    # so that independent runs using the same inputfile, with different seeds,
    # can be run in the same directory. We have to modify the "qinput" structure
    if append_nodeid:
        oldtag = qinput["outputfiletag"]
        runid = "ffffffff"
        if not dry_run:
            runid = str(uuid.uuid4())[:8]
        qinput["outputfiletag"] = qinput["outputfiletag"] + "_{}".format(runid)

    # This is the new tag for the writing below
    tag = qinput["outputfiletag"]

    # Set the seed and write the inputfile to tag.json
    if seed is None:
        qinput["seed"] = random.randint(1, 2000000000)
    else:
        qinput["seed"] = seed

    # Check the current setup of inputs
    checkinputs(qinput)

    # Write the current qinput structure to the json file which will
    # be read by the program
    qinputtojson(qinput)

    #
    # Prepare the shell script for running on perlmutter
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
          (qinput["outputfiletag"]), file=fh)
    # Write the command that actually runds the program
    print("srun -n %d --cpu_bind=cores -c %d %s -input %s " %
          (tasks, cpuspertask, program, qinput["outputfiletag"]+'.json'), end=' ', file=fh)

    # This additional options are  added to the srun command
    print("-log_view", end=' ', file=fh)
    for opt in moreopts:
        print(opt, end=' ', file=fh)
    print(file=fh)

    # # Do any post processing of the run
    print('date  "+%%x %%T" >> %s_time.out' %
          (qinput["outputfiletag"]), file=fh)
    fh.close()

    # Submit the shell script to the batch system, executing
    #  
    # sbatch tag.sh
    if not dry_run:
        subprocess.run(['sbatch', filenamesh])

    # There was a side effect that the outputfiletag got modified
    # This should be undone for transparency
    if append_nodeid:
        qinput["outputfiletag"] = oldtag
    # return to the root directory
    dstack.popd()


if __name__ == "__main__":
    q1 = get_qinput() 
    q1["outputfiletag"] = get_qfilename(q1, "foo")
    print(q1)
