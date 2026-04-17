#!/usr/bin/env python3
import json
import math
import os
import random
import subprocess
import sys
import uuid

import dstack

random.seed()

data = {
    # lattice dimension
    "NX": 32,

    # Time stepping
    "finaltime": 10,
    "initialtime": 0,
    "deltat": 0.24,

    # Action of O4 mdoel parameters
    "mass0": -4.70052,
    "dmassdt": 0,
    "lambda": 4.0,
    "H": 0.003,
    "chi": 5.0,

    # Transport Coefficients
    "gamma": 1.0,
    "diffusion": 0.3333333,
    # evolverType: selects the time-evolution algorithm. Options:
    #
    #   "PV2HBSplit23"       -- Default. Predictor-Verlet + heat-bath split
    #                           stepper, 2nd/3rd order, fixed sequence
    #                           "ABBABBABBC". A = ideal (Verlet) substep,
    #                           B = heat-bath substep, C = diffusion substep.
    #
    #   "PV2HBSplitGeneral"  -- Same stepper but with a user-defined step
    #                           sequence and toggles for each substep type,
    #                           controlled by the "pv2hb_split_general" block.
    #
    #   "SuperSplitStep"     -- Split stepper for the superfluid sector.
    #                           Use with superfluidmode=True. Requires a
    #                           "SuperSplitStep" block with keys
    #                           "step_sequence" (str, default "AB") and
    #                           "use_implicit_step" (bool, default false).
    "evolverType": "PV2HBSplit23",
    # Full control over the stepper. This is when evolverType is set to
    # "PV2HBSplitGeneral"
    "pv2hb_split_general": {
        "steps": "ABBABBABBC",
        "include_ideal": True,
        "include_heatbath": True,
        "include_diffusion": True,
    },
    "seed": 122335456,
    "restart": False,
    "outputfiletag": "grun",
    "saveFrequency": 3,
    "thermalization_time": 0.0,
    # Options are [ "default", "restart", "quench_mode" ,  "randomspins"]
    "initialization": "default",
    # for quenched initial conditions
    "quench_mode": False,
    "quench_mode_mass0": -4.70052,
    # For running multi-events
    "eventmode": False,
    "nevents": 1,
    # parameters for superfluid events
    "superfluidmode": False,
    "f2_constant": 5.0,
}

# This is a flag to toggle input checking
CheckInputs = True

# Check that the setup is sensible


def checkinputs():
    if not CheckInputs:
        return
    if data["mass0"] > 0:
        raise SystemExit("The parameters mass0 should be negative")
    if data["dmassdt"] > 0:
        raise SystemExit("The parameters dmassdt should be negative")
    if data["chi"] != 5.0:
        raise SystemExit("Chi should be five")
    if data["superfluidmode"]:
        if data["f2_constant"] != 5.0:
            raise SystemExit("f2 constant should be five")
        if data["evolverType"] != "SuperSplitStep":
            raise SystemExit('In superfluidmode evolverType must be "SuperSplitStep"')
    else:
        if data["evolverType"] not in ("PV2HBSplit23", "PV2HBSplitGeneral"):
            raise SystemExit(
                'evolverType must be "PV2HBSplit23" or "PV2HBSplitGeneral" '
                "for non-superfluid runs"
            )


# dump the data into a .json file
def datatojson():
    with open(data["outputfiletag"] + ".json", "w") as outfile:
        json.dump(data, outfile, indent=4)


# Canonicalize the names for a given set of parameters


def get_kzfilename(tag):
    name = "%s_N%03d_m%08d_h%06d_tkz%06d" % (
        tag,
        data["NX"],
        round(100000 * data["mass0"]),
        round(1000000 * data["H"]),
        round(1.0 / data["dmassdt"]),
    )
    return name


# Canonicalize the names for a given set of parameters


def get_qkzfilename(tag):
    name = "%s_N%03d_m%08d_h%06d_q" % (
        tag,
        data["NX"],
        round(100000 * data["mass0"]),
        round(1000000 * data["H"]),
    )
    return name


def getdefault_filename(tag):
    tag = data["outputfiletag"]
    name = "%s_N%03d_m%08d_h%06d_c%05d" % (
        tag,
        data["NX"],
        round(100000 * data["mass0"]),
        round(1000000 * data["H"]),
        round(100 * data["chi"]),
    )
    return name


# Find the program looking in the environment variable for the path
def find_program(program_name="SuperPions.exe"):
    path = os.environ.get("MODELGEXEPATH")
    if path is None:
        print("Unable to find the path MODELGEXEPATH")
        sys.exit(1)

    abspath = os.path.join(path, program_name)
    if os.path.exists(abspath):
        print("Found the executable {}".format(abspath))
    else:
        print("Unable to find the executable {}".format(abspath))
    return abspath


#########################################################################
# Runs on perlmutter
#########################################################################
def prlmrun(
    time=2,
    debug=False,
    dry_run=True,
    moreopts=["-log_view"],
    seed=None,
    nnodes=1,
    nodeid=False,
    interactive=False,
):
    prgm = find_program()

    # Create a run directory "name"  if does not exist, and cd to it
    dstack.pushd(data["outputfiletag"], mkdir=True)

    # If nodeid is True then append a random 8 digit hex number
    # to the tag labelling the run. This is so that independent runs using the
    # same inputfile, with different seeds, can be run in the same directory
    oldtag = data["outputfiletag"]
    if nodeid:
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
    filenamesh = tag + ".sh"
    filenamestdout = tag + ".stdout"

    fh = open(filenamesh, "w")

    tasks = int(nnodes * 128)
    cpuspertask = int(2 * 128 / (tasks / nnodes))
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
        print("#SBATCH -t {}".format(int(math.ceil(time * 60.0))), file=fh)
        print("#SBATCH -N {}".format(nnodes), file=fh)
        print("#SBATCH --ntasks={}".format(tasks), file=fh)
        print("#SBATCH --cpus-per-task={}".format(cpuspertask), file=fh)

    # Set up the shell environment
    print("", file=fh)
    print("export HDF5_DISABLE_VERSION_CHECK=2", file=fh)
    print("", file=fh)
    print("#run the application:", file=fh)
    print('date  "+%%x %%T" > %s_time.out' % (data["outputfiletag"]), file=fh)
    # Write the command that actually runds the program
    print(
        "srun -n %d --cpu_bind=cores -c %d %s -input %s "
        % (tasks, cpuspertask, prgm, data["outputfiletag"] + ".json"),
        end=" ",
        file=fh,
    )
    # This additional options are  added to the srun command
    for opt in moreopts:
        print(opt, end=" ", file=fh)
    print(file=fh)

    # # Do any post processing of the run
    # programpy = find_program(program_name="x2k.py")
    # print("python {} {}.json".format(programpy,data["outputfiletag"]), file=fh)
    # print(file=fh)

    print('date  "+%%x %%T" >> %s_time.out' % (data["outputfiletag"]), file=fh)
    fh.close()

    # Submit the shell script
    if not dry_run:
        if interactive:
            with open(filenamestdout, "w") as outfile:
                subprocess.run(["sh", filenamesh], stdout=outfile, check=True)
        else:
            subprocess.run(["sbatch", filenamesh])

    # There was a side effect that the outputfiletag got modified
    # This should be undone for transparency
    if nodeid:
        data["outputfiletag"] = oldtag
    # return to the root directory
    dstack.popd()

# runs the actual command current value of data  with mpiexec

########################################################################
# runs the program with current value of data  and mpiexec on local
# mac.
########################################################################


def run(
    program_name="SuperPions.exe",
    moreopts=[],
    dry_run=True,
    time=0,
    seed=None,
    ncpus="2",
    log_view=True,
    mpiexec="mpiexec",
):

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
    opts = [mpiexec, "-n", ncpus, prgm, "-input", tag + ".json"]
    if log_view:
        opts.append("-log_view")
    opts.extend(moreopts)
    print(opts)
    if not dry_run:
        subprocess.run(opts)

    # Go back to the working directory
    dstack.popd()


if __name__ == "__main__":
    print(getdefault_filename("foo"))
