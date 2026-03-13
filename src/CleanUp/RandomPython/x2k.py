import subprocess
import os
import sys
import json
import argparse

# The purpose of the this script is to post-process the output from
# the Model G Code.
#
# Each event has an output file that needs to be fourier transformed.
# and rotated into the vev-plane.
#
# Each run produces many events.
#
# The run has an input file foo.json
#
# Calling:
#
#  x2k.py foo.json
#
# Finds the serial event list and does the jobs. This involves opening
# the outputfile foo.h5 and producing the processed output foo_out.h5
#
#
# usage: x2k.py [-h] name.json [name.json ...]

# Computes the fourier transform of a run or runs

# positional arguments:
#   name.json   name.json or a list of jsonfiles

# optional arguments:
#   -h, --help  show this help message and exit
def find_program(program_name):
    # find the program
    # Assumes that MODELGPATH is set in the environemn
    path = os.path.abspath(os.path.join(os.getenv('MODELGPATH'),program_name))
    print(path)
    return path

def x2k(filename) :
    # Run the x2k program for wallx wally and wallz computing 
    # the fourier transform in the x, y, and z directions
    program = find_program("x2k.exe")
    cmd = program + " " + filename + " wallx"
    print('x2k: processing:', filename)
    result = subprocess.run(cmd, shell=True, capture_output=True)
    cmd = program + " " + filename + " wally"
    result = subprocess.run(cmd, shell=True, capture_output=True)
    cmd = program + " " + filename + " wallz"
    result = subprocess.run(cmd, shell=True, capture_output=True)

def x2k_rotated(filename) :
    # Run the x2k_roated program for wallx wally and wallz
    # This rotates the O4 field in field space
    print('x2k_rotated: processing:', filename)
    program = find_program("x2k_rotated.exe")
    cmd = program + " " + filename + " wallx_k"
    result = subprocess.run(cmd, shell=True, capture_output=True)
    cmd = program + " " + filename + " wally_k"
    result = subprocess.run(cmd, shell=True, capture_output=True)
    cmd = program + " " + filename + " wallz_k"
    result = subprocess.run(cmd, shell=True, capture_output=True)

#Processes foo.json
def process_file(name) :
    with open(name,"r") as fjson: 
        data = json.load(fjson)
        for event in range(0,data["nevents"]) :
            event_name = name.replace(".json", "_{:04d}.h5".format(event))
            print(event_name)
            if os.path.isfile(event_name):
                x2k(event_name)
                x2k_rotated(event_name.replace(".h5","_out.h5"))
            else:
                print("Can't find file ", event_name)

# Main program
if __name__ == '__main__': 
    parser = argparse.ArgumentParser(description="Computes the fourier transform of a run or runs")
    parser.add_argument("files", metavar="name.json",nargs='+', help='name.json or a list of jsonfiles')
    values = parser.parse_args(sys.argv[1:])

    jsonfiles = [k for k in values.files if '.json' in k]
    for file in jsonfiles:
        process_file(file)

