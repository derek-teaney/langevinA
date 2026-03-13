import h5py
import sys
import argparse
import os.path 
import numpy as np
import subprocess
import dstack
import re


def update1(filename) :
    """ Updates a generic file with fourier structure, and if does not exist,
    the timeout structure 
    """
    if not h5exists(filename, "timeout"):
        maketimeout(filename) 
    x2kall(filename) 

def update_corrx_0(filename):
    if not h5exists(filename, "wallx") :
        rename_dset(filename, "corrx", "wallx") 
    x2kall(filename)  
    if not h5exists(filename, "timeout") :
        maketimeout(filename, dt=0.72)

def update_corrx_modelA(filename):
    rename_dset(filename, "corrx", "wallx") 
    x2kall(filename)  
    if not h5exists(filename, "timeout") :
        maketimeout(filename, dt=0.8)
    
################################################################################
#
#    
def rename_dset(filename, oldset, newset):
    """Renames a dataset in a file"""
    if not os.path.exists(filename):
        print("Unable to find file: {}".format(filename))
        return

    file = h5py.File(filename, 'r+')
    if oldset in file:
        file[newset] = file[oldset]
        del file[oldset]
    else:
        print("Dataset {} does not exist in file {}".format(filename, oldset))

def x2k(filename, dset):
    """Calls x2k on a file. Using the c++ x2k.exe routine """
    this_dir, this_filename = os.path.split(os.path.abspath(__file__))
    cmd = this_dir + "/x2k.exe"
    subprocess.run([cmd, filename, dset])

def h5exists(filename, dset) :
    """Check if a dataset exists in the file"""
    file = h5py.File(filename, 'r')
    return dset in file

def x2kall(filename) :
    """Calls x2k for the walls"""
    file = h5py.File(filename, 'r+')
    walls = []
    if "wallx" in file:
        walls.append("wallx")
    if "wally" in file:
        walls.append("wally")
    if "wallz" in file:
        walls.append("wallz")
    file.close()
    for wall in walls:
        x2k(filename, wall)
        
def maketimeout(filename, dt=0.72, mass0=None, timeoutname = "timeout") :
    """Makes the timeout data structure

    Arguments:
        filename (str): 
        dt (float): The timestep, default 0.72
        mass0 (float): value of the mass, stripped from the name
        timeoutname (str): name of the timeout array, default="timeout"

    Outputs:
        Writes the timeout array into the file 
    """
    file = h5py.File(filename, 'r+') 
    if timeoutname in file:
        print("maketimeout: timeout data structure already exists")
        return

    phi = file["phi"]

    # Construct the time steps based on dt
    timeout = np.zeros((phi.shape[0],2), dtype=np.float64) 
    timeout[:, 0] = np.arange(0,phi.shape[0]*dt, dt, dtype=np.float64)

    # deduce mass0
    this_dir, this_filename = os.path.split(filename)
    match = re.search(r'_m-(\d+)_', this_filename)
    if match:
        mass0 = - int(match.group(1))/100000
        timeout[:,1] = np.float64(mass0)
    else:
        print("Unable to deduce mass0, using zero") 

    print("Making the timeout data structure with dt = {} and mass = {} in file {}".format(dt,mass0,filename)) 
    
    file.create_dataset(timeoutname, maxshape=(None,2), data=timeout)
    
def _testmaketimeout() :
    fname = "x2ktestdata2_N032_m-0470052_h003000_c00500" + "/" + "x2ktestdata2_N032_m-0470052_h003000_c00500.h5"

    file = h5py.File(fname, 'r+') 
    if "timeouttest" in file:
        del file["timeouttest"] 
    file.close()

    maketimeout(fname, timeoutname="timeouttest")
    
    subprocess.run(["h5dump", "-d", "/timeout", "-d", "/timeouttest", fname]) 

def _testx2kall() :
    fname1 = "x2ktestdata2_N032_m-0470052_h003000_c00500" + "/" + "x2ktestdata2_N032_m-0470052_h003000_c00500.h5"

    file = h5py.File(fname1, 'r+') 
    if "wallx_k" in file:
        del file["wallx_k"] 
        del file["wally_k"] 
        del file["wallz_k"] 
    file.close()

    fname2 = "x2ktestdata2_N032_m-0470052_h003000_c00500" + "/" + "x2ktestdata2_N032_m-0470052_h003000_c00500.h5"


    x2kall(fname1)    
    subprocess.run(["h5diff", fname1, fname2, "/wallx_k"])

def _testupdate():
    fname = "x2ktestdata2_N032_m-0470052_h003000_c00500" + "/" + "x2ktestdata2_N032_m-0470052_h003000_c00500.h5"

    if h5exists(fname, "wallx_k") :
        f = h5py.File(fname, 'r+') 
        del f["wallx_k"]
        del f["wally_k"]
        del f["wallz_k"]
        del f["timeout"]
        f.close()

    update1(fname)

#_testupdate()
if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filename")
    parser.add_argument("-v", "--version", choices=["default", "corrx", "corrx_modelA"], default="default")
    args = parser.parse_args()


    if args.version == "default":
        update1(args.filename)
    elif args.version == "corrx":
        update_corrx_0(args.filename)
    elif args.version == "corrx_modelA":
        update_corrx_modelA(args.filename)
    
    

    

    

    