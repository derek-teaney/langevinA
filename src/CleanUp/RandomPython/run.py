import os 
import subprocess
import random
import matplotlib as mpl 
import matplotlib.pyplot as plt 
import h5py
from ipywidgets import interact, interactive, fixed, interact_manual
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
random.seed()

data = {
	"H" : 0.01,
	"NX" : 16,
	"NY" : 16,
	"NZ" : 16,
	"LX" : (16-1),
	"LY" : (16-1),
	"LZ" : (16-1),
	"Ndof": 4,
	"finaltime" : 15,
	"intialtime" : 0,
	"deltat"     : 0.01,
	"mass"  : -10,
	"lambda " : 5/4,
	"gamma"   : 1,
	"seed" : 1,
    "filename" : "output"
}

def addoption(opt, value,opts):
    key="-o4_data_" + opt
    sval=str(value)
    opts.extend([key,sval])


def run(moreopts=[]) :
    path =os.path.abspath('')
#   path = os.path.abspath(os.path.dirname(__file__))
    prgm = path + "/ModelA-Beuler.exe"
#    opts = ["/opt/local/lib/petsc/lib/petsc/bin/petscmpiexec","-n", "8",prgm]
    opts = ["mpiexec","-n", "1",prgm]
    data["seed"] = random.randint(1,2000000)
    for datum in data:
        addoption(datum,data[datum],opts) 
    opts.extend(moreopts)
    print(opts)
    subprocess.run(opts)



for h in range(0,1,1) :
    data["H"]=float(h)
    data["lambda"]=5
    data["filename"]="H"+str(h)
    data["finaltime"]=30
    run([])

