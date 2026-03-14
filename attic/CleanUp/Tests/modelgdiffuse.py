

#!/usr/bin/env python3
import argparse
import sys
import os
import pprint
import math


def main():
    """ Runs the program to make a test of diffusion. The current test just
    initializes a gaussian in x and evolves it for a period of time. 

    The expected solution solution is plotted with modelgdiffuse_plot.py
    """
    args = sys.argv[1:]

    paths = '/Users/derekteaney/common/superfluid/superpaper5/langevinA_061325_superquench/langevinA/ModelG'
    mpiexec = "mpiexec-mpich-clang17"
    os.environ['MODELGPATH'] = paths
    sys.path.append(paths + '/Run')

    import grunner

    prgrm = grunner.find_program('SuperPions.exe')
    pprint.pprint(grunner.data)
    grunner.data["NX"] = 32
    grunner.data["deltat"] = 1./12.
    grunner.data["outputfiletag"] = "modelgdiffuse"
    grunner.data["finaltime"] = 40
    grunner.data["initialization"] = "gaussians"
    grunner.data["writeFrequency"] = 1
    grunner.data["gaussians"] = {"sigmax": 3.,
                                 "sigmay": 3*10e8,
                                 "sigmaz": 3*10e8,
                                 "amplitude": 1,
                                 "theta": 0.0,
                                 "phi": 0.0}
    grunner.data["evolverType"] = "ModelGDiffusionStep"
    grunner.data["ModelGDiffusionStep"] = {"use_implicit_step": False}
    grunner.run(dry_run=False, seed=123, ncpus="2",
                mpiexec=mpiexec, moreopts=sys.argv[1:])


if __name__ == "__main__":
    main()
